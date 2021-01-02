use clap::Clap;
use genomicsqlite::ConnectionMethods;
use log::{debug, info, log_enabled, warn};
use rusqlite::{params, OpenFlags, OptionalExtension, NO_PARAMS};
use std::collections::{BTreeMap, BinaryHeap, HashMap};

use crate::util::Result;
use crate::{bad_command, connectivity, load, util, view};

#[derive(Clap)]
pub struct Opts {
    /// gfab filename
    pub gfab: String,

    /// output filename
    pub outfile: String,

    /// desired segment(s) or path(s)
    #[clap(name = "SEGMENT")]
    pub segments: Vec<String>,

    /// <SEGMENT>s are actually paths to get all segments of
    #[clap(long)]
    pub path: bool,

    /// <SEGMENT>s are reference ranges like chr7:1,234-5,678 to locate in segment mappings
    #[clap(long)]
    pub reference: bool,

    /// Follow links from specified segment(s) to get full connected component(s)
    #[clap(long)]
    pub connected: bool,

    /// Follow links from specified segment(s) out to given radius
    #[clap(long)]
    pub radius: Option<u32>,

    /// Follow links from specified segment(s) out to given radius, each link weighted by segment
    /// sequence length
    #[clap(long)]
    pub radius_bp: Option<u32>,

    /// Follow links from specified segment(s), but don't cross through graph cutpoints (combine
    /// with --radius or --radius-bp to stop at N'th cutpoint)
    #[clap(long)]
    pub to_cutpoints: bool,

    /// Write .gfa instead of .gfab to outfile (- for stdout)
    #[clap(long)]
    pub view: bool,

    /// Omit segment sequences
    #[clap(long)]
    pub no_sequences: bool,

    /// Omit index of subgraph.gfab connectivity
    #[clap(long)]
    pub no_connectivity: bool,

    /// compression level (-5 to 22) for output .gfab
    #[clap(long, default_value = "6")]
    pub compress: i8,
}

pub fn main(opts: &Opts) -> Result<()> {
    if opts.segments.len() == 0 {
        bad_command!("specify one or more desired subgraph segments on the command line");
    }
    util::check_gfabase_filename_schema(&opts.gfab)?;
    if !opts.view {
        sub_gfab(opts)
    } else {
        sub_gfa(opts)
    }
}

fn sub_gfab(opts: &Opts) -> Result<()> {
    if opts.outfile == "" || opts.outfile == "-" {
        bad_command!("output .gfab filename required; cannot output .gfab to stdout")
    }
    if !opts.outfile.ends_with(".gfab") {
        warn!("output filename should end in .gfab");
    }

    // create output database
    let mut db = load::new_db(&opts.outfile, opts.compress, 1024)?;

    // attach input database
    let mut dbopts_in = json::object::Object::new();
    dbopts_in.insert("immutable", json::JsonValue::from(true));
    let attach_sql = db.genomicsqlite_attach_sql(&opts.gfab, "input", &dbopts_in)?;
    db.execute_batch(&attach_sql)?;
    util::check_gfabase_schema(&db, "input.")?;

    let sub_segment_count: i64;
    {
        let txn = db.transaction()?;

        compute_subgraph(&txn, opts, "input.")?;
        load::create_tables(&txn)?;
        sub_segment_count =
            txn.query_row("SELECT count(1) FROM temp.sub_segments", NO_PARAMS, |row| {
                row.get(0)
            })?;

        if sub_segment_count == 0 {
            warn!("no segments matched the command-line criteria")
        } else {
            info!(
                "copying {} segments (and links & paths touching only those segments)",
                sub_segment_count
            );
            if !opts.no_sequences {
                txn.execute_batch(
                    "INSERT INTO gfa1_segment_sequence(segment_id, sequence_twobit)
                     SELECT segment_id, sequence_twobit FROM input.gfa1_segment_sequence
                     WHERE segment_id IN temp.sub_segments",
                )?;
            }
            txn.execute_batch(include_str!("query/sub.sql"))?;
        }

        load::create_indexes(&txn, !opts.no_connectivity)?;

        debug!("flushing {} ...", &opts.outfile);
        txn.commit()?
    }

    if log_enabled!(log::Level::Debug) {
        load::summary(&db)?;
    }
    db.close().map_err(|(_, e)| e)?;
    if sub_segment_count == 0 {
        return Err(util::Error::EmptyGfab);
    }
    info!("🗹 done");
    Ok(())
}

fn sub_gfa(opts: &Opts) -> Result<()> {
    let mut dbopts = json::object::Object::new();
    dbopts.insert("immutable", json::JsonValue::from(true));

    // open db
    let mut db = genomicsqlite::open(
        &opts.gfab,
        OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
        &dbopts,
    )?;

    // open output writer
    let mut writer_box = view::writer(Some(&opts.outfile))?;
    let out = &mut *writer_box;

    let txn = db.transaction()?;
    compute_subgraph(&txn, opts, "")?;
    txn.execute_batch(
        // desired paths are those without a segment not in temp.sub_segments
        "CREATE TABLE temp.sub_paths(path_id INTEGER PRIMARY KEY);
         INSERT INTO temp.sub_paths(path_id)
            SELECT path_id FROM gfa1_path
            WHERE path_id NOT IN
                (SELECT DISTINCT path_id FROM gfa1_path_element
                 WHERE segment_id NOT IN temp.sub_segments)",
    )?;

    view::write_header(&txn, out)?;
    view::write_segments(
        &txn,
        "WHERE segment_id IN temp.sub_segments",
        !opts.no_sequences,
        out,
    )?;
    view::write_links(
        &txn,
        "WHERE from_segment IN temp.sub_segments AND to_segment IN temp.sub_segments",
        out,
    )?;
    view::write_paths(&txn, "WHERE path_id IN temp.sub_paths", out)?;
    out.flush()?;

    txn.rollback()?;
    Ok(())
}

// populate temp.sub_segments with the segment IDs of the desired subgraph
fn compute_subgraph(db: &rusqlite::Connection, opts: &Opts, input_schema: &str) -> Result<()> {
    let mut check_start_segments = false;
    db.execute_batch("CREATE TABLE temp.start_segments(segment_id INTEGER PRIMARY KEY)")?;
    {
        let mut insert_segment = if !opts.reference {
            db.prepare("INSERT OR REPLACE INTO temp.start_segments(segment_id) VALUES(?)")?
        } else {
            // GRI query
            db.prepare(&format!(
                "INSERT OR REPLACE INTO temp.start_segments(segment_id)
                     SELECT segment_id FROM {}gfa1_segment_mapping
                        WHERE _rowid_ in genomic_range_rowids(
                            '{}gfa1_segment_mapping',
                            parse_genomic_range_sequence(?1),
                            parse_genomic_range_begin(?1),
                            parse_genomic_range_end(?1))",
                input_schema, input_schema
            ))?
        };
        let mut find_segment_by_name = db.prepare(&format!(
            "SELECT segment_id FROM {}gfa1_segment_meta WHERE name=?",
            input_schema
        ))?;
        let mut find_path_by_name = db.prepare(&format!(
            "SELECT path_id FROM {}gfa1_path WHERE name=?",
            input_schema
        ))?;
        let mut insert_path = db.prepare(&format!(
            "INSERT OR REPLACE INTO temp.start_segments(segment_id)
             SELECT segment_id FROM {}gfa1_path_element WHERE path_id=?",
            input_schema
        ))?;
        for segment in &opts.segments {
            if opts.reference {
                if insert_segment.execute(params![segment])? < 1 {
                    bad_command!("no segments found overlapping {}", segment);
                }
            } else if let Some(segment_id) = load::name_to_id(segment) {
                if !opts.path {
                    insert_segment.execute(params![segment_id]).map(|_| ())?;
                    check_start_segments = true;
                } else if insert_path.execute(params![segment_id])? < 1 {
                    bad_command!("unknown path {}", segment_id);
                }
            } else if !opts.path {
                let maybe_segment_id: Option<i64> = find_segment_by_name
                    .query_row(params![segment], |row| row.get(0))
                    .optional()?;
                if let Some(segment_id) = maybe_segment_id {
                    insert_segment.execute(params![segment_id])?;
                } else {
                    bad_command!("unknown segment {}", segment)
                }
            } else {
                let maybe_path_id: Option<i64> = find_path_by_name
                    .query_row(params![segment], |row| row.get(0))
                    .optional()?;
                let mut inserted = 0;
                if let Some(path_id) = maybe_path_id {
                    inserted = insert_path.execute(params![path_id])?;
                }
                if inserted == 0 {
                    bad_command!("unknown path {}", segment)
                }
            }
        }
    }

    if check_start_segments {
        // check existence of all the specified segments
        db.execute_batch(&format!(
            "CREATE TABLE temp.unknown_segments(segment_id INTEGER PRIMARY KEY);
             INSERT INTO temp.unknown_segments
                SELECT temp.start_segments.segment_id AS segment_id
                FROM temp.start_segments LEFT JOIN {}gfa1_segment_meta USING (segment_id)
                WHERE {}gfa1_segment_meta.segment_id IS NULL",
            input_schema, input_schema
        ))?;
        {
            let mut missing_examples = Vec::new();
            let mut missing_examples_stmt =
                db.prepare("SELECT segment_id FROM temp.unknown_segments LIMIT 10")?;
            let mut missing_examples_cursor = missing_examples_stmt.query(NO_PARAMS)?;
            while let Some(row) = missing_examples_cursor.next()? {
                let missing_segment_id: i64 = row.get(0)?;
                missing_examples.push(missing_segment_id.to_string())
            }
            if missing_examples.len() > 0 {
                bad_command!(
                    "desired segment IDs aren't present in {} such as: {}",
                    &opts.gfab,
                    missing_examples.join(" ")
                );
            }
        }
    }

    let mut opts_radius = opts.radius;
    if opts.to_cutpoints {
        if !connectivity::has_index(db, input_schema)? {
            panic!("input .gfab lacks connectivity index required for --to-cutpoints")
        }
        if opts_radius.is_none() && opts.radius_bp.is_none() {
            opts_radius = Some(1)
        }
    }

    if opts.connected {
        if opts.radius.is_some() || opts.radius_bp.is_some() {
            warn!("--connected overriding --radius[-bp]");
        }
        if opts.to_cutpoints {
            warn!("--connected overriding --to-cutpoints");
        }
        // sub_segments = connected components including start_segments
        if connectivity::has_index(db, input_schema)? {
            let connected_sql = format!(
                "CREATE TABLE temp.sub_segments(segment_id INTEGER PRIMARY KEY);
                 INSERT INTO temp.sub_segments(segment_id)
                        SELECT segment_id FROM {s}gfa1_connectivity WHERE component_id IN
                            (SELECT DISTINCT component_id FROM {s}gfa1_connectivity
                             WHERE segment_id IN temp.start_segments)
                    UNION ALL
                        SELECT segment_id FROM
                            temp.start_segments LEFT JOIN {s}gfa1_connectivity USING(segment_id)
                        WHERE component_id IS NULL",
                s = input_schema
            );
            db.execute_batch(&connected_sql)?;
        } else {
            warn!("input .gfab lacks connectivity index; discovering connected component(s) on-the-fly...");
            let connected_sql = include_str!("query/connected.sql").to_string()
                + "\nALTER TABLE temp.connected_segments RENAME TO sub_segments";
            db.execute_batch(&connected_sql)?;
        }
    } else if let Some(radius) = opts_radius {
        if opts.radius_bp.is_some() {
            warn!("--radius overriding --radius-bp")
        }
        connected_radius(db, input_schema, radius as i64, false, opts.to_cutpoints)?;
    } else if let Some(radius_bp) = opts.radius_bp {
        connected_radius(db, input_schema, radius_bp as i64, true, opts.to_cutpoints)?;
    } else {
        // sub_segments = start_segments
        db.execute_batch("ALTER TABLE temp.start_segments RENAME TO sub_segments")?;
    }

    Ok(())
}

// Compute all segments reachable by following up to `radius` links (in any direction) from a set
// of start segments. Optionally, each link can be weighted by the sequence length of the segment
// reached by following it. Then `radius` is roughly a bp distance (but It's Complicated, because
// this doesn't consider the direction of the links).
//
// Also, we might only weight links reaching graph cutpoints. "Stopping" at a cutpoint means not
// "crossing through" the cutpoint segment -- we still follow the cutpoint's other links going in
// the same direction (relative to the cutpoint) as the link through which we arrived.
//
//  IN: segment IDs in temp.start_segments
// OUT: segment IDs in temp.connected_segments
fn connected_radius(
    db: &rusqlite::Connection,
    input_schema: &str,
    radius: i64,
    sequence_length_weighted: bool,
    to_cutpoints: bool,
) -> Result<()> {
    // max-heap on (remaining_radius, 0-segment_id, reverse_link)
    // reverse_link tracks the direction of the link we followed to get to segment (true = segment
    // is the 'from' end of the link)
    let mut queue: BinaryHeap<(i64, i64, bool)> = BinaryHeap::new();
    // (segment_id,reverse_link) -> max remaining_radius of completed visits
    // reverse_link may be true only for cutpoints
    let mut visited: BTreeMap<(i64, bool), i64> = BTreeMap::new();
    let mut cutpoints: HashMap<i64, bool> = HashMap::new();

    // fill queue with starting segments & initial radius
    {
        let mut start_segments = db.prepare("SELECT segment_id FROM temp.start_segments")?;
        let mut start_segments_cursor = start_segments.query(NO_PARAMS)?;
        while let Some(row) = start_segments_cursor.next()? {
            let segment_id: i64 = row.get(0)?;
            queue.push((radius, 0 - segment_id, false));
        }
    }

    {
        // queries for incoming & outgoing links of segment, and "weight" of each link (either
        // sequence length or 1)
        let outgoing_links_sql = if sequence_length_weighted {
            format!(
                "SELECT to_segment, sequence_length
                 FROM {s}gfa1_link LEFT JOIN {s}gfa1_segment_meta ON to_segment = segment_id
                 WHERE from_segment = ?",
                s = input_schema
            )
        } else {
            format!(
                "SELECT to_segment, 1 FROM {}gfa1_link WHERE from_segment = ?",
                input_schema
            )
        };
        let mut outgoing_links_query = db.prepare(&outgoing_links_sql)?;
        let incoming_links_sql = if sequence_length_weighted {
            format!(
                "SELECT from_segment, sequence_length
                 FROM {s}gfa1_link LEFT JOIN {s}gfa1_segment_meta ON from_segment = segment_id
                 WHERE to_segment = ?",
                s = input_schema
            )
        } else {
            format!(
                "SELECT from_segment, 1 FROM {}gfa1_link WHERE to_segment = ?",
                input_schema
            )
        };
        let mut incoming_links_query = db.prepare(&incoming_links_sql)?;
        let mut is_cutpoint_query = db.prepare(&format!(
            "SELECT is_cutpoint FROM {}gfa1_connectivity WHERE segment_id = ?",
            input_schema
        ))?;

        // take from queue until empty
        while let Some((remaining_radius, neg_segment_id, reverse_link)) = queue.pop() {
            let segment = 0 - neg_segment_id;
            let is_cutpoint = to_cutpoints
                && if let Some(&known) = cutpoints.get(&segment) {
                    known
                } else {
                    let ans = is_cutpoint_query
                        .query_row(params![segment], |row| row.get(0))
                        .optional()?
                        .unwrap_or(0 as i64)
                        != 0;
                    cutpoints.insert(segment, ans);
                    ans
                };

            // skip if we've already visited segment with at least remaining_radius
            if visited
                .get(&(segment, is_cutpoint && reverse_link))
                .map_or(false, |&visited_radius| visited_radius >= remaining_radius)
            {
                continue;
            }
            // record this visit
            visited.insert((segment, is_cutpoint && reverse_link), remaining_radius);

            if remaining_radius > 0 {
                // if we still have radius, enqueue the segment's neighbors with appropriate deduction
                let mut links_cursor = outgoing_links_query.query(params![segment])?;
                while let Some(row) = links_cursor.next()? {
                    let neighbor: i64 = row.get(0)?;
                    let weight = if !to_cutpoints || is_cutpoint {
                        row.get(1)?
                    } else {
                        0 as i64
                    };
                    queue.push((remaining_radius - weight, 0 - neighbor, false));
                }
                links_cursor = incoming_links_query.query(params![segment])?;
                while let Some(row) = links_cursor.next()? {
                    let neighbor: i64 = row.get(0)?;
                    let weight = if !to_cutpoints || is_cutpoint {
                        row.get(1)?
                    } else {
                        0 as i64
                    };
                    queue.push((remaining_radius - weight, 0 - neighbor, true));
                }
            } else if is_cutpoint {
                // follow only the links with the same direction according to reverse_link
                let mut links_cursor = if reverse_link {
                    outgoing_links_query.query(params![segment])?
                } else {
                    incoming_links_query.query(params![segment])?
                };
                while let Some(row) = links_cursor.next()? {
                    let neighbor: i64 = row.get(0)?;
                    let weight = if !to_cutpoints || is_cutpoint {
                        row.get(1)?
                    } else {
                        0 as i64
                    };
                    queue.push((remaining_radius - weight, 0 - neighbor, !reverse_link));
                }
            }
        }
    }

    // fill temp.sub_segments with all visited segment IDs
    {
        db.execute_batch("CREATE TABLE temp.sub_segments(segment_id INTEGER PRIMARY KEY)")?;
        let mut insert_sub = db.prepare("INSERT INTO temp.sub_segments(segment_id) VALUES(?)")?;
        for (segment, _) in visited.keys() {
            insert_sub.execute(params![segment])?;
        }
    }

    Ok(())
}
