use clap::Clap;
use genomicsqlite::ConnectionMethods;
use log::{debug, info, log_enabled, warn};
use rusqlite::{params, OpenFlags, OptionalExtension, NO_PARAMS};
use std::collections::{BinaryHeap, HashMap};
use std::{cmp, io};

use crate::util::Result;
use crate::{bad_command, connectivity, load, util, view};

#[derive(Clap)]
pub struct Opts {
    /// gfab filename
    pub gfab: String,

    /// output filename [omit or - for standard output, implying --view]
    #[clap(short, default_value = "-")]
    pub outfile: String,

    /// desired segments/paths/ranges
    #[clap(name = "SEGMENT")]
    pub segments: Vec<String>,

    /// SEGMENTs are Path names to get all segments of
    #[clap(long)]
    pub path: bool,

    /// SEGMENTs are reference sequence ranges like chr7:1,234-5,678 to locate in segment mappings
    #[clap(long)]
    pub range: bool,

    /// Treat SEGMENTs as text names even if they look like integer IDs
    #[clap(long)]
    pub always_names: bool,

    /// Expand from specified segments to complete (undirected) connected component(s)
    #[clap(long)]
    pub connected: bool,

    /// Expand from specified segments to subgraph connected by crossing fewer than R cutpoints
    #[clap(long, name = "R", default_value = "0")]
    pub cutpoints: u64,

    /// Modifies --cutpoints to count only cutpoints with segment sequence at least L nucleotides
    /// (implies --cutpoints >= 1)
    #[clap(long, name = "L", default_value = "0")]
    pub cutpoints_nt: u64,

    /// Write .gfa instead of .gfab to outfile
    #[clap(long)]
    pub view: bool,

    /// Launch Bandage on .gfa outfile (implies --view; create temporary file if outfile unspecified)
    #[clap(long)]
    pub bandage: bool,

    /// For each segment with reference mappings, set gr:Z tag with one guessed range summarizing the mappings (implies --view)
    #[clap(long)]
    pub guess_ranges: bool,

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
// TODO: add help headings e.g. #[clap(long, help_heading = Some("SEGMENT"))]
// pending release of fix for https://github.com/clap-rs/clap/issues/2279

pub fn main(opts: &Opts) -> Result<()> {
    if opts.segments.len() == 0 {
        bad_command!("specify one or more desired subgraph segments on the command line");
    }
    util::check_gfab_filename_schema(&opts.gfab)?;
    if opts.view || opts.bandage || opts.guess_ranges || opts.outfile == "-" {
        sub_gfa(opts)
    } else {
        sub_gfab(opts)
    }
}

fn sub_gfab(opts: &Opts) -> Result<()> {
    assert_ne!(opts.outfile, "-");
    if !opts.view && !opts.outfile.ends_with(".gfab") {
        warn!("output filename should end in .gfab");
    }
    if opts.view && opts.outfile != "-" && !opts.outfile.ends_with(".gfa") {
        warn!("output filename should end in .gfa");
    }

    // create output database
    let mut db = load::new_db(&opts.outfile, opts.compress, 1024)?;

    // attach input database
    let mut dbopts_in = json::object::Object::new();
    dbopts_in.insert("immutable", json::JsonValue::from(true));
    let attach_sql = db.genomicsqlite_attach_sql(&opts.gfab, "input", &dbopts_in)?;
    db.execute_batch(&attach_sql)?;
    let gfab_version = util::check_gfab_schema(&db, "input.")?;
    debug!("gfabase v{} created {}", gfab_version, opts.gfab);
    util::check_gfab_version(&gfab_version)?;

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

    let mut maybe_guesser = if opts.guess_ranges {
        Some(view::SegmentRangeGuesser::new(
            &txn,
            "WHERE segment_id IN temp.sub_segments",
        )?)
    } else {
        None
    };

    if opts.outfile == "-" && !opts.bandage && atty::is(atty::Stream::Stdout) {
        // interactive mode: pipe into less -S
        view::less(|less_in| sub_gfa_write(&txn, &mut maybe_guesser, !opts.no_sequences, less_in))?
    } else {
        let mut output_gfa = String::from(&opts.outfile);
        if opts.bandage && output_gfa == "-" {
            output_gfa = view::bandage_temp_filename()?
        }

        {
            let mut writer_box = view::writer(&output_gfa)?;
            sub_gfa_write(
                &txn,
                &mut maybe_guesser,
                !opts.no_sequences,
                &mut *writer_box,
            )?;
        }

        if opts.bandage {
            if let Some(ref mut guesser) = maybe_guesser {
                guesser.write_bandage_csv(&output_gfa)?
            }
            view::bandage(&output_gfa)?
        }
    }
    std::mem::drop(maybe_guesser);

    txn.rollback()?;
    Ok(())
}

fn sub_gfa_write(
    db: &rusqlite::Connection,
    maybe_guesser: &mut Option<view::SegmentRangeGuesser>,
    sequences: bool,
    out: &mut dyn io::Write,
) -> Result<()> {
    let mut tag_editor = |segment_id: i64, tags: &mut json::JsonValue| -> Result<()> {
        if let Some(ref mut guesser) = maybe_guesser {
            if let Some(gr) = guesser.get(segment_id)? {
                tags.insert("gr:Z", gr).unwrap()
            }
        }
        Ok(())
    };

    view::write_header(db, out)?;
    view::write_segments(
        db,
        "WHERE segment_id IN temp.sub_segments",
        sequences,
        &mut tag_editor,
        out,
    )?;
    view::write_links(
        db,
        "WHERE from_segment IN temp.sub_segments AND to_segment IN temp.sub_segments",
        out,
    )?;
    view::write_paths(&db, "WHERE path_id IN temp.sub_paths", out)
}

// populate temp.sub_segments with the segment IDs of the desired subgraph
fn compute_subgraph(db: &rusqlite::Connection, opts: &Opts, input_schema: &str) -> Result<()> {
    let mut check_start_segments = false;
    db.execute_batch("CREATE TABLE temp.start_segments(segment_id INTEGER PRIMARY KEY)")?;
    {
        let mut insert_segment = if !opts.range {
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
            if opts.range {
                if insert_segment.execute(params![segment])? < 1 {
                    bad_command!("no segments found overlapping {}", segment);
                }
            } else if !opts.always_names && load::name_to_id(segment).is_some() {
                let id = load::name_to_id(segment).unwrap();
                if !opts.path {
                    insert_segment.execute(params![id]).map(|_| ())?;
                    check_start_segments = true;
                } else if insert_path.execute(params![id])? < 1 {
                    bad_command!("unknown path {}", id);
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

    let mut cutpoints = opts.cutpoints;
    if opts.cutpoints_nt > 0 {
        cutpoints = cmp::max(1, cutpoints);
    }

    if opts.connected {
        if cutpoints > 0 {
            warn!("--connected overriding --cutpoints");
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
            warn!("input .gfab lacks connectivity index; discovering connected component(s) on-the-fly (requires SQLite >= 3.34.0)...");
            let connected_sql = include_str!("query/connected.sql").to_string()
                + "\nALTER TABLE temp.connected_segments RENAME TO sub_segments";
            db.execute_batch(&connected_sql)?;
        }
    } else if cutpoints > 0 {
        if !connectivity::has_index(db, input_schema)? {
            panic!("input .gfab lacks connectivity index required for --cutpoints")
        }
        expand_to_cutpoints(db, input_schema, cutpoints as i64, opts.cutpoints_nt as i64)?
    } else {
        // sub_segments = start_segments
        db.execute_batch("ALTER TABLE temp.start_segments RENAME TO sub_segments")?;
    }

    Ok(())
}

// Compute all segments connected to a set of start segments without traversing a cutpoint segment
// (generally: traverse up to R-1 cutpoint segments of at least L bp each, defaults R=1 L=0).
// "Traversing a cutpoint segment" means arriving on one side and departing the other; so when we
// "stop" at a cutpoint segment, we still follow other edges touching the side at which we arrived.
// Otherwise, treat non-cutpoint segments as unsided and all edges as undirected.
//
//  IN: segment IDs in temp.start_segments
// OUT: segment IDs in temp.sub_segments
fn expand_to_cutpoints(
    db: &rusqlite::Connection,
    input_schema: &str,
    radius: i64,
    cutpoint_length: i64,
) -> Result<()> {
    // query for segment's neighbors, also indicating which sides of segment & neighbor are linked
    let neighbors_sql = format!(
        "   SELECT (NOT from_reverse) AS rhs, to_segment AS neighbor, to_reverse AS neighbor_rhs
            FROM {s}gfa1_link WHERE from_segment = ?1
         UNION
            SELECT to_reverse AS rhs, from_segment AS neighbor, (NOT from_reverse) AS neighbor_rhs
            FROM {s}gfa1_link WHERE to_segment = ?1",
        s = input_schema
    );
    let mut neighbors_query = db.prepare(&neighbors_sql)?;

    let mut is_cutpoint_query = db.prepare(&format!(
        "SELECT is_cutpoint FROM {}gfa1_connectivity WHERE segment_id = ?",
        input_schema
    ))?;
    let mut sequence_length_query = db.prepare(&format!(
        "SELECT sequence_length FROM {}gfa1_segment_meta WHERE segment_id = ?",
        input_schema
    ))?;
    let mut is_cutpoint_memo: HashMap<i64, bool> = HashMap::new();

    // max-heap on (remaining_radius, segment_id, rhs)
    // Scheduled visit to (left/right-hand side of) segment, with remaining_radius "energy" left.
    // Prioritizing visits energy-first promotes locality
    let mut queue: BinaryHeap<(i64, i64, bool)> = BinaryHeap::new();
    // (segment, cutpoint_rhs) -> max remaining_radius of previous visits
    let mut visited: HashMap<(i64, bool), i64> = HashMap::new();

    // seed the queue with starting segments (both sides, in case a starting segment is itself a
    // cutpoint)
    {
        let mut start_segments = db.prepare("SELECT segment_id FROM temp.start_segments")?;
        let mut start_segments_cursor = start_segments.query(NO_PARAMS)?;
        while let Some(row) = start_segments_cursor.next()? {
            let segment_id: i64 = row.get(0)?;
            queue.push((radius - 1, segment_id, false));
            queue.push((radius - 1, segment_id, true));
        }
    }

    // take from queue until empty
    while let Some((remaining_radius, segment, rhs)) = queue.pop() {
        let is_cutpoint = if let Some(&ans) = is_cutpoint_memo.get(&segment) {
            ans
        } else {
            let mut ans = is_cutpoint_query
                .query_row(params![segment], |row| row.get(0))
                .optional()?
                .unwrap_or(0 as i64)
                != 0;
            if ans && cutpoint_length != 0 {
                let len: i64 =
                    sequence_length_query.query_row(params![segment], |row| row.get(0))?;
                ans = len >= cutpoint_length;
            }
            is_cutpoint_memo.insert(segment, ans);
            ans
        };

        // skip if we've already visited segment with at least remaining_radius
        if visited
            .get(&(segment, is_cutpoint && rhs))
            .map_or(false, |&visited_radius| visited_radius >= remaining_radius)
        {
            continue;
        }
        // record this visit
        visited.insert((segment, is_cutpoint && rhs), remaining_radius);

        // schedule visits to segment's neighbors:
        //   those linked to same side of segment, always
        //   those linked to opposite side -- unless we're at a cutpoint & out of energy
        let mut neighbors_cursor = neighbors_query.query(params![segment])?;
        while let Some(row) = neighbors_cursor.next()? {
            let link_rhs: bool = row.get(0)?;
            let neighbor: i64 = row.get(1)?;
            let neighbor_rhs: bool = row.get(2)?;
            if rhs == link_rhs || !is_cutpoint || remaining_radius > 0 {
                let cost = if is_cutpoint && rhs != link_rhs { 1 } else { 0 };
                queue.push((remaining_radius - cost, neighbor, neighbor_rhs))
            }
        }
    }

    // fill temp.sub_segments with all visited segment IDs
    db.execute_batch("CREATE TABLE temp.sub_segments(segment_id INTEGER PRIMARY KEY)")?;
    let mut insert_sub =
        db.prepare("INSERT OR IGNORE INTO temp.sub_segments(segment_id) VALUES(?)")?;
    for (segment, _) in visited.keys() {
        insert_sub.execute(params![segment])?;
    }

    Ok(())
}
