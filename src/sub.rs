use clap::Clap;
use genomicsqlite::ConnectionMethods;
use log::{debug, info, warn};
use rusqlite::{params, OpenFlags, OptionalExtension, NO_PARAMS};

use crate::util::Result;
use crate::{bad_command, load, util, view};

#[derive(Clap)]
pub struct Opts {
    /// gfab filename
    pub gfab: String,

    /// output filename
    pub outfile: String,

    /// follow links from specified segment(s) to get full connected component(s)
    #[clap(long)]
    pub connected: bool,

    /// <SEGMENT>s are actually paths to get all segments of
    #[clap(long)]
    pub path: bool,

    /// <SEGMENT>s are reference ranges like chr7:1,234-5,678 to locate in segment mappings
    #[clap(long)]
    pub reference: bool,

    /// desired segment(s) or path(s)
    #[clap(name = "SEGMENT")]
    pub segments: Vec<String>,

    // instead of .gfab, output .gfa to outfile (- for stdout)
    #[clap(long)]
    pub view: bool,
}

pub fn main(opts: &Opts) -> Result<()> {
    if opts.segments.len() == 0 {
        bad_command!("specify one or more desired subgraph segments on the command line");
    }
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
    let mut db = load::new_db(&opts.outfile, 6, 1024)?;

    // attach input database
    let mut dbopts_in = json::object::Object::new();
    dbopts_in.insert("immutable", json::JsonValue::from(true));
    let attach_sql = db.genomicsqlite_attach_sql(&opts.gfab, "input", &dbopts_in)?;
    db.execute_batch(&attach_sql)?;

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
            txn.execute_batch(include_str!("query/sub.sql"))?;
        }

        load::create_indexes(&txn)?;

        debug!("flushing {} ...", &opts.outfile);
        txn.commit()?
    }

    load::summary(&db)?;
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
    let db = genomicsqlite::open(
        &opts.gfab,
        OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
        &dbopts,
    )?;

    // open output writer
    let mut writer_box = view::writer(Some(&opts.outfile))?;
    let out = &mut *writer_box;

    compute_subgraph(&db, opts, "")?;
    db.execute_batch(
        // desired paths are those without a segment not in temp.sub_segments
        "CREATE TABLE temp.sub_paths(path_id INTEGER PRIMARY KEY);
         INSERT INTO temp.sub_paths(path_id)
            SELECT path_id FROM gfa1_path
            WHERE path_id NOT IN
                (SELECT DISTINCT path_id FROM gfa1_path_element
                 WHERE segment_id NOT IN temp.sub_segments)",
    )?;

    view::write_segments(&db, "WHERE segment_id IN temp.sub_segments", out)?;
    view::write_links(
        &db,
        "WHERE from_segment IN temp.sub_segments AND to_segment IN temp.sub_segments",
        out,
    )?;
    view::write_paths(&db, "WHERE path_id IN temp.sub_paths", out)?;
    out.flush()?;

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
            db.prepare(
                "INSERT OR REPLACE INTO temp.start_segments(segment_id)
                 SELECT segment_id FROM gfa1_segment_mapping
                    WHERE _rowid_ in genomic_range_rowids(
                        'gfa1_segment_mapping',
                        parse_genomic_range_sequence(?1),
                        parse_genomic_range_begin(?1),
                        parse_genomic_range_end(?1))",
            )?
        };
        let mut find_segment_by_name =
            db.prepare("SELECT segment_id FROM gfa1_segment_meta WHERE name=?")?;
        let mut find_path_by_name = db.prepare("SELECT path_id FROM gfa1_path WHERE name=?")?;
        let mut insert_path = db.prepare(
            "INSERT OR REPLACE INTO temp.start_segments(segment_id)
             SELECT segment_id FROM gfa1_path_element WHERE path_id=?",
        )?;
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

    if opts.connected {
        // sub_segments = connected components including start_segments
        let connected_sql = include_str!("query/connected_old.sql").to_string()
            + "\nALTER TABLE temp.connected_segments RENAME TO sub_segments";
        debug!("computing connected component(s)...");
        db.execute_batch(&connected_sql)?;
    } else {
        // sub_segments = start_segments
        db.execute_batch("ALTER TABLE temp.start_segments RENAME TO sub_segments")?;
    }

    Ok(())
}
