use clap::Clap;
use genomicsqlite::ConnectionMethods;
use log::{debug, info, warn};
use rusqlite::{params, NO_PARAMS};

use crate::util::Result;
use crate::{bad_command, load, util, view};

#[derive(Clap)]
pub struct Opts {
    /// gfab filename
    pub gfab: String,

    /// output filename (required for default .gfab output)
    #[clap(default_value = "")]
    pub outfile: String,

    /// expand desired segment(s) to full connected component(s)
    #[clap(long)]
    pub connected: bool,

    /// <SEGMENT>s are reference ranges like chr7:1,234-5,678 to locate in rGFA
    #[clap(long)]
    pub reference: bool,

    /// desired segment(s)
    #[clap(name = "SEGMENT")]
    pub segments: Vec<String>,
    /*
    // skip .gfab output & render .gfa to outfile (or stdout if omitted)
    #[clap(long)]
    pub view: bool,
    */
}

pub fn main(opts: &Opts) -> Result<()> {
    let sub_segment_count: i64;

    if opts.outfile == "" || opts.outfile == "-" {
        bad_command!("missing required outfile argument");
    }
    if opts.segments.len() == 0 {
        bad_command!("specify one or more desired segments on the command line");
    }

    // create output database
    let mut db = load::new_db(&opts.outfile, 6, 1024)?;

    // attach input database
    let mut dbopts_in = json::object::Object::new();
    dbopts_in.insert("immutable", json::JsonValue::from(true));
    let attach_sql = db.genomicsqlite_attach_sql(&opts.gfab, "input", &dbopts_in)?;
    db.execute_batch(&attach_sql)?;

    {
        let txn = db.transaction()?;

        // populate temp ID table with the user-specified segments
        txn.execute_batch("CREATE TABLE temp.start_segments(segment_id INTEGER PRIMARY KEY)")?;
        {
            let mut insert_segment = if !opts.reference {
                txn.prepare("INSERT OR REPLACE INTO temp.start_segments(segment_id) VALUES(?)")?
            } else {
                // GRI query
                txn.prepare(
                    "INSERT OR REPLACE INTO temp.start_segments(segment_id)
                     SELECT segment_id FROM gfa1_segment_mapping
                        WHERE _rowid_ in genomic_range_rowids(
                            'gfa1_segment_mapping',
                            parse_genomic_range_sequence(?1),
                            parse_genomic_range_begin(?1),
                            parse_genomic_range_end(?1))",
                )?
            };
            for segment in &opts.segments {
                if opts.reference {
                    if insert_segment.execute(params![segment])? < 1 {
                        bad_command!("no segments found overlapping {}", segment);
                    }
                } else if let Some(segment_id) = load::name_to_id(segment) {
                    insert_segment.execute(params![segment_id]).map(|_| ())?;
                } else {
                    // FIXME try to find by segment name
                    bad_command!("don't know what to do with {}", segment);
                }
            }
        }

        // check existence of all the specified segments
        txn.execute_batch(
            "CREATE TABLE temp.unknown_segments(segment_id INTEGER PRIMARY KEY);
             INSERT INTO temp.unknown_segments
                 SELECT temp.start_segments.segment_id AS segment_id
                 FROM temp.start_segments LEFT JOIN input.gfa1_segment_meta USING (segment_id)
                 WHERE input.gfa1_segment_meta.segment_id IS NULL",
        )?;
        {
            let mut missing_examples = Vec::new();
            let mut missing_examples_stmt =
                txn.prepare("SELECT segment_id FROM temp.unknown_segments LIMIT 10")?;
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

        if opts.connected {
            // sub_segments = connected components including start_segments
            let connected_sql = include_str!("query/connected_old.sql").to_string()
                + "\nALTER TABLE temp.connected_segments RENAME TO sub_segments";
            debug!("computing connected component(s)...");
            txn.execute_batch(&connected_sql)?;
        } else {
            // sub_segments = start_segments
            txn.execute_batch("ALTER TABLE temp.start_segments RENAME TO sub_segments")?;
        }

        sub_segment_count =
            txn.query_row("SELECT count(1) FROM temp.sub_segments", NO_PARAMS, |row| {
                row.get(0)
            })?;

        load::create_tables(&txn)?;

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
