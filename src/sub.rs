use clap::Clap;
use genomicsqlite::ConnectionMethods;
use json::object;
use log::{debug, error, info, warn};
use rusqlite::{params, OpenFlags, Statement, NO_PARAMS};

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
    let prefix = "";

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
    debug!("{}", attach_sql);
    db.execute_batch(&attach_sql)?;

    {
        let txn = db.transaction()?;

        // populate temp ID table with the user-specified segments
        txn.execute_batch("CREATE TABLE temp.start_segments(segment_id INTEGER PRIMARY KEY)")?;
        {
            let mut insert_segment = if !opts.reference {
                txn.prepare("INSERT INTO temp.start_segments(segment_id) VALUES(?)")?
            } else {
                txn.prepare(&format!(
                    "INSERT INTO temp.start_segments(segment_id)
                     SELECT segment_id FROM {}gfa1_reference
                        WHERE _rowid_ in genomic_range_rowids(
                            '{}gfa1_reference',
                            parse_genomic_range_sequence(?1),
                            parse_genomic_range_begin(?1),
                            parse_genomic_range_end(?1))",
                    prefix, prefix
                ))?
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
        txn.execute_batch(&format!(
            "CREATE TABLE temp.unknown_segments(segment_id INTEGER PRIMARY KEY);
             INSERT INTO temp.unknown_segments
                 SELECT temp.start_segments.segment_id AS segment_id
                 FROM temp.start_segments LEFT JOIN input.{}gfa1_segment_meta USING (segment_id)
                 WHERE input.{}gfa1_segment_meta.segment_id IS NULL",
            prefix, prefix
        ))?;
        let n_unknown: i64 = txn.query_row(
            "SELECT count(1) FROM temp.unknown_segments",
            NO_PARAMS,
            |row| row.get(0),
        )?;
        if n_unknown > 0 {
            // TODO: error! log up to 10 of the missing segment IDs
            bad_command!(
                "{} of {} desired segment IDs aren't present in {}",
                n_unknown,
                opts.segments.len(),
                &opts.gfab
            );
        }

        if opts.connected {
            // sub_segments = connected components including start_segments
            let mut connected_sql = include_str!("query/connected_old.sql").to_string();
            connected_sql = util::simple_placeholder(&connected_sql, "prefix", prefix)
                + "\nALTER TABLE temp.connected_segments RENAME TO sub_segments";
            info!("computing connected component(s)...");
            txn.execute_batch(&connected_sql)?;
        } else {
            // sub_segments = start_segments
            txn.execute_batch("ALTER TABLE temp.start_segments RENAME TO sub_segments")?;
        }

        let rgfa = view::is_rgfa(&txn, "input.", prefix)?;
        load::create_tables(&txn, prefix, rgfa)?;

        info!("copying segments & links...");
        let mut sub_sql = include_str!("query/sub.sql").to_string();
        sub_sql = util::simple_placeholder(&sub_sql, "prefix", prefix);
        txn.execute_batch(&sub_sql)?;

        if rgfa {
            // TODO: copy _gri_refseq
            info!("copying rGFA coordinates...");
            txn.execute_batch(&format!(
                "INSERT INTO {}gfa1_reference(segment_id,rid,position,length,rank)
                 SELECT segment_id, rid, position, length, rank FROM input.{}gfa1_reference
                    WHERE segment_id IN temp.sub_segments
                    ORDER BY segment_id",
                prefix, prefix
            ))?;
        }

        load::create_indexes(&txn, prefix, rgfa)?;

        info!("flushing {} ...", &opts.outfile);
        txn.commit()?
    }

    load::summary(&db)?;
    db.close().map_err(|(_, e)| e)?;
    info!("🗹 done");
    Ok(())
}
