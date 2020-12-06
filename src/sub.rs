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

    // load the desired segment IDs into a temp table
    {
        let txn = db.transaction()?;
        txn.execute_batch(
            "CREATE TABLE temp.desired_segments(segment_id INTEGER PRIMARY KEY);
                           CREATE TABLE temp.unknown_segments(segment_id INTEGER PRIMARY KEY);",
        )?;
        {
            let mut insert_segment =
                txn.prepare("INSERT INTO temp.desired_segments(segment_id) VALUES(?)")?;
            for segment in &opts.segments {
                match load::segment_name_to_id(segment) {
                    None => bad_command!("don't know what to do with {}", segment), // FIXME
                    Some(segment_id) => insert_segment.execute(params![segment_id]).map(|_| ())?,
                }
            }
        }
        txn.commit()?
    }

    // attach input database
    let mut dbopts_in = json::object::Object::new();
    dbopts_in.insert("immutable", json::JsonValue::from(true));
    let attach_sql = db.genomicsqlite_attach_sql(&opts.gfab, "input", &dbopts_in)?;
    debug!("{}", attach_sql);
    db.execute_batch(&attach_sql)?;

    // check existence of all desired segment IDs
    db.execute_batch(&format!(
        "INSERT INTO temp.unknown_segments
         SELECT temp.desired_segments.segment_id AS segment_id
         FROM temp.desired_segments LEFT JOIN input.{}gfa1_segment_meta USING (segment_id)
         WHERE input.{}gfa1_segment_meta.segment_id IS NULL",
        prefix, prefix
    ))?;
    let n_unknown: i64 = db.query_row(
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

    let rgfa = view::is_rgfa(&db, "input.", prefix)?;
    {
        let txn = db.transaction()?;
        // TODO: if rgfa, copy _gri_refseq
        load::create_tables(&txn, prefix, rgfa)?;

        info!("copying segment sequences...");
        txn.execute_batch(&format!(
            "INSERT INTO {}gfa1_segment_sequence(segment_id,sequence_twobit)
             SELECT segment_id, sequence_twobit FROM input.{}gfa1_segment_sequence
                WHERE segment_id IN temp.desired_segments
                ORDER BY segment_id",
            prefix, prefix
        ))?;

        info!("copying segment metadata...");
        txn.execute_batch(&format!(
            "INSERT INTO {}gfa1_segment_meta(segment_id,name,tags_json)
             SELECT segment_id, name, tags_json FROM input.{}gfa1_segment_meta
                WHERE segment_id IN temp.desired_segments
                ORDER BY segment_id",
            prefix, prefix
        ))?;

        if rgfa {
            // TODO: copy _gri_refseq
            info!("copying rGFA coordinates...");
            txn.execute_batch(&format!(
                "INSERT INTO {}gfa1_reference(segment_id,rid,position,length,rank)
                 SELECT segment_id, rid, position, length, rank FROM input.{}gfa1_reference
                    WHERE segment_id IN temp.desired_segments
                    ORDER BY segment_id",
                prefix, prefix
            ))?;
        }

        info!("copying relevant links...");
        txn.execute_batch(&format!(
            "INSERT INTO {}gfa1_link(from_segment, from_reverse, to_segment, to_reverse, cigar, tags_json)
             SELECT from_segment, from_reverse, to_segment, to_reverse, cigar, tags_json FROM input.{}gfa1_link
                WHERE from_segment IN temp.desired_segments
                AND to_segment IN temp.desired_segments
                ORDER BY from_segment, to_segment",
            prefix, prefix
        ))?;

        load::create_indexes(&txn, prefix, rgfa)?;

        info!("flushing {} ...", &opts.outfile);
        txn.commit()?
    }

    load::summary(&db)?;
    db.close().map_err(|(_, e)| e)?;
    info!("🗹 done");
    Ok(())
}
