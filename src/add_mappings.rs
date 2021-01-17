use clap::Clap;
use log::{debug, info, warn};
use rusqlite::{params, OpenFlags, OptionalExtension, NO_PARAMS};

use crate::bad_command;
use crate::load;
use crate::util;
use crate::util::Result;

#[derive(Clap)]
pub struct Opts {
    /// Assembly .gfab filename (to modify in-place; copy first if needed)
    pub gfab: String,
    /// Uncompressed .paf filename [omit or - for standard input]
    #[clap(default_value = "-")]
    pub mappings: String,

    /// First delete all existing segment mappings from .gfab
    #[clap(long)]
    pub replace: bool,

    /// Ignore mappings with lower quality score
    #[clap(long, default_value = "0")]
    pub quality: u64,

    /// Ignore mappings with shorter alignment block length
    #[clap(long, default_value = "0")]
    pub length: u64,

    /// Treat mappings' query names as text segment names even if they look like integer IDs
    #[clap(long)]
    pub always_names: bool,

    /// Ignore (instead of error on) mappings whose query doesn't match a segment in the .gfab
    #[clap(long)]
    pub ignore_unknown: bool,
}

pub fn main(opts: &Opts) -> Result<()> {
    util::check_gfab_filename_schema(&opts.gfab)?;
    if opts.mappings == "-" && atty::is(atty::Stream::Stdin) {
        bad_command!("pipe in .paf data or supply uncompressed filename")
    }

    // formulate GenomicSQLite configuration JSON
    let mut dbopts = json::object::Object::new();
    dbopts.insert("unsafe_load", json::JsonValue::from(true));

    // open db
    let mut db = genomicsqlite::open(
        &opts.gfab,
        OpenFlags::SQLITE_OPEN_READ_WRITE | OpenFlags::SQLITE_OPEN_NO_MUTEX,
        &dbopts,
    )?;
    let gfab_version = util::check_gfab_schema(&db, "")?;
    debug!("gfabase v{} created {}", gfab_version, opts.gfab);
    util::check_gfab_version(&gfab_version)?;

    {
        // open transaction & apply schema
        let txn = db.transaction()?;
        insert_paf(&txn, &opts)?;
        debug!("flushing {} ...", &opts.gfab);
        txn.commit()?;
    }

    Ok(())
}

macro_rules! invalid_paf {
    ($($arg:tt)*) => (util::Error::InvalidPaf(format!($($arg)*)))
}

pub fn insert_paf(db: &rusqlite::Connection, opts: &Opts) -> Result<()> {
    // create temp table
    db.execute_batch(
        "CREATE TABLE temp.segment_mapping_hold(
            segment_id INTEGER NOT NULL,
            refseq_name TEXT NOT NULL COLLATE UINT,
            refseq_begin INTEGER NOT NULL,
            refseq_end INTEGER NOT NULL,
            tags_json TEXT
        );",
    )?;

    let mut segment_id_check =
        db.prepare("SELECT segment_id FROM gfa1_segment_meta WHERE segment_id = ?")?;
    let mut segment_name_to_id =
        db.prepare("SELECT segment_id FROM gfa1_segment_meta WHERE name = ?")?;
    let mut insert_mapping = db.prepare("INSERT INTO temp.segment_mapping_hold(segment_id, refseq_name, refseq_begin, refseq_end, tags_json) VALUES(?,?,?,?,?)")?;

    // iterate tsv records
    let mut insert_count = 0;
    let mut all_count = 0;
    let mut unknown_count = 0;
    let insert_paf1 = |tsv: &Vec<&str>| -> Result<()> {
        all_count += 1;
        if tsv.len() < 12 {
            return Err(invalid_paf!("malformed PAF line: {}", tsv.join("\n")));
        }
        // check if mapping passes filters
        if opts.length > 0 {
            let aln_length: u64 = tsv[10].parse().map_err(|_| {
                invalid_paf!("malformed alignment block length: {}", tsv.join("\n"))
            })?;
            if aln_length < opts.length {
                return Ok(());
            }
        }
        if opts.quality > 0 {
            let map_q: u64 = tsv[11]
                .parse()
                .map_err(|_| invalid_paf!("malformed mapQ: {}", tsv.join("\n")))?;
            if map_q < opts.quality {
                return Ok(());
            }
        }
        // look up segment ID
        let mut maybe_segment_id = None;
        if !opts.always_names {
            if let Some(id) = load::name_to_id(tsv[0]) {
                maybe_segment_id = Some(id)
            }
        }
        if let Some(id) = maybe_segment_id {
            maybe_segment_id = segment_id_check
                .query_row(params![id], |row| row.get(0))
                .optional()?
        } else {
            maybe_segment_id = segment_name_to_id
                .query_row(params![tsv[0]], |row| row.get(0))
                .optional()?
        }
        if maybe_segment_id.is_none() {
            if opts.ignore_unknown {
                unknown_count += 1;
                return Ok(());
            }
            return Err(invalid_paf!("query name isn't a known segment: {}", tsv[0]));
        }
        let segment_id = maybe_segment_id.unwrap();
        let segment_begin: u64 = tsv[2]
            .parse()
            .map_err(|_| invalid_paf!("malformed query start: {}", tsv.join("\n")))?;
        let segment_end: u64 = tsv[3]
            .parse()
            .map_err(|_| invalid_paf!("malformed query end: {}", tsv.join("\n")))?;
        // parse target range
        let target_name = tsv[5];
        // TODO: handle GAF path if target_name starts with '>' or '<'
        let target_begin: u64 = tsv[7]
            .parse()
            .map_err(|_| invalid_paf!("malformed target start: {}", tsv.join("\n")))?;
        let target_end: u64 = tsv[8]
            .parse()
            .map_err(|_| invalid_paf!("malformed target end: {}", tsv.join("\n")))?;
        if target_begin > target_end {
            return Err(invalid_paf!("target begin > end: {}", tsv.join("\n")))?;
        }
        // prepare tags
        let mut tags_json = json::object::Object::new();
        tags_json.insert("sb:i", json::JsonValue::from(segment_begin));
        tags_json.insert("se:i", json::JsonValue::from(segment_end));
        tags_json.insert("so:Z", json::JsonValue::from(tsv[4]));
        // insert into temp table
        insert_mapping.execute(params![
            segment_id,
            target_name,
            target_begin as i64,
            target_end as i64,
            tags_json.dump()
        ])?;
        insert_count += 1;
        Ok(())
    };
    util::iter_tsv_no_comments(insert_paf1, &opts.mappings, Some('#' as u8))?;
    if unknown_count > 0 {
        warn!(
            "ignored {} mappings with unknown query names",
            unknown_count
        )
    }

    // delete existing mappings if desired
    debug!("buffered {} of {} mappings", insert_count, all_count);
    if opts.replace {
        let deleted = db.execute("DELETE FROM gfa1_segment_mapping", NO_PARAMS)?;
        if deleted > 0 {
            warn!("deleted {} existing mappings", deleted)
        }
    }
    // sort temp table into gfab
    debug!("sorting mappings...");
    db.execute_batch(
        "INSERT INTO gfa1_segment_mapping(segment_id, refseq_name, refseq_begin, refseq_end, tags_json)
            SELECT segment_id, refseq_name, refseq_begin, refseq_end, tags_json
            FROM temp.segment_mapping_hold NOT INDEXED
            ORDER BY refseq_name COLLATE UINT, refseq_begin, refseq_end",
    )?;
    info!("inserted {} of {} mappings", insert_count, all_count);
    return Ok(());
}
