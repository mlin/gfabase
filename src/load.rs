use clap::Clap;
use genomicsqlite::ConnectionMethods;
use json::object;
use log::{debug, error, info, warn};
use rusqlite::{params, OpenFlags, OptionalExtension, Statement, Transaction, NO_PARAMS};
use std::collections::HashMap;

use crate::invalid_gfa;
use crate::util;
use crate::util::Result;

#[derive(Clap)]
pub struct Opts {
    /// uncompressed GFA file/stream (use - for stdin)
    pub gfa: String,

    /// destination gfab filename
    pub gfab: String,

    /// add genomic range index from rGFA input
    #[clap(long)]
    pub rgfa: bool,

    /// Zstandard compression level (-7 to 22)
    #[clap(long, default_value = "6")]
    pub compress: i8,
}

pub fn main(opts: &Opts) -> Result<()> {
    let prefix = "";

    // formulate GenomicSQLite configuration JSON
    let mut db = new_db(&opts.gfab, opts.compress, 1024)?;

    {
        // open transaction & apply schema
        let txn = db.transaction()?;
        create_tables(&txn, prefix, opts.rgfa)?;

        // add a temp table to hold segment metadata, which we'll copy into the main db file after
        // writing all the segment sequences; this ensures the metadata is stored ~contiguously
        // instead of interspersed among the (typically much larger) sequence data.
        txn.execute_batch(
            "CREATE TABLE temp.segment_meta_hold(
                segment_id INTEGER PRIMARY KEY, name TEXT,
                sequence_length INTEGER, tags_json TEXT
            )",
        )?;

        // intake GFA records
        info!("processing GFA1 records...");
        insert_gfa1(&opts.gfa, &txn, prefix, opts.rgfa)?;
        info!("writing segment metadata...");
        txn.execute_batch(&format!(
            "INSERT INTO {}gfa1_segment_meta(segment_id, name, sequence_length, tags_json)
             SELECT segment_id, name, sequence_length, tags_json FROM temp.segment_meta_hold",
            prefix
        ))?;
        info!("insertions complete");

        // indexing
        create_indexes(&txn, prefix, opts.rgfa)?;

        // done
        info!("flushing {} ...", &opts.gfab);
        txn.commit()?;
    }

    summary(&db)?;
    db.close().map_err(|(_, e)| e)?;
    info!("🗹 done");
    Ok(())
}

pub fn new_db(filename: &str, compress: i8, page_cache_MiB: i32) -> Result<rusqlite::Connection> {
    // formulate GenomicSQLite configuration JSON
    let dbopts = match object! {
        unsafe_load: true,
        inner_page_KiB: 64,
        outer_page_KiB: 2,
        zstd_level: compress,
        page_cache_MiB: page_cache_MiB
    } {
        json::JsonValue::Object(o) => o,
        _ => {
            assert!(false);
            json::object::Object::new()
        }
    };

    // create db
    util::delete_existing_file(filename)?;
    let db = genomicsqlite::open(
        filename,
        OpenFlags::SQLITE_OPEN_CREATE
            | OpenFlags::SQLITE_OPEN_READ_WRITE
            | OpenFlags::SQLITE_OPEN_NO_MUTEX,
        &dbopts,
    )?;
    db.execute_batch("PRAGMA foreign_keys = OFF")?;
    Ok(db)
}

pub fn create_tables(db: &rusqlite::Connection, prefix: &str, rgfa: bool) -> Result<()> {
    let mut gfa_sql = include_str!("schema/GFA1.sql").to_string();
    gfa_sql = util::simple_placeholder(&gfa_sql, "prefix", prefix);
    db.execute_batch(&gfa_sql)?;

    if rgfa {
        let mut rgfa_sql = include_str!("schema/rGFA1.sql").to_string();
        rgfa_sql = util::simple_placeholder(&rgfa_sql, "prefix", prefix);
        db.execute_batch(&rgfa_sql)?;
        info!("created [r]GFA1 tables");
    } else {
        info!("created GFA1 tables");
    }

    Ok(())
}

pub fn create_indexes(db: &rusqlite::Connection, prefix: &str, rgfa: bool) -> Result<()> {
    info!("indexing:");
    create_indexes_exec(db, &include_str!("schema/GFA1.index.sql"), prefix)?;

    if rgfa {
        create_indexes_exec(db, &include_str!("schema/rGFA1.index.sql"), prefix)?;

        let gri_sql = db.create_genomic_range_index_sql(
            &format!("{}gfa1_reference", prefix),
            "rid",
            "position",
            "position+length",
        )?;
        info!("\tGenomic Range Indexing ...");
        db.execute_batch(&gri_sql)?;
    }

    info!("\tANALYZE ...");
    db.execute_batch("PRAGMA analysis_limit = 1000; ANALYZE main")?;

    Ok(())
}

fn create_indexes_exec(db: &rusqlite::Connection, sql: &str, prefix: &str) -> Result<()> {
    for index_spec in sql.split(";") {
        let mut index_sql = index_spec.trim().to_string();
        if !index_sql.is_empty() {
            index_sql = util::simple_placeholder(&index_sql, "prefix", prefix);
            info!("\t{} ...", index_sql.splitn(2, " ON").next().unwrap());
            db.execute_batch(&index_sql)?;
        }
    }
    Ok(())
}

fn insert_gfa1(filename: &str, txn: &Transaction, prefix: &str, rgfa: bool) -> Result<()> {
    // prepared statements
    let mut stmt_insert_segment_meta =
        txn.prepare("INSERT INTO temp.segment_meta_hold(segment_id,name,sequence_length,tags_json) VALUES(?,?,?,?)")?;
    let mut stmt_insert_segment_sequence = txn.prepare(&format!(
        "INSERT INTO {}gfa1_segment_sequence(segment_id,sequence_twobit) VALUES(?,nucleotides_twobit(?))",
        prefix
    ))?;
    let mut stmt_insert_link = txn.prepare(&format!(
        "INSERT INTO {}gfa1_link(from_segment,from_reverse,to_segment,to_reverse,cigar,tags_json) VALUES(?,?,?,?,?,?)",
        prefix
    ))?;
    let mut stmt_insert_rgfa = None;
    if rgfa {
        stmt_insert_rgfa = Some(txn.prepare(&format!(
            "INSERT INTO {}gfa1_reference(segment_id,rid,position,length,rank) VALUES(?,?,?,?,?)",
            prefix
        ))?)
    }
    let mut stmt_insert_path = txn.prepare(&format!(
        "INSERT INTO {}gfa1_path(path_id,name,tags_json) VALUES(?,?,?)",
        prefix
    ))?;
    let mut stmt_insert_path_element = txn.prepare(&format!(
        "INSERT INTO {}gfa1_path_element(path_id,ordinal,segment_id,reverse,cigar_vs_next) VALUES(?,?,?,?,?)",
        prefix
    ))?;

    let mut segments_by_name = HashMap::new();

    // closure to process one record
    let dispatch = |tsv: &Vec<&str>| -> Result<()> {
        match tsv[0] {
            "S" => insert_gfa1_segment(
                tsv,
                txn,
                &mut stmt_insert_segment_meta,
                &mut stmt_insert_segment_sequence,
                stmt_insert_rgfa.as_mut(),
                &mut segments_by_name,
            ),
            "L" => insert_gfa1_link(tsv, &mut stmt_insert_link, &segments_by_name),
            "P" => insert_gfa1_path(
                tsv,
                txn,
                &mut stmt_insert_path,
                &mut stmt_insert_path_element,
                &segments_by_name,
            ),
            _ => invalid_gfa!("GFA record type {} not yet supported", tsv[0]),
        }
    };

    // iterate tsv records
    util::iter_tsv_no_comments(dispatch, filename, Some('#' as u8))
}

fn insert_gfa1_segment(
    tsv: &Vec<&str>,
    txn: &Transaction,
    stmt_meta: &mut Statement,
    stmt_sequence: &mut Statement,
    maybe_stmt_rgfa: Option<&mut Statement>,
    segments_by_name: &mut HashMap<String, i64>,
) -> Result<()> {
    assert_eq!(tsv[0], "S");
    if tsv.len() < 2 {
        invalid_gfa!("malformed S line");
    }

    let rowid = name_to_id(tsv[1]);
    let name = if rowid.is_some() { None } else { Some(tsv[1]) };
    let maybe_sequence = if tsv.len() > 2 && tsv[2] != "*" {
        Some(tsv[2])
    } else {
        None
    };
    let mut tags_json = prepare_tags_json(tsv, 3)?;

    let ln_tag = tags_json.remove("LN:i").map(|j| j.as_i64()).flatten();
    let sequence_len = match (maybe_sequence, ln_tag) {
        (Some(seq), Some(lni)) if lni != (seq.len() as i64) => {
            invalid_gfa!(
                "segment with inconsistent sequence length and LN tag: {}",
                tsv[1]
            )
        }
        (Some(seq), _) => Some(seq.len() as i64),
        (None, Some(lni)) => Some(lni),
        (None, None) => None,
    };

    let rgfa_tags = if maybe_stmt_rgfa.is_some() {
        let sn = tags_json
            .remove("SN:Z")
            .map(|j| j.as_str().map(|s| String::from(s)))
            .flatten();
        let so = tags_json.remove("SO:i").map(|j| j.as_i64()).flatten();
        let sr = tags_json.remove("SR:i").map(|j| j.as_i64()).flatten();
        if sn.is_none() || so.is_none() || sr.is_none() || sequence_len.is_none() {
            invalid_gfa!(
                "segment missing required rGFA tags (SN:Z SO:i SR:i LN:i): {}",
                tsv[1]
            )
        }
        Some((sn.unwrap(), so.unwrap(), sr.unwrap()))
    } else {
        None
    };

    let tags_json_text = tags_json.dump();
    stmt_meta.execute(params![rowid, name, sequence_len, tags_json_text])?;
    let rowid_actual = txn.last_insert_rowid();

    if let Some(nm) = name {
        segments_by_name.insert(String::from(nm), rowid_actual);
    }
    if let Some(seq) = maybe_sequence {
        stmt_sequence.execute(params![rowid_actual, seq])?;
    }
    if let Some(stmt_rgfa) = maybe_stmt_rgfa {
        let (sn, so, sr) = rgfa_tags.unwrap();
        stmt_rgfa.execute(params!(rowid_actual, sn, so, sequence_len.unwrap(), sr))?;
    }

    Ok(())
}

fn insert_gfa1_link(
    tsv: &Vec<&str>,
    stmt: &mut Statement,
    segments_by_name: &HashMap<String, i64>,
) -> Result<()> {
    assert_eq!(tsv[0], "L");
    if tsv.len() < 5 {
        invalid_gfa!("malformed L line: {}", tsv.join("\t"));
    }

    let (from_segment, from_reverse) = segment_and_orientation(tsv[1], tsv[2], segments_by_name)?;
    let (to_segment, to_reverse) = segment_and_orientation(tsv[3], tsv[4], segments_by_name)?;
    let cigar = if tsv.len() > 5 && tsv[5] != "*" {
        Some(tsv[5])
    } else {
        None
    };
    let tags_json = prepare_tags_json(tsv, 6)?;
    let tags_json_text = tags_json.dump();
    stmt.execute(params![
        from_segment,
        from_reverse,
        to_segment,
        to_reverse,
        cigar,
        if tags_json_text.trim() != "{}" {
            Some(tags_json_text)
        } else {
            None
        }
    ])?;
    Ok(())
}

fn insert_gfa1_path(
    tsv: &Vec<&str>,
    txn: &Transaction,
    stmt_path: &mut Statement,
    stmt_ele: &mut Statement,
    segments_by_name: &HashMap<String, i64>,
) -> Result<()> {
    assert_eq!(tsv[0], "P");
    if tsv.len() < 3 {
        invalid_gfa!("malformed P line: {}", tsv.join("\t"));
    }

    let rowid = name_to_id(tsv[1]);
    let name = if rowid.is_some() { None } else { Some(tsv[1]) };

    let maybe_cigars: Option<Vec<&str>> = if tsv.len() > 3 && tsv[3] != "*" {
        Some(tsv[3].split(',').collect())
    } else {
        None
    };
    let tags_json = prepare_tags_json(tsv, 4)?;
    let tags_json_text = tags_json.dump();

    stmt_path.execute(params![
        rowid,
        name,
        if tags_json_text.trim() != "{}" {
            Some(tags_json_text)
        } else {
            None
        }
    ])?;
    let rowid_actual = txn.last_insert_rowid();

    let segs: Vec<&str> = tsv[2].split(',').collect();
    if let Some(cigars) = &maybe_cigars {
        if cigars.len() + 1 != segs.len() {
            invalid_gfa!("incorrect overlaps: {}", tsv.join("\t"));
        }
    }
    for ord in 0..segs.len() {
        let ele = segs[ord];
        if ele.len() < 2 {
            invalid_gfa!("malformed path: {}", tsv[2]);
        }
        let (segment_id, reverse) = segment_and_orientation(
            &ele[..(ele.len() - 1)],
            &ele[(ele.len() - 1)..ele.len()],
            segments_by_name,
        )?;
        let cigar = match (&maybe_cigars, ord) {
            (Some(cigars), i) if i < segs.len() - 1 => Some(cigars[i]),
            _ => None,
        };
        stmt_ele.execute(params![
            rowid_actual,
            ord as i64,
            segment_id,
            reverse,
            cigar
        ])?;
    }
    Ok(())
}

fn segment_and_orientation(
    segment: &str,
    orientation: &str,
    segments_by_name: &HashMap<String, i64>,
) -> Result<(i64, i64)> {
    let segment_id = if let Some(id) = name_to_id(segment) {
        id
    } else if let Some(idr) = segments_by_name.get(segment) {
        *idr
    } else {
        invalid_gfa!("unknown segment in link/path: {}", segment)
    };
    let reverse = match orientation {
        "+" => 0,
        "-" => 1,
        _ => {
            invalid_gfa!("malformed segment orientation: {}", orientation)
        }
    };
    Ok((segment_id, reverse))
}

pub fn name_to_id(name: &str) -> Option<i64> {
    let namelen = name.len();
    match name.parse() {
        Ok(i) => Some(i),
        Err(_) if namelen > 1 => name[1..namelen].parse().ok(),
        Err(_) => None,
    }
}

fn prepare_tags_json(tsv: &Vec<&str>, offset: usize) -> Result<json::object::Object> {
    let mut ans = json::object::Object::new();
    if tsv.len() > offset {
        for cursor in offset..tsv.len() {
            let fields: Vec<&str> = tsv[cursor].splitn(3, ':').collect();
            if fields.len() != 3 || fields[0].is_empty() || fields[1].is_empty() {
                invalid_gfa!("malformed tag: {}", tsv[cursor]);
            };
            let v = match fields[1] {
                "A" | "Z" | "H" => json::JsonValue::from(fields[2]),
                "i" => {
                    let iv: i64 = fields[2].parse().or_else(|_| {
                        invalid_gfa!("malformed tag integer: {}", tsv[cursor]);
                    })?;
                    json::JsonValue::from(iv)
                }
                "f" => {
                    let fv: f64 = fields[2].parse().or_else(|_| {
                        invalid_gfa!("malformed tag float: {}", tsv[cursor]);
                    })?;
                    json::JsonValue::from(fv)
                }
                // TODO: B & J
                _ => {
                    invalid_gfa!("tag type not yet supported: {}", tsv[cursor]);
                }
            };
            ans.insert(&format!("{}:{}", fields[0], fields[1]), v)
        }
    }
    Ok(ans)
}

pub fn summary(db: &rusqlite::Connection) -> Result<()> {
    info!("tables & row counts:");
    let mut stmt_tables = db.prepare("SELECT name FROM sqlite_master WHERE type='table'")?;
    let mut tables = stmt_tables.query(NO_PARAMS)?;
    while let Some(row) = tables.next()? {
        let table: String = row.get(0)?;
        let ct: i64 = db.query_row(
            &format!("SELECT count(1) FROM {}", table),
            NO_PARAMS,
            |ctr| ctr.get(0),
        )?;
        info!("\t{}\t{}", table, ct);
    }
    if let Some(e) = db
        .query_row("PRAGMA foreign_key_check", NO_PARAMS, |row| {
            Ok(util::Error::InvalidGfab {
                message: String::from("foreign key integrity violation"),
                table: row.get(0)?,
                rowid: row.get(1)?,
            })
        })
        .optional()?
    {
        return Err(e);
    }
    Ok(())
}
