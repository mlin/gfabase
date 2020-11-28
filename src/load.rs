use clap::Clap;
use genomicsqlite::ConnectionMethods;
use log::{debug, error, info, warn};
use rusqlite::{params, OpenFlags, Statement, Transaction, NO_PARAMS};
use std::collections::HashMap;
use std::path::Path;

use crate::invalid_gfa;
use crate::util;
use crate::util::Result;

#[derive(Clap)]
pub struct Opts {
    /// destination gfab filename
    pub gfab: String,
    /// input GFA file (omit for standard input)
    pub gfa: Option<String>,

    /// add genomic range index from rGFA input
    #[clap(long)]
    pub rgfa: bool,
}

pub fn main(opts: &Opts) -> Result<()> {
    let prefix = "";

    // formulate GenomicSQLite configuration JSON
    let mut dbopts = json::object::Object::new();
    dbopts.insert("unsafe_load", json::JsonValue::from(true));
    dbopts.insert("inner_page_KiB", json::JsonValue::from(64));
    dbopts.insert("outer_page_KiB", json::JsonValue::from(2));

    // delete existing file, if any
    {
        let p = Path::new(&opts.gfab);
        if p.is_file() {
            warn!("delete existing file {}", p.to_str().unwrap());
            std::fs::remove_file(p)?
        }
    }
    // create db
    let mut db = genomicsqlite::open(
        &opts.gfab,
        OpenFlags::SQLITE_OPEN_CREATE | OpenFlags::SQLITE_OPEN_READ_WRITE,
        &dbopts,
    )?;

    // open transaction & apply schema
    let txn = db.transaction()?;
    let mut gfa_sql = include_str!("schema/GFA1.sql").to_string();
    gfa_sql = util::simple_replace_all(&gfa_sql, "prefix", prefix);
    txn.execute_batch(&gfa_sql)?;
    info!("created GFA1 schema in {}", &opts.gfab);

    if opts.rgfa {
        let mut rgfa_sql = include_str!("schema/rGFA1.sql").to_string();
        rgfa_sql = util::simple_replace_all(&rgfa_sql, "prefix", prefix);
        txn.execute_batch(&rgfa_sql)?;
        info!("created rGFA schema too");
    }

    // intake GFA records
    info!("processing GFA1 records...");
    insert_gfa1(
        opts.gfa.as_ref().unwrap_or(&String::from("")),
        &txn,
        prefix,
        opts.rgfa,
    )?;
    info!("insertions complete; indexing...");

    // create indexes
    create_indexes(&include_str!("schema/GFA1.index.sql"), prefix, &txn)?;

    if opts.rgfa {
        create_indexes(&include_str!("schema/rGFA1.index.sql"), prefix, &txn)?;

        let gri_sql = txn.create_genomic_range_index_sql(
            &format!("{}gfa1_reference", prefix),
            "rid",
            "position",
            "position+length",
        )?;
        info!("\tGenomic Range Indexing ...");
        txn.execute_batch(&gri_sql)?;
    }

    // report stats
    info!("tables & row counts:");
    {
        let mut stmt_tables = txn.prepare("SELECT name FROM sqlite_master WHERE type='table'")?;
        let mut tables = stmt_tables.query(NO_PARAMS)?;
        while let Some(row) = tables.next()? {
            let table: String = row.get(0)?;
            let ct: i64 = txn.query_row(
                &format!("SELECT count(1) FROM {}", table),
                NO_PARAMS,
                |ctr| ctr.get(0),
            )?;
            info!("\t{}\t{}", table, ct);
        }
    }

    // done
    info!("flushing {}", &opts.gfab);
    txn.commit()?;
    info!("\t🗹");
    Ok(())
}

fn create_indexes(sql: &str, prefix: &str, txn: &Transaction) -> Result<()> {
    for index_spec in sql.split(";") {
        let mut index_sql = index_spec.trim().to_string();
        if !index_sql.is_empty() {
            index_sql = util::simple_replace_all(&index_sql, "prefix", prefix);
            info!("\t{} ...", index_sql.splitn(2, " ON").next().unwrap());
            txn.execute_batch(&index_sql)?;
        }
    }
    Ok(())
}

fn insert_gfa1(filename: &str, txn: &Transaction, prefix: &str, rgfa: bool) -> Result<()> {
    // prepared statements
    let mut stmt_insert_segment = txn.prepare(&format!(
        "INSERT INTO {}gfa1_segment(_rowid_,name,tags_json,sequence_twobit) VALUES(?,?,?,nucleotides_twobit(?))",
        prefix
    ))?;
    let mut stmt_insert_link = txn.prepare(&format!(
        "INSERT INTO {}gfa1_link(from_segment,from_reverse,to_segment,to_reverse,cigar,tags_json) VALUES(?,?,?,?,?,?)",
        prefix
    ))?;
    let mut stmt_insert_rgfa = None;
    if rgfa {
        stmt_insert_rgfa = Some(txn.prepare(&format!(
            "INSERT INTO {}gfa1_reference(segment,rid,position,length,rank) VALUES(?,?,?,?,?)",
            prefix
        ))?)
    }

    let mut segments_by_name = HashMap::new();

    // closure to process one record
    let dispatch = |tsv: &Vec<&str>| -> Result<()> {
        match tsv[0] {
            "S" => insert_gfa1_segment(
                tsv,
                txn,
                &mut stmt_insert_segment,
                stmt_insert_rgfa.as_mut(),
                &mut segments_by_name,
            ),
            "L" => insert_gfa1_link(tsv, &mut stmt_insert_link, &segments_by_name),
            _ => invalid_gfa!("GFA record type {} not yet supported", tsv[0]),
        }
    };

    // iterate tsv records
    util::iter_tsv_no_comments(dispatch, filename, Some('#' as u8))
}

fn insert_gfa1_segment(
    tsv: &Vec<&str>,
    txn: &Transaction,
    stmt: &mut Statement,
    maybe_stmt_rgfa: Option<&mut Statement>,
    segments_by_name: &mut HashMap<String, i64>,
) -> Result<()> {
    assert_eq!(tsv[0], "S");
    if tsv.len() < 2 {
        invalid_gfa!("malformed S line");
    }

    let rowid = name_to_rowid(tsv[1]);
    let name = if rowid.is_some() { None } else { Some(tsv[1]) };
    let sequence = if tsv.len() > 2 { Some(tsv[2]) } else { None };
    let mut tags_json = prepare_tags_json(tsv, 3)?;

    let rgfa_tags = if maybe_stmt_rgfa.is_some() {
        let sn = tags_json
            .remove("SN:Z")
            .map(|j| String::from(j.as_str().unwrap()));
        let so = tags_json.remove("SO:i").map(|j| j.as_i64().unwrap());
        let sr = tags_json.remove("SR:i").map(|j| j.as_i64().unwrap());
        if sn.is_none() || so.is_none() || sr.is_none() {
            invalid_gfa!(
                "segment missing required rGFA tags (SN:Z SO:i SR:i): {}",
                tsv[1]
            )
        }
        Some((sn.unwrap(), so.unwrap(), sr.unwrap()))
    } else {
        None
    };

    let tags_json_text = tags_json.dump();
    stmt.execute(params![rowid, name, tags_json_text, sequence])?;

    if let Some(nm) = name {
        segments_by_name.insert(String::from(nm), txn.last_insert_rowid());
    };

    if let Some(stmt_rgfa) = maybe_stmt_rgfa {
        let (sn, so, sr) = rgfa_tags.unwrap();
        stmt_rgfa.execute(params!(
            txn.last_insert_rowid(),
            sn,
            so,
            sequence.unwrap_or_default().len() as i64,
            sr
        ))?;
    };

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
        tags_json_text
    ])?;
    Ok(())
}

fn segment_and_orientation(
    segment: &str,
    orientation: &str,
    segments_by_name: &HashMap<String, i64>,
) -> Result<(i64, i64)> {
    let segment_id = if let Some(id) = name_to_rowid(segment) {
        id
    } else if let Some(idr) = segments_by_name.get(segment) {
        *idr
    } else {
        invalid_gfa!("unknown segment in link: {}", segment)
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

fn name_to_rowid(name: &str) -> Option<i64> {
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
                    let fv: f32 = fields[2].parse().or_else(|_| {
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
