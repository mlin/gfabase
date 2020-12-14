use clap::Clap;
use genomicsqlite::ConnectionMethods;
use json::JsonValue;
use log::{debug, error, info, warn};
use rusqlite::{params, OpenFlags, Statement, NO_PARAMS};
use std::{fs, io};

use crate::util;
use crate::util::Result;

#[derive(Clap)]
pub struct Opts {
    /// gfab filename
    pub gfab: String,
    /// output GFA file [stdout]
    pub gfa: Option<String>,
}

pub fn main(opts: &Opts) -> Result<()> {
    // formulate GenomicSQLite configuration JSON
    let mut dbopts = json::object::Object::new();
    dbopts.insert("immutable", json::JsonValue::from(true));

    // open db
    let db = genomicsqlite::open(
        &opts.gfab,
        OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
        &dbopts,
    )?;

    // open output writer
    let mut writer_box: Box<dyn io::Write> = match opts.gfa.as_ref().map(String::as_str) {
        None | Some("-") | Some("") => Box::new(io::BufWriter::new(io::stdout())),
        Some(p) => Box::new(io::BufWriter::new(fs::File::create(p)?)),
    };
    let writer = &mut *writer_box;
    write_segments(&db, writer)?;
    write_links(&db, writer)?;
    write_paths(&db, writer)?;
    writer.flush()?;

    Ok(())
}

fn write_segments(db: &rusqlite::Connection, writer: &mut dyn io::Write) -> Result<()> {
    let segments_query_sql = "SELECT
            segment_id, coalesce(name, cast(segment_id AS TEXT)), sequence_length,
            coalesce(tags_json, '{}'), sequence
            FROM gfa1_segment";
    let mut segments_query = db.prepare(segments_query_sql)?;
    let mut segments_cursor = segments_query.query(NO_PARAMS)?;
    while let Some(segrow) = segments_cursor.next()? {
        let rowid: i64 = segrow.get(0)?;
        let name: String = segrow.get(1)?;
        let maybe_sequence_length: Option<i64> = segrow.get(2)?;
        let tags_json: String = segrow.get(3)?;
        let sequence: Option<String> = segrow.get(4)?;
        writer.write_fmt(format_args!(
            "S\t{}\t{}",
            name,
            sequence.as_ref().map(String::as_str).unwrap_or("*")
        ))?;
        if let Some(sequence_length) = maybe_sequence_length {
            writer.write_fmt(format_args!("\tLN:i:{}", sequence_length))?;
        }
        write_tags("gfa1_segments_meta", rowid, &tags_json, writer)?;
        writer.write(b"\n")?;
    }
    Ok(())
}

fn write_links(db: &rusqlite::Connection, writer: &mut dyn io::Write) -> Result<()> {
    let mut links_query = db.prepare(
        // this two-layer join resolves the two segment IDs to names (if any)
        "SELECT
            link_id, from_segment_name, from_reverse,
            coalesce(gfa1_segment_meta.name, cast(to_segment AS TEXT)) AS to_segment_name,
            to_reverse, cigar, link_tags_json
        FROM
            (SELECT
                gfa1_link._rowid_ AS link_id,
                coalesce(gfa1_segment_meta.name, cast(from_segment AS TEXT)) AS from_segment_name,
                from_reverse, to_segment, to_reverse, coalesce(cigar, '*') AS cigar,
                coalesce(gfa1_link.tags_json, '{}') AS link_tags_json
            FROM
                gfa1_link LEFT JOIN gfa1_segment_meta ON from_segment = segment_id
            ORDER BY from_segment, to_segment)
            LEFT JOIN gfa1_segment_meta ON to_segment = segment_id",
    )?;
    let mut links_cursor = links_query.query(NO_PARAMS)?;
    while let Some(linkrow) = links_cursor.next()? {
        let link_id: i64 = linkrow.get(0)?;
        let from_segment: String = linkrow.get(1)?;
        let from_reverse: i8 = linkrow.get(2)?;
        let to_segment: String = linkrow.get(3)?;
        let to_reverse: i8 = linkrow.get(4)?;
        let cigar: String = linkrow.get(5)?;
        let tags_json: String = linkrow.get(6)?;
        writer.write_fmt(format_args!(
            "L\t{}\t{}\t{}\t{}\t{}",
            from_segment,
            if from_reverse == 0 { '+' } else { '-' },
            to_segment,
            if to_reverse == 0 { '+' } else { '-' },
            cigar
        ))?;
        write_tags("gfa1_link", link_id, &tags_json, writer)?;
        writer.write(b"\n")?;
    }
    Ok(())
}

fn write_paths(db: &rusqlite::Connection, writer: &mut dyn io::Write) -> Result<()> {
    let mut paths_query = db.prepare(
        "SELECT path_id, coalesce(name, cast(path_id AS TEXT)), coalesce(tags_json, '{}')
         FROM gfa1_path ORDER BY path_id",
    )?;
    let mut elements_query = db.prepare(
        "SELECT
            coalesce(name, cast(segment_id AS TEXT)) AS segment_name, reverse, cigar_vs_previous
         FROM gfa1_path_element LEFT JOIN gfa1_segment_meta USING(segment_id)
         WHERE path_id=? ORDER BY path_id, ordinal",
    )?;
    let mut paths_cursor = paths_query.query(NO_PARAMS)?;
    while let Some(pathrow) = paths_cursor.next()? {
        let path_id: i64 = pathrow.get(0)?;
        let name: String = pathrow.get(1)?;
        let tags_json: String = pathrow.get(2)?;

        let mut elts_csv = Vec::new();
        let mut cigars_csv = Vec::new();
        let mut elts_cursor = elements_query.query(params![path_id])?;
        while let Some(eltrow) = elts_cursor.next()? {
            let segment_name: String = eltrow.get(0)?;
            let reverse: i64 = eltrow.get(1)?;
            let maybe_cigar: Option<String> = eltrow.get(2)?;
            elts_csv.push(segment_name + if reverse == 0 { "+" } else { "-" });
            if let Some(cigar) = maybe_cigar {
                cigars_csv.push(cigar);
            }
        }

        writer.write_fmt(format_args!(
            "P\t{}\t{}\t{}",
            &name,
            &elts_csv.join(","),
            if cigars_csv.len() > 0 {
                cigars_csv.join(",")
            } else {
                String::from("*")
            }
        ))?;
        write_tags("gfa1_path", path_id, &tags_json, writer)?;
        writer.write(b"\n")?;
    }
    Ok(())
}

fn write_tags(table: &str, rowid: i64, tags_json: &str, writer: &mut dyn io::Write) -> Result<()> {
    let invalid = || util::Error::InvalidGfab {
        message: String::from("invalid tags_json"),
        table: String::from(table),
        rowid: rowid,
    };
    let tags = json::parse(&tags_json).map_err(|_| invalid())?;
    for (k, v) in tags.entries() {
        let kfields: Vec<&str> = k.split(':').collect();
        if kfields.len() != 2 {
            return Err(invalid());
        }

        let vstr = match kfields[1] {
            "A" | "Z" | "H" => JsonValue::as_str(v).ok_or_else(invalid)?.to_string(),
            "i" => JsonValue::as_i64(v).ok_or_else(invalid)?.to_string(),
            "f" => JsonValue::as_f64(v).ok_or_else(invalid)?.to_string(),
            // TODO: B & J
            _ => return Err(invalid()),
        };
        writer.write_fmt(format_args!("\t{}:{}", k, vstr))?;
    }
    Ok(())
}
