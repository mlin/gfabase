use clap::Clap;
use genomicsqlite::ConnectionMethods;
use json::JsonValue;
use log::{debug, error, info, warn};
use rusqlite::{params, OpenFlags, Statement, NO_PARAMS};
use std::io::Write;
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

// segment selection criteria:
// - connected component of {segment(s), path}
//   - selectable outgoing/incoming/bidirectional
//   - max degree
// - rGFA reference range
//   - reference only
//   - connected component (out/in/bi)
//   - follow rank>0 arcs until rejoining rank==0 (to/from/by)
//
// idea: optionally include rank>0 segments in rGFA GRI, based on the widest linear range reachable
// by passing other rank>0 segments only. Then we could find alleles that "contain" a query segment
// even if they don't start or end within. Build recursively starting from nodes directly connected
// to reference segments. This would allow us to serve linear range queries without link-traversal
// loops.
// To seed coordinates, propagate (i) the begin_pos of each R0 segment along each of its outgoing
// links; (ii) the end_pos of each R0 segment along each of its incoming links. We can only know
// the coordinates for a R>0 segment once all of its edges have been thusly filled in.
// (speed heuristic -- "page in" the links table in slabs)

pub fn main(opts: &Opts) -> Result<()> {
    let prefix = "";

    // formulate GenomicSQLite configuration JSON
    let mut dbopts = json::object::Object::new();
    dbopts.insert("immutable", json::JsonValue::from(true));

    // open db
    let db = genomicsqlite::open(&opts.gfab, OpenFlags::SQLITE_OPEN_READ_ONLY, &dbopts)?;

    // open output writer
    let mut writer: Box<dyn io::Write> = match opts.gfa.as_ref().map(String::as_str) {
        None | Some("-") | Some("") => Box::new(io::BufWriter::new(io::stdout())),
        Some(p) => Box::new(io::BufWriter::new(fs::File::create(p)?)),
    };

    // output segments
    let segments_table = String::from(prefix) + "gfa1_segment";
    let mut segments_query = db.prepare(&format!(
        "SELECT _rowid_, name, tags_json, twobit_dna(sequence_twobit) FROM {} ORDER BY _rowid_",
        segments_table
    ))?;
    let mut segments_cursor = segments_query.query(NO_PARAMS)?;
    while let Some(segrow) = segments_cursor.next()? {
        let rowid: i64 = segrow.get(0)?;
        let name: Option<String> = segrow.get(1)?;
        let tags_json: String = segrow.get(2)?;
        let sequence: Option<String> = segrow.get(3)?;
        // TODO: add rGFA tags back in, if needed
        writer.write_fmt(format_args!(
            "S\t{}\t{}",
            name.unwrap_or(rowid.to_string()),
            sequence.as_ref().map(String::as_str).unwrap_or("*")
        ))?;
        write_tags(&segments_table, rowid, &tags_json, &mut writer)?;
        writer.write(b"\n")?;
    }

    writer.flush()?;
    Ok(())
}

fn write_tags(
    table: &str,
    rowid: i64,
    tags_json: &str,
    writer: &mut Box<dyn io::Write>,
) -> Result<()> {
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
            "f" => JsonValue::as_f32(v).ok_or_else(invalid)?.to_string(),
            // TODO: B & J
            _ => return Err(invalid()),
        };
        writer.write_fmt(format_args!("\t{}:{}", k, vstr))?;
    }
    Ok(())
}
