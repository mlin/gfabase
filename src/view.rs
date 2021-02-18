use clap::Clap;
use io::Write;
use json::JsonValue;
use log::{info, warn};
use num_format::{Locale, ToFormattedString};
use rusqlite::{params, OpenFlags, OptionalExtension, NO_PARAMS};
use std::{env, fs, io, path, process};

use crate::bad_command;
use crate::util;
use crate::util::Result;

#[derive(Clap)]
pub struct Opts {
    /// gfab filename or http[s] URL
    pub input_gfab: String,
    /// output GFA filename [omit or - for standard output]
    #[clap(short, default_value = "-")]
    pub output_gfa: String,
    /// Omit segment sequences
    #[clap(long)]
    pub no_sequences: bool,
    /// Launch Bandage on output file (temporary file, if unspecified)
    #[clap(long)]
    pub bandage: bool,
    /// For each segment with reference mappings, set gr:Z tag with one guessed range summarizing the mappings
    #[clap(long)]
    pub guess_ranges: bool,

    /// log extra progress reports
    #[clap(short, long)]
    pub verbose: bool,

    /// log errors only
    #[clap(short, long)]
    pub quiet: bool,
}

pub fn main(opts: &Opts) -> Result<()> {
    // formulate GenomicSQLite configuration JSON
    let mut dbopts = json::object::Object::new();
    dbopts.insert("immutable", json::JsonValue::from(true));

    // open db
    let (_gfab_version, mut db) = util::open_gfab(
        &opts.input_gfab,
        OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
        &dbopts,
    )?;

    {
        let txn = db.transaction()?;
        let mut maybe_guesser = if opts.guess_ranges {
            Some(SegmentRangeGuesser::new(&txn, "")?)
        } else {
            None
        };
        let mut tag_editor = |segment_id: i64, tags: &mut json::JsonValue| -> Result<()> {
            if let Some(ref mut guesser) = maybe_guesser {
                if let Some(gr) = guesser.get(segment_id)? {
                    tags.insert("gr:Z", gr).unwrap()
                }
            }
            Ok(())
        };

        if opts.output_gfa == "-" && !opts.bandage && atty::is(atty::Stream::Stdout) {
            // interactive mode: pipe into less -S
            less(|less_in| {
                write_header(&txn, less_in)
                    .and_then(|_| {
                        write_segments(&txn, "", !opts.no_sequences, &mut tag_editor, less_in)
                    })
                    .and_then(|_| write_links(&txn, "", less_in))
                    .and_then(|_| write_paths(&txn, "", less_in))
            })?
        } else {
            let mut output_gfa = String::from(&opts.output_gfa);
            if opts.bandage && output_gfa == "-" {
                output_gfa = bandage_temp_filename()?
            }

            {
                let mut writer_box = writer(&output_gfa)?;
                let out = &mut *writer_box;
                write_header(&txn, out)?;
                write_segments(&txn, "", !opts.no_sequences, &mut tag_editor, out)?;
                write_links(&txn, "", out)?;
                write_paths(&txn, "", out)?
            }

            if opts.bandage {
                if let Some(ref mut guesser) = maybe_guesser {
                    guesser.write_bandage_csv(&output_gfa)?
                }
                bandage(&output_gfa)?
            }
        }
    }

    Ok(())
}

pub fn writer(gfa_filename: &str) -> Result<Box<dyn io::Write>> {
    if gfa_filename.is_empty() || gfa_filename == "-" {
        return Ok(Box::new(io::BufWriter::new(io::stdout())));
    }
    if !gfa_filename.ends_with(".gfa") {
        warn!("output filename should end with .gfa")
    }
    Ok(Box::new(io::BufWriter::new(fs::File::create(
        gfa_filename,
    )?)))
}

/// Start `less -S` and call `write` with its standard input pipe.
/// Tolerate BrokenPipe errors (user exited before viewing all data)
pub fn less<F>(write: F) -> Result<()>
where
    F: FnOnce(&mut dyn io::Write) -> Result<()>,
{
    if which::which("less").is_err() {
        return write(&mut io::stdout());
    }
    let mut child = process::Command::new("less")
        .arg("-S")
        .stdin(process::Stdio::piped())
        .spawn()?;
    {
        let mut less_in = child.stdin.take().unwrap();

        match write(&mut less_in) {
            Ok(()) => (),
            Err(util::Error::IoError(err)) if err.kind() == io::ErrorKind::BrokenPipe => (),
            Err(e) => return Err(e),
        }
    }
    child.wait()?;
    Ok(())
}

pub fn bandage(gfa: &str) -> Result<()> {
    info!("Bandage load {} --draw", gfa);
    if process::Command::new("Bandage")
        .arg("load")
        .arg(gfa)
        .arg("--draw")
        .spawn()
        .is_err()
    {
        bad_command!("failed to launch Bandage; make sure it's in PATH")
    }
    Ok(())
}

pub fn bandage_temp_filename() -> Result<String> {
    if !atty::is(atty::Stream::Stdout) {
        bad_command!("supply -o filename.gfa on which to launch Bandage")
    }
    Ok(String::from(
        path::Path::new(&env::var("TMPDIR").unwrap_or(String::from("/tmp")))
            .join(format!(
                "gfabase-bandage-{}.gfa",
                chrono::Local::now().to_rfc3339_opts(chrono::SecondsFormat::Micros, true)
            ))
            .to_str()
            .unwrap(),
    ))
}

pub fn write_header(db: &rusqlite::Connection, writer: &mut dyn io::Write) -> Result<()> {
    let tags_json: String = db.query_row(
        "SELECT tags_json FROM gfa1_header WHERE _rowid_ = 1",
        NO_PARAMS,
        |row| row.get(0),
    )?;
    writer.write(b"H")?;
    write_tags("gfa1_header", 1, &tags_json, writer)?;
    writer.write(b"\n")?;
    Ok(())
}

pub fn write_segments(
    db: &rusqlite::Connection,
    where_clause: &str,
    with_sequences: bool,
    mut tag_editor: impl FnMut(i64, &mut json::JsonValue) -> Result<()>,
    writer: &mut dyn io::Write,
) -> Result<()> {
    let segments_query_sql = String::from(if with_sequences {
        "SELECT
                segment_id, coalesce(name, cast(segment_id AS TEXT)), sequence_length,
                coalesce(tags_json, '{}'), sequence
             FROM gfa1_segment "
    } else {
        "SELECT
                segment_id, coalesce(name, cast(segment_id AS TEXT)),
                sequence_length, coalesce(tags_json, '{}')
             FROM gfa1_segment_meta "
    }) + where_clause;
    let mut segments_query = db.prepare(&segments_query_sql)?;
    let mut segments_cursor = segments_query.query(NO_PARAMS)?;
    while let Some(segrow) = segments_cursor.next()? {
        let rowid: i64 = segrow.get(0)?;
        let name: String = segrow.get(1)?;
        let maybe_sequence_length: Option<i64> = segrow.get(2)?;
        let tags_json: String = segrow.get(3)?;
        let sequence: Option<String> = if with_sequences { segrow.get(4)? } else { None };
        writer.write_fmt(format_args!(
            "S\t{}\t{}",
            name,
            sequence.as_ref().map(String::as_str).unwrap_or("*")
        ))?;
        if let Some(sequence_length) = maybe_sequence_length {
            writer.write_fmt(format_args!("\tLN:i:{}", sequence_length))?;
        }
        write_tags_with_editor(
            "gfa1_segments_meta",
            rowid,
            &tags_json,
            &mut tag_editor,
            writer,
        )?;
        writer.write(b"\n")?;
    }
    Ok(())
}

pub fn write_links(
    db: &rusqlite::Connection,
    where_clause: &str,
    writer: &mut dyn io::Write,
) -> Result<()> {
    let links_query_sql = format!(
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
                coalesce(gfa1_link.tags_json, '{{}}') AS link_tags_json
            FROM
                gfa1_link LEFT JOIN gfa1_segment_meta ON from_segment = segment_id
            {}
            ORDER BY from_segment, to_segment)
            LEFT JOIN gfa1_segment_meta ON to_segment = segment_id",
        where_clause
    );
    let mut links_query = db.prepare(&links_query_sql)?;
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

pub fn write_paths(
    db: &rusqlite::Connection,
    where_clause: &str,
    writer: &mut dyn io::Write,
) -> Result<()> {
    let paths_query_sql = format!(
        "SELECT path_id, coalesce(name, cast(path_id AS TEXT)), coalesce(tags_json, '{{}}')
         FROM gfa1_path {} ORDER BY path_id",
        where_clause
    );
    let mut paths_query = db.prepare(&paths_query_sql)?;
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

fn write_tags_with_editor(
    table: &str,
    rowid: i64,
    tags_json: &str,
    mut editor: impl FnMut(i64, &mut json::JsonValue) -> Result<()>,
    writer: &mut dyn io::Write,
) -> Result<()> {
    let invalid = || util::Error::InvalidGfab {
        message: String::from("invalid tags_json"),
        table: String::from(table),
        rowid: rowid,
    };
    let mut tags = json::parse(&tags_json).map_err(|_| invalid())?;
    editor(rowid, &mut tags)?;
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

fn write_tags(table: &str, rowid: i64, tags_json: &str, writer: &mut dyn io::Write) -> Result<()> {
    write_tags_with_editor(table, rowid, tags_json, |_, _| Ok(()), writer)
}

// Helpers roughly guessing a genomic range for a segment based on its PAF mappings. Selects the
// chromosome with the most coverage in the mappings, then the min and max mapped position on that
// chromosome.
pub struct SegmentRangeGuesser<'a> {
    getter: rusqlite::Statement<'a>,
    csv_query: rusqlite::Statement<'a>,
}

impl<'a> SegmentRangeGuesser<'_> {
    pub fn new(
        db: &'a rusqlite::Connection,
        where_clause: &str,
    ) -> Result<SegmentRangeGuesser<'a>> {
        // analyze mappings to generate temp.segment_range_guess
        db.execute(
            "CREATE TABLE temp.segment_range_guess(
                segment_id INTEGER PRIMARY KEY,
                refseq_name TEXT NOT NULL,
                refseq_begin INTEGER NOT NULL, refseq_end INTEGER NOT NULL)",
            NO_PARAMS,
        )?;
        let sql = format!(
            "WITH summary AS
                (SELECT
                    segment_id, refseq_name,
                    min(refseq_begin) AS min_begin, max(refseq_end) AS max_end,
                    max(refseq_end) - min(refseq_begin) AS coverage,
                    sum(refseq_end - refseq_begin) AS coverage2
                FROM gfa1_segment_mapping
                {}
                GROUP BY segment_id, refseq_name)
             INSERT INTO temp.segment_range_guess(segment_id, refseq_name, refseq_begin, refseq_end)
                SELECT segment_id, refseq_name, min_begin, max_end
                FROM
                    (SELECT
                        segment_id, refseq_name, min_begin, max_end,
                        row_number() OVER (PARTITION BY segment_id ORDER BY coverage DESC, coverage2 DESC)
                            AS coverage_rank
                    FROM summary)
                WHERE coverage_rank = 1",
            where_clause
        );
        let n = db.execute(&sql, NO_PARAMS)?;
        info!("guessed ranges for {} segments", n);
        // prepare queries on temp.segment_range_guess
        Ok(SegmentRangeGuesser {
            getter: db.prepare(
                "SELECT refseq_name, refseq_begin, refseq_end
                 FROM temp.segment_range_guess WHERE segment_id = ?",
            )?,
            csv_query: db.prepare(
                "SELECT
                    coalesce(name, cast(segment_id AS TEXT)),
                    refseq_name, refseq_begin, refseq_end
                 FROM temp.segment_range_guess LEFT JOIN gfa1_segment_meta USING(segment_id)",
            )?,
        })
    }

    pub fn get(&mut self, segment_id: i64) -> Result<Option<String>> {
        let maybe_row: Option<(String, i64, i64)> = self
            .getter
            .query_row(params![segment_id], |row| {
                Ok((row.get(0)?, row.get(1)?, row.get(2)?))
            })
            .optional()?;
        if let Some((refseq_name, refseq_begin, refseq_end)) = maybe_row {
            return Ok(Some(format!(
                "~{}:{}-{}",
                refseq_name,
                (refseq_begin + 1).to_formatted_string(&Locale::en),
                refseq_end.to_formatted_string(&Locale::en)
            )));
        }
        Ok(None)
    }

    pub fn write_bandage_csv(&mut self, gfa_filename: &str) -> Result<()> {
        // write a CSV file with the guessed ranges that Bandage can show as labels
        let csv_filename = String::from(gfa_filename.strip_suffix(".gfa").unwrap_or(gfa_filename))
            + ".guessed_ranges.csv";
        {
            let mut writer = io::BufWriter::new(fs::File::create(&csv_filename)?);
            writer.write_fmt(format_args!("Name,Guessed range\n"))?;
            let mut cursor = self.csv_query.query(NO_PARAMS)?;
            while let Some(row) = cursor.next()? {
                let name: String = row.get(0)?;
                let refseq_name: String = row.get(1)?;
                let refseq_begin: i64 = row.get(2)?;
                let refseq_end: i64 = row.get(3)?;
                writer.write_fmt(format_args!(
                    "\"{}\",\"~{}:{}-{}\"\n",
                    name,
                    refseq_name,
                    (refseq_begin + 1).to_formatted_string(&Locale::en),
                    refseq_end.to_formatted_string(&Locale::en)
                ))?;
            }
        }
        info!("wrote CSV with guessed segment ranges to {}", csv_filename);
        Ok(())
    }
}
