use clap::Clap;
use genomicsqlite::ConnectionMethods;
use json::object;
use log::{debug, info, log_enabled, warn};
use num_format::{Locale, ToFormattedString};
use rusqlite::{params, OpenFlags, OptionalExtension, Statement, Transaction};
use std::cmp;
use std::collections::{HashMap, HashSet};

use crate::bad_command;
use crate::connectivity;
use crate::invalid_gfa;
use crate::util;
use crate::util::Result;

#[derive(Clap)]
pub struct Opts {
    /// Destination gfab filename
    #[clap(short)]
    pub output_gfab: String,

    /// Uncompressed GFA file/pipe [omit or - for standard input]
    #[clap(default_value = "-")]
    pub input_gfa: String,

    /// Always store segment/path name text (don't attempt to parse integer ID)
    #[clap(long)]
    pub always_names: bool,

    /// Omit index of graph connectivity (saves loading time & memory / disables certain queries)
    #[clap(long)]
    pub no_connectivity: bool,

    /// Omit segment sequences
    #[clap(long)]
    pub no_sequences: bool,

    /// Disable two-bit encoding for segment sequences (preserves lowercase nucleotides and U's / less efficient)
    #[clap(long)]
    pub no_twobit: bool,

    /// Memory budget (GiB)
    #[clap(long, default_value = "4")]
    pub memory_gbytes: u32,

    /// Compression level (-5 to 22)
    #[clap(long, default_value = "6")]
    pub compress: i8,

    /// log extra progress reports
    #[clap(short, long)]
    pub verbose: bool,

    /// log errors only
    #[clap(short, long)]
    pub quiet: bool,
}

pub fn main(opts: &Opts) -> Result<()> {
    if !opts.output_gfab.ends_with(".gfab") {
        warn!("output filename should end in .gfab")
    }
    if opts.input_gfa == "-" && atty::is(atty::Stream::Stdin) {
        bad_command!("pipe in .gfa data or supply input filename")
    }

    // formulate GenomicSQLite configuration JSON
    let mut db = new_db(
        &opts.output_gfab,
        opts.compress,
        std::cmp::max(1024, opts.memory_gbytes * 400),
    )?;

    let records_processed;
    {
        // open transaction & apply schema
        let txn = db.transaction()?;
        create_tables(&txn)?;

        // add temp tables for metadata, which we'll copy into the main db file after writing all
        // the segment sequences; this ensures the metadata is stored ~contiguously instead of
        // interspersed among the (typically much larger) sequence data.
        txn.execute_batch(
            "CREATE TABLE temp.segment_meta_hold(
                segment_id INTEGER PRIMARY KEY, name TEXT,
                sequence_length INTEGER, tags_json TEXT
            );
            CREATE TABLE temp.segment_mapping_hold(
                segment_id INTEGER NOT NULL,
                refseq_name TEXT NOT NULL,
                refseq_begin INTEGER NOT NULL,
                refseq_end INTEGER NOT NULL
            );
            CREATE TABLE temp.path_hold(
                path_id INTEGER PRIMARY KEY,
                name TEXT COLLATE UINT,
                tags_json TEXT
            );
            CREATE TABLE temp.walk_hold(
                walk_id INTEGER PRIMARY KEY,
                sample TEXT COLLATE UINT,
                hap_idx INTEGER NOT NULL,
                refseq_name TEXT NOT NULL COLLATE UINT,
                refseq_begin INTEGER NOT NULL,
                refseq_end INTEGER NOT NULL,
                min_segment_id INTEGER NOT NULL,
                max_segment_id INTEGER NOT NULL,
                tags_json TEXT
            )",
        )?;

        // intake GFA records
        debug!("processing GFA1 records...");
        records_processed = insert_gfa1(&opts.input_gfa, &txn, &opts)?;
        if records_processed == 0 {
            warn!("no input records processed")
        } else {
            info!("processed {} GFA1 record(s)", records_processed);
            debug!("writing metadata tables for Segments, Paths, and Walks...");
            // copy metadata as planned
            txn.execute_batch(
                "INSERT INTO gfa1_segment_meta(segment_id, name, sequence_length, tags_json)
                    SELECT segment_id, name, sequence_length, tags_json
                    FROM temp.segment_meta_hold;
                INSERT INTO gfa1_segment_mapping(segment_id, refseq_name, refseq_begin, refseq_end)
                    SELECT segment_id, refseq_name, refseq_begin, refseq_end
                    FROM temp.segment_mapping_hold ORDER BY segment_id;
                INSERT INTO gfa1_path(path_id, name, tags_json)
                    SELECT path_id, name, tags_json
                    FROM temp.path_hold;
                INSERT INTO gfa1_walk(walk_id, sample, hap_idx, refseq_name, refseq_begin, refseq_end,
                                      min_segment_id, max_segment_id, tags_json)
                    SELECT
                        walk_id, sample, hap_idx, refseq_name, refseq_begin, refseq_end,
                        min_segment_id, max_segment_id, tags_json
                    FROM temp.walk_hold",
            )?;
            debug!("insertions complete");
        }

        // indexing
        create_indexes(&txn, !opts.no_connectivity)?;

        // done
        debug!("flushing {} ...", &opts.output_gfab);
        txn.commit()?;
    }

    if log_enabled!(log::Level::Debug) {
        summary(&db)?;
    }
    db.close().map_err(|(_, e)| e)?;
    if records_processed > 0 {
        info!("🗹 done");
        Ok(())
    } else {
        Err(util::Error::EmptyGfab)
    }
}

pub fn new_db(
    filename: &str,
    compress: i8,
    page_cache_mebibytes: u32,
) -> Result<rusqlite::Connection> {
    // formulate GenomicSQLite configuration JSON
    let dbopts = match object! {
        unsafe_load: true,
        inner_page_KiB: 64,
        outer_page_KiB: 2,
        zstd_level: compress,
        page_cache_MiB: page_cache_mebibytes
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
    // "gfab" = 0x67666162
    db.execute_batch("PRAGMA foreign_keys = OFF; PRAGMA application_id = 0x67666162")?;
    Ok(db)
}

pub fn create_tables(db: &rusqlite::Connection) -> Result<()> {
    db.execute_batch(include_str!("schema/GFA1.sql"))?;
    debug!("created GFA1 tables");
    Ok(())
}

pub fn create_indexes(db: &rusqlite::Connection, connectivity: bool) -> Result<()> {
    info!("indexing...");

    let ddl = include_str!("schema/GFA1.index.sql");
    for index_spec in ddl.split(";") {
        let index_sql = index_spec.trim();
        if !index_sql.is_empty() {
            debug!("\t{} ...", index_sql.splitn(2, " ON").next().unwrap());
            db.execute_batch(index_sql)?;
        }
    }

    // add GRIs
    debug!("\tindexing segment mappings & walks by genomic range ...");
    for table in vec!["gfa1_segment_mapping", "gfa1_walk"] {
        let gri_sql =
            db.create_genomic_range_index_sql(table, "refseq_name", "refseq_begin", "refseq_end")?;
        db.execute_batch(&gri_sql)?;
    }

    if connectivity {
        debug!("\tindexing graph connectivity ...");
        connectivity::index(db)?;
    }

    debug!("\tANALYZE ...");
    db.execute_batch("PRAGMA analysis_limit = 1000; ANALYZE main")?;

    Ok(())
}

fn insert_gfa1(filename: &str, txn: &Transaction, opts: &Opts) -> Result<usize> {
    // prepared statements
    let mut stmt_insert_segment_meta =
        txn.prepare("INSERT INTO temp.segment_meta_hold(segment_id,name,sequence_length,tags_json) VALUES(?,?,?,?)")?;
    let mut stmt_insert_segment_sequence = txn.prepare(&format!(
        "INSERT INTO gfa1_segment_sequence(segment_id,sequence_twobit) VALUES(?,{})",
        if !opts.no_twobit {
            "nucleotides_twobit(?)"
        } else {
            "?"
        }
    ))?;
    let mut stmt_insert_link = txn.prepare(
        "INSERT INTO gfa1_link(from_segment,from_reverse,to_segment,to_reverse,cigar,tags_json) VALUES(?,?,?,?,?,?)"
    )?;
    let mut stmt_insert_segment_mapping = txn.prepare(
        "INSERT INTO temp.segment_mapping_hold(segment_id,refseq_name,refseq_begin,refseq_end) VALUES(?,?,?,?)"
    )?;
    let mut stmt_parse_rr = txn.prepare(
        "SELECT
            parse_genomic_range_sequence(?1),
            parse_genomic_range_begin(?1),
            parse_genomic_range_end(?1)",
    )?;
    let mut stmt_insert_path =
        txn.prepare("INSERT INTO temp.path_hold(path_id,name,tags_json) VALUES(?,?,?)")?;
    let mut stmt_insert_path_element = txn.prepare(
        "INSERT INTO gfa1_path_element(path_id,ordinal,segment_id,reverse,cigar_vs_previous) VALUES(?,?,?,?,?)"
    )?;
    let mut stmt_insert_walk =
        txn.prepare("INSERT INTO temp.walk_hold(sample,hap_idx,refseq_name,refseq_begin,refseq_end,min_segment_id,max_segment_id,tags_json) VALUES(?,?,?,?,?,?,?,?)")?;
    let mut stmt_insert_walk_steps =
        txn.prepare("INSERT INTO gfa1_walk_steps(walk_id,steps_jsarray) VALUES(?,?)")?;

    let mut segments_by_name = HashMap::new();
    let mut records: usize = 0;
    let mut header_records: usize = 0;
    let mut maybe_header = None;
    let mut sequence_char_warning = opts.no_twobit;

    // closure to process one record
    let mut other_record_types = HashSet::new();
    let dispatch = |line_num: usize, tsv: &Vec<&str>| -> Result<()> {
        match tsv[0] {
            "S" => {
                records += 1;
                insert_gfa1_segment(
                    line_num,
                    tsv,
                    txn,
                    !opts.no_sequences,
                    opts.always_names,
                    &mut stmt_insert_segment_meta,
                    &mut stmt_insert_segment_sequence,
                    &mut stmt_insert_segment_mapping,
                    &mut stmt_parse_rr,
                    &mut segments_by_name,
                    &mut sequence_char_warning,
                )
            }
            "L" => {
                records += 1;
                insert_gfa1_link(line_num, tsv, &mut stmt_insert_link, &segments_by_name)
            }
            "P" => {
                records += 1;
                insert_gfa1_path(
                    line_num,
                    tsv,
                    txn,
                    opts.always_names,
                    &mut stmt_insert_path,
                    &mut stmt_insert_path_element,
                    &segments_by_name,
                )
            }
            "C" => {
                panic!("gfabase hasn't yet implemented GFA Containment records; please bug the maintainers");
            }
            "H" => {
                records += 1;
                header_records += 1;
                if maybe_header.is_none() {
                    maybe_header = Some(prepare_tags_json(line_num, tsv, 1)?);
                } else if header_records == 2 {
                    warn!("ignored additional header (H) record(s) after the first");
                }
                Ok(())
            }
            "W" => {
                records += 1;
                insert_gfa1_walk(
                    line_num,
                    tsv,
                    txn,
                    opts.always_names,
                    &mut stmt_insert_walk,
                    &mut stmt_insert_walk_steps,
                    &segments_by_name,
                )
            }
            other => {
                if !other_record_types.contains(other) {
                    warn!("ignored record(s) with RecordType = {}", other);
                    other_record_types.insert(String::from(other));
                }
                Ok(())
            }
        }
    };

    // iterate tsv records
    util::iter_tsv_no_comments(dispatch, filename, Some('#' as u8))?;

    let mut header = maybe_header.unwrap_or(object::Object::new());
    header.insert(
        "PG:Z",
        json::JsonValue::from(format!("gfabase-v{}", env!("CARGO_PKG_VERSION"))),
    );
    txn.execute(
        "INSERT INTO gfa1_header(_rowid_, tags_json) VALUES(1, ?)",
        params![header.dump()],
    )?;

    Ok(records)
}

fn insert_gfa1_segment(
    line_num: usize,
    tsv: &Vec<&str>,
    txn: &Transaction,
    sequences: bool,
    always_names: bool,
    stmt_meta: &mut Statement,
    stmt_sequence: &mut Statement,
    stmt_mapping: &mut Statement,
    stmt_parse_rr: &mut Statement,
    segments_by_name: &mut HashMap<String, i64>,
    sequence_char_warning: &mut bool,
) -> Result<()> {
    assert_eq!(tsv[0], "S");
    if tsv.len() < 2 {
        invalid_gfa!("malformed S line");
    }

    let rowid = if !always_names {
        name_to_id(tsv[1])
    } else {
        None
    };
    let name = if rowid.is_some() { None } else { Some(tsv[1]) };
    let maybe_sequence = if tsv.len() > 2 && tsv[2] != "*" {
        Some(tsv[2])
    } else {
        None
    };
    let mut tags_json = prepare_tags_json(line_num, tsv, 3)?;

    // remove tag LN:i if present because we'll keep a dedicated column for this info (make sure
    // it's consistent)
    let ln_tag = tags_json.remove("LN:i").map(|j| j.as_i64()).flatten();
    let maybe_sequence_len = match (maybe_sequence, ln_tag) {
        (Some(seq), Some(lni)) if lni != (seq.len() as i64) => {
            invalid_gfa!(
                "(Ln {}) segment with inconsistent sequence length and LN tag: {}",
                line_num,
                tsv[1]
            )
        }
        (Some(seq), _) => Some(seq.len() as i64),
        (None, Some(lni)) => Some(lni),
        (None, None) => None,
    };

    let tags_json_text = tags_json.dump();
    stmt_meta.execute(params![rowid, name, maybe_sequence_len, tags_json_text])?;
    let rowid_actual = txn.last_insert_rowid();

    if let Some(nm) = name {
        segments_by_name.insert(String::from(nm), rowid_actual);
    }
    if sequences {
        if let Some(seq) = maybe_sequence {
            if !*sequence_char_warning {
                for ch in seq.chars() {
                    match ch {
                        'a' | 'c' | 'g' | 't' | 'u' | 'U' => {
                            warn!("segment sequences contain 'U' and/or lowercase nucleotides, which may not be preserved in the .gfab encoding (example segment_id = {})", rowid_actual);
                            *sequence_char_warning = true;
                            break;
                        }
                        _ => (),
                    }
                }
            }
            stmt_sequence.execute(params![rowid_actual, seq])?;
        }
    }

    // add a mapping from rGFA tags, if present
    let sn = tags_json
        .get("SN:Z")
        .map(|j| j.as_str().map(|s| String::from(s)))
        .flatten();
    let so = tags_json.get("SO:i").map(|j| j.as_i64()).flatten();
    match (sn, so, maybe_sequence_len) {
        (Some(refseq_name), Some(refseq_begin), Some(sequence_len)) => {
            stmt_mapping.execute(params!(
                rowid_actual,
                refseq_name,
                refseq_begin,
                refseq_begin + sequence_len,
            ))?;
        }
        _ => (),
    }

    // add a mapping from rr:Z tag (if present)
    if let Some(rr) = tags_json
        .get("rr:Z")
        .map(|j| j.as_str().map(|s| String::from(s)))
        .flatten()
    {
        stmt_parse_rr
            .query_row(params![rr], |row| {
                let refseq_name: String = row.get(0)?;
                let refseq_begin: i64 = row.get(1)?;
                let refseq_end: i64 = row.get(2)?;
                stmt_mapping.execute(params!(rowid_actual, refseq_name, refseq_begin, refseq_end,))
            })
            .map_err(|_| {
                util::Error::InvalidGfa(format!(
                    "(Ln {}) unable to parse rr:Z as genomic range (e.g. chr1:2,345-6,789): {}",
                    line_num, rr
                ))
            })?;
    }

    Ok(())
}

fn insert_gfa1_link(
    line_num: usize,
    tsv: &Vec<&str>,
    stmt: &mut Statement,
    segments_by_name: &HashMap<String, i64>,
) -> Result<()> {
    assert_eq!(tsv[0], "L");
    if tsv.len() < 5 {
        invalid_gfa!("malformed L line: {}", tsv.join("\t"));
    }

    let (from_segment, from_reverse) =
        segment_and_orientation(line_num, tsv[1], tsv[2], segments_by_name)?;
    let (to_segment, to_reverse) =
        segment_and_orientation(line_num, tsv[3], tsv[4], segments_by_name)?;
    let cigar = if tsv.len() > 5 && tsv[5] != "*" {
        Some(tsv[5])
    } else {
        None
    };
    let tags_json = prepare_tags_json(line_num, tsv, 6)?;
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
    line_num: usize,
    tsv: &Vec<&str>,
    txn: &Transaction,
    always_names: bool,
    stmt_path: &mut Statement,
    stmt_ele: &mut Statement,
    segments_by_name: &HashMap<String, i64>,
) -> Result<()> {
    assert_eq!(tsv[0], "P");
    if tsv.len() < 3 {
        invalid_gfa!("(Ln {}) malformed P line: {}", line_num, tsv.join("\t"));
    }

    let rowid = if !always_names {
        name_to_id(tsv[1])
    } else {
        None
    };
    let name = if rowid.is_some() { None } else { Some(tsv[1]) };

    let maybe_cigars: Option<Vec<&str>> = if tsv.len() > 3 && tsv[3] != "*" {
        Some(tsv[3].split(',').collect())
    } else {
        None
    };
    let tags_json = prepare_tags_json(line_num, tsv, 4)?;
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
            invalid_gfa!(
                "(Ln {}) incorrect overlaps: {} CIGARs != {} segments - 1",
                line_num,
                cigars.len(),
                segs.len()
            );
        }
    }
    for ord in 0..segs.len() {
        let ele = segs[ord];
        if ele.len() < 2 {
            invalid_gfa!("(Ln {}) malformed path: {}", line_num, tsv[2]);
        }
        let (segment_id, reverse) = segment_and_orientation(
            line_num,
            &ele[..(ele.len() - 1)],
            &ele[(ele.len() - 1)..ele.len()],
            segments_by_name,
        )?;
        let cigar = match (&maybe_cigars, ord) {
            (Some(cigars), i) if i > 0 => Some(cigars[i - 1]),
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

fn insert_gfa1_walk(
    line_num: usize,
    tsv: &Vec<&str>,
    txn: &Transaction,
    always_names: bool,
    stmt_walk: &mut Statement,
    stmt_steps: &mut Statement,
    segments_by_name: &HashMap<String, i64>,
) -> Result<()> {
    assert_eq!(tsv[0], "W");
    if tsv.len() < 7 {
        invalid_gfa!("(Ln {}) malformed W line: {}", line_num, tsv.join("\t"));
    }

    let sample = tsv[1];
    let hap_idx: u64 = match tsv[2].parse() {
        Ok(i) => i,
        Err(_) => invalid_gfa!("(Ln {}) malformed HapIdx: {}", line_num, tsv[2]),
    };
    let refseq_name = tsv[3];
    let refseq_begin: u64 = match tsv[4].parse() {
        Ok(i) => i,
        Err(_) => invalid_gfa!("(Ln {}) malformed SeqStart: {}", line_num, tsv[4]),
    };
    let refseq_end: u64 = match tsv[5].parse() {
        Ok(i) => i,
        Err(_) => invalid_gfa!("(Ln {}) malformed SeqEnd: {}", line_num, tsv[5]),
    };
    let tags_json = prepare_tags_json(line_num, tsv, 7)?;
    let tags_json_text = tags_json.dump();

    // JSON-encode the walk steps
    let mut first_step = true;
    let mut prev_segment_id = -1;
    let mut prev_reverse = false;
    let mut min_segment_id = i64::MAX;
    let mut max_segment_id = i64::MIN;
    let mut steps_json_text = String::from("[");
    for pre_step in tsv[6].split('>') {
        if first_step && pre_step.is_empty() {
            continue;
        }
        let mut reverse = false;
        for segment_name in pre_step.split('<') {
            if first_step {
                if segment_name.is_empty() {
                    continue;
                }
                steps_json_text.push_str("{\"");
            } else {
                steps_json_text.push_str(",{\"");
            }

            let maybe_segment_id = if always_names {
                None
            } else {
                name_to_id(segment_name)
            };
            let segment_id: i64 = match maybe_segment_id {
                Some(id) => id,
                None => match segments_by_name.get(segment_name) {
                    Some(id) => *id,
                    None => {
                        invalid_gfa!("(Ln {}) unknown segment name: {}", line_num, segment_name)
                    }
                },
            };

            let segment_id_text = format!("s\":{}", segment_id);
            // delta-encode the segment id, if that's shorter than writing it out
            let maybe_delta_text = if first_step {
                None
            } else {
                let (sign, delta) = if segment_id >= prev_segment_id {
                    ("+", segment_id - prev_segment_id)
                } else {
                    ("-", prev_segment_id - segment_id)
                };
                let delta_text = format!("{}\":{}", sign, delta);
                if delta_text.len() < segment_id_text.len() {
                    Some(delta_text)
                } else {
                    None
                }
            };
            steps_json_text.push_str(&maybe_delta_text.unwrap_or(segment_id_text));

            // write orientation if it's flipping wrt previous
            if first_step || prev_reverse != reverse {
                steps_json_text.push_str(if reverse { ",\"r\":1}" } else { ",\"r\":0}" })
            } else {
                steps_json_text.push_str("}")
            }

            prev_segment_id = segment_id;
            min_segment_id = cmp::min(min_segment_id, segment_id);
            max_segment_id = cmp::max(max_segment_id, segment_id);
            prev_reverse = reverse;
            // after the nested splits on '>' then '<', the first element of pre_step is forward
            // and remaining elements are reverse
            reverse = true;
            first_step = false;
        }
    }
    steps_json_text.push_str("]");

    if first_step {
        invalid_gfa!("(Ln {}) empty walk", line_num)
    }

    stmt_walk.execute(params![
        sample,
        hap_idx as i64,
        refseq_name,
        refseq_begin as i64,
        refseq_end as i64,
        min_segment_id,
        max_segment_id,
        if tags_json_text.trim() != "{}" {
            Some(tags_json_text)
        } else {
            None
        }
    ])?;
    let walk_id = txn.last_insert_rowid();

    stmt_steps.execute(params![walk_id, steps_json_text])?;
    Ok(())
}

fn segment_and_orientation(
    line_num: usize,
    segment: &str,
    orientation: &str,
    segments_by_name: &HashMap<String, i64>,
) -> Result<(i64, i64)> {
    let segment_id = if let Some(id) = name_to_id(segment) {
        id
    } else if let Some(idr) = segments_by_name.get(segment) {
        *idr
    } else {
        invalid_gfa!(
            "(Ln {}) unknown segment in link/path: {}",
            line_num,
            segment
        )
    };
    let reverse = match orientation {
        "+" => 0,
        "-" => 1,
        _ => {
            invalid_gfa!(
                "(Ln {}) malformed segment orientation: {}",
                line_num,
                orientation
            )
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

fn prepare_tags_json(
    line_num: usize,
    tsv: &Vec<&str>,
    offset: usize,
) -> Result<json::object::Object> {
    let mut ans = json::object::Object::new();
    if tsv.len() > offset {
        for cursor in offset..tsv.len() {
            if tsv[cursor].is_empty() && cursor == tsv.len() - 1 {
                continue;
            }
            let fields: Vec<&str> = tsv[cursor].splitn(3, ':').collect();
            if fields.len() != 3 || fields[0].is_empty() || fields[1].is_empty() {
                invalid_gfa!("(Ln {}) malformed tag: {}", line_num, tsv[cursor]);
            };
            let v = match fields[1] {
                "A" | "Z" | "H" => json::JsonValue::from(fields[2]),
                "i" => {
                    let iv: i64 = fields[2].parse().or_else(|_| {
                        invalid_gfa!("(Ln {}) malformed tag integer: {}", line_num, tsv[cursor]);
                    })?;
                    json::JsonValue::from(iv)
                }
                "f" => {
                    let fv: f64 = fields[2].parse().or_else(|_| {
                        invalid_gfa!("(Ln {}) malformed tag float: {}", line_num, tsv[cursor]);
                    })?;
                    json::JsonValue::from(fv)
                }
                // TODO: B & J
                _ => {
                    invalid_gfa!(
                        "(Ln {}) tag type not yet supported: {}",
                        line_num,
                        tsv[cursor]
                    );
                }
            };
            ans.insert(&format!("{}:{}", fields[0], fields[1]), v)
        }
    }
    Ok(ans)
}

pub fn summary(db: &rusqlite::Connection) -> Result<()> {
    debug!("tables & row counts:");
    let mut stmt_tables = db.prepare("SELECT name FROM sqlite_master WHERE type='table'")?;
    let mut tables = stmt_tables.query([])?;
    while let Some(row) = tables.next()? {
        let table: String = row.get(0)?;
        let ct: i64 = db.query_row(&format!("SELECT count(1) FROM {}", table), [], |ctr| {
            ctr.get(0)
        })?;
        debug!("\t{}\t{}", table, ct.to_formatted_string(&Locale::en));
    }
    if let Some(e) = db
        .query_row("PRAGMA foreign_key_check", [], |row| {
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
    if connectivity::has_index(db, "")? {
        debug!("undirected graph connectivity:");
        db.query_row(
            "SELECT count(component_id), max(size), sum(size), sum(cuts), sum(bp), sum(cuts_bp) FROM
                (SELECT
                        component_id, count(segment_id) AS size,
                        sum(is_cutpoint) AS cuts,
                        sum(sequence_length) AS bp,
                        sum(is_cutpoint*sequence_length) AS cuts_bp
                 FROM gfa1_connectivity INNER JOIN gfa1_segment_meta USING(segment_id)
                 GROUP BY component_id)",
            [],
            |row| {
                let count: i64 = row.get(0)?;
                if count > 0 {
                    let maxsize: i64 = row.get(1)?;
                    let sumsize: i64 = row.get(2)?;
                    let cuts: i64 = row.get(3)?;
                    let bp: i64 = row.get(4)?;
                    let cuts_bp: i64 = row.get(5)?;
                    debug!("\t{} connected components (of at least 2 segments)", count.to_formatted_string(&Locale::en));
                    debug!(
                        "\t{} segments in these components, totaling {} bp",
                        sumsize.to_formatted_string(&Locale::en), bp.to_formatted_string(&Locale::en)
                    );
                    debug!(
                        "\t{} cutpoint segments, totaling {} bp",
                        cuts.to_formatted_string(&Locale::en), cuts_bp.to_formatted_string(&Locale::en)
                    );
                    debug!("\t{} segments in largest component", maxsize.to_formatted_string(&Locale::en))
                } else {
                    warn!("graph has no links")
                }
                Ok(())
            },
        )?;
    }
    Ok(())
}
