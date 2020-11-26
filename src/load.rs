use rusqlite::{params, Statement, Transaction};
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

type Result<T> = std::result::Result<T, Box<dyn Error>>;

fn ioerr<E>(kind: std::io::ErrorKind, error: E) -> Box<dyn Error>
where
    E: Into<Box<dyn Error + Send + Sync>>,
{
    Box::new(std::io::Error::new(kind, error))
}

fn fold_tsv_no_comments<F, X>(
    mut f: F,
    x0: X,
    filename_or_stdin: Option<&str>,
    comment: Option<u8>,
) -> Result<X>
where
    F: FnMut(X, &Vec<&str>) -> Result<X>,
{
    let reader: Box<dyn BufRead> = match filename_or_stdin {
        // https://stackoverflow.com/a/49964042/13393076
        None | Some("-") => Box::new(BufReader::new(io::stdin())),
        Some(filename) => Box::new(BufReader::new(File::open(filename)?)),
    };

    let mut x = x0;
    for readline in reader.lines() {
        let line = readline?;
        if !line.is_empty() && comment.map_or(true, |ch| line.as_bytes()[0] != ch) {
            x = f(x, &line.split('\t').collect())?
        }
    }
    Ok(x)
}

fn iter_tsv_no_comments<F>(
    mut f: F,
    filename_or_stdin: Option<&str>,
    comment: Option<u8>,
) -> Result<()>
where
    F: FnMut(&Vec<&str>) -> Result<()>,
{
    fold_tsv_no_comments(|(), tsv| f(tsv), (), filename_or_stdin, comment)
}

pub fn insert_gfa1(filename_or_stdin: Option<&str>, txn: &Transaction, prefix: &str) -> Result<()> {
    // prepared statements
    let mut stmt_insert_segment = txn.prepare(&format!(
        "INSERT INTO {}gfa1_segment(_rowid_,name,tags_json,sequence_twobit) VALUES(?,?,?,nucleotides_twobit(?))",
        prefix
    ))?;
    let mut stmt_insert_link = txn.prepare(&format!(
        "INSERT INTO {}gfa1_link(from_segment,from_reverse,to_segment,to_reverse,cigar,tags_json) VALUES(?,?,?,?,?,?)",
        prefix
    ))?;

    let mut segments_by_name = HashMap::new();

    // closure to process one record
    let dispatch = |tsv: &Vec<&str>| -> Result<()> {
        match tsv[0] {
            "S" => insert_gfa1_segment(tsv, txn, &mut stmt_insert_segment, &mut segments_by_name),
            "L" => insert_gfa1_link(tsv, &mut stmt_insert_link, &segments_by_name),
            _ => Ok(()),
        }
    };

    // iterate tsv records
    iter_tsv_no_comments(dispatch, filename_or_stdin, Some('#' as u8))
}

pub fn insert_gfa1_segment(
    tsv: &Vec<&str>,
    txn: &Transaction,
    stmt: &mut Statement,
    segments_by_name: &mut HashMap<String, i64>,
) -> Result<()> {
    assert_eq!(tsv[0], "S");
    if tsv.len() < 2 {
        return Err(ioerr(std::io::ErrorKind::InvalidInput, "malformed S line"));
    }

    let rowid = name_to_rowid(tsv[1]);
    let name = if rowid.is_some() { None } else { Some(tsv[1]) };
    let sequence = if tsv.len() > 2 { Some(tsv[2]) } else { None };
    let tags_json = prepare_tags_json(tsv, 3)?;
    let tags_json_text = tags_json.dump();
    stmt.execute(params![rowid, name, tags_json_text, sequence])?;

    if let Some(nm) = name {
        segments_by_name.insert(String::from(nm), txn.last_insert_rowid());
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
        return Err(ioerr(
            std::io::ErrorKind::InvalidInput,
            String::from("malformed L line: ") + &tsv.join("\t"),
        ));
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
    } else {
        *segments_by_name.get(segment).ok_or_else(|| {
            ioerr(
                std::io::ErrorKind::InvalidInput,
                String::from("unknown segment in link: ") + segment,
            )
        })?
    };
    let reverse = match orientation {
        "+" => 0,
        "-" => 1,
        _ => {
            return Err(ioerr(
                std::io::ErrorKind::InvalidInput,
                String::from("malformed segment orientation: ") + orientation,
            ));
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
                return Err(ioerr(
                    std::io::ErrorKind::InvalidInput,
                    String::from("malformed tag: ") + tsv[cursor],
                ));
            };
            let v = match fields[1] {
                "A" | "Z" | "H" => json::JsonValue::from(fields[2]),
                "i" => {
                    let iv: i64 = fields[2].parse().map_err(|_| {
                        ioerr(
                            std::io::ErrorKind::InvalidInput,
                            String::from("malformed tag integer: ") + tsv[cursor],
                        )
                    })?;
                    json::JsonValue::from(iv)
                }
                "f" => {
                    let fv: f32 = fields[2].parse().map_err(|_| {
                        ioerr(
                            std::io::ErrorKind::InvalidInput,
                            String::from("malformed tag float: ") + tsv[cursor],
                        )
                    })?;
                    json::JsonValue::from(fv)
                }
                // TODO: B & J
                _ => {
                    return Err(ioerr(
                        std::io::ErrorKind::Other,
                        String::from("tag type not yet supported: ") + tsv[cursor],
                    ))
                }
            };
            ans.insert(&format!("{}:{}", fields[0], fields[1]), v)
        }
    }
    Ok(ans)
}
