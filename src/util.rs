use crate::version::GFAB_VERSION_REQ;
use io::BufRead;
use log::{debug, warn};
use std::path::Path;
use std::{fs, io};
use thiserror::Error;

// helpful: https://nick.groenen.me/posts/rust-error-handling/
#[derive(Error, Debug)]
pub enum Error {
    #[error(transparent)]
    IoError(#[from] std::io::Error),

    #[error(transparent)]
    DbError(#[from] rusqlite::Error),

    #[error("[bad command] {0}")]
    BadCommand(String),

    #[error("[invalid GFA input] {0}")]
    InvalidGfa(String),

    #[error("[invalid PAF input] {0}")]
    InvalidPaf(String),

    #[error("[invalid .gfab][table = {table:?}, rowid = {rowid:?}] {message:?}")]
    InvalidGfab {
        message: String,
        table: String,
        rowid: i64,
    },

    #[error("file isn't .gfab format (or corrupt)")]
    NotGfab,

    #[error("[.gfab version incompatible] version: {version:?}, required: {required:?}")]
    IncompatibleGfab { version: String, required: String },

    #[error("empty .gfab")]
    EmptyGfab,
}

pub type Result<T> = std::result::Result<T, Error>;

#[macro_export]
macro_rules! invalid_gfa {
    ($($arg:tt)*) => ({
        return Err(util::Error::InvalidGfa(format!($($arg)*)));
    })
}

#[macro_export]
macro_rules! bad_command {
    ($($arg:tt)*) => ({
        return Err(util::Error::BadCommand(format!($($arg)*)));
    })
}

/// Fold over tab-separated lines of the file, excluding lines starting with specified comment
/// character, if any, e.g. `Some('#' as u8)`. Set `filename` empty to read standard input.
pub fn fold_tsv_no_comments<F, X>(mut f: F, x0: X, filename: &str, comment: Option<u8>) -> Result<X>
where
    F: FnMut(X, &Vec<&str>) -> Result<X>,
{
    // https://stackoverflow.com/a/49964042/13393076
    let reader: Box<dyn io::BufRead> = if filename.is_empty() || filename == "-" {
        Box::new(io::BufReader::new(io::stdin()))
    } else {
        Box::new(io::BufReader::new(fs::File::open(filename)?))
    };

    let mut x = x0;
    for readline in reader.lines() {
        let line = readline?;
        if comment.map_or(true, |ch| (line.is_empty() || line.as_bytes()[0] != ch)) {
            x = f(x, &line.split('\t').collect())?
        }
    }
    Ok(x)
}

// Iterate `f` over tab-separated lines of the file
pub fn iter_tsv_no_comments<F>(mut f: F, filename: &str, comment: Option<u8>) -> Result<()>
where
    F: FnMut(&Vec<&str>) -> Result<()>,
{
    fold_tsv_no_comments(|(), tsv| f(tsv), (), filename, comment)
}

pub fn delete_existing_file(filename: &str) -> Result<()> {
    debug!("deleting if present: {}", filename);
    let p = Path::new(filename);
    match std::fs::remove_file(p) {
        Ok(()) => {
            warn!("deleted existing file {}", filename);
            Ok(())
        }
        Err(ioerr) if ioerr.kind() == io::ErrorKind::NotFound => Ok(()),
        Err(ioerr) => Err(Error::IoError(ioerr)),
    }
}

pub fn check_gfab_schema(db: &rusqlite::Connection, schema: &str) -> Result<semver::Version> {
    let pg_result: rusqlite::Result<String> = db.query_row(
        &format!(
            "SELECT json_extract(tags_json, '$.PG:Z') FROM {}gfa1_header WHERE _rowid_ = 1",
            schema
        ),
        rusqlite::NO_PARAMS,
        |row| row.get(0),
    );
    if let Ok(pg) = pg_result {
        if pg.starts_with("gfabase-v") {
            if let Ok(v) = semver::Version::parse(&pg[9..]) {
                return Ok(v);
            }
        }
    }
    Err(Error::NotGfab)
}

pub fn check_gfab_version(gfab_version: &semver::Version) -> Result<()> {
    let req = semver::VersionReq::parse(GFAB_VERSION_REQ).unwrap();
    if req.matches(gfab_version) {
        return Ok(());
    }
    Err(Error::IncompatibleGfab {
        version: gfab_version.to_string(),
        required: req.to_string(),
    })
}

pub fn check_gfab_filename_schema(filename: &str) -> Result<semver::Version> {
    // not "safe", but usually gives more-helpful error message:
    if !Path::new(filename).is_file() {
        return Err(Error::IoError(io::Error::new(
            io::ErrorKind::NotFound,
            "File not found",
        )));
    }
    match genomicsqlite::open(
        filename,
        rusqlite::OpenFlags::SQLITE_OPEN_READ_ONLY | rusqlite::OpenFlags::SQLITE_OPEN_NO_MUTEX,
        &json::object::Object::new(),
    ) {
        Ok(db) => check_gfab_schema(&db, ""),
        _ => Err(Error::NotGfab),
    }
}
