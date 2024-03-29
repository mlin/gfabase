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

    #[error(transparent)]
    Utf8Error(#[from] std::str::Utf8Error),

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
    F: FnMut(usize, X, &Vec<&str>) -> Result<X>,
{
    // https://stackoverflow.com/a/49964042/13393076
    let reader: Box<dyn io::BufRead> = if filename.is_empty() || filename == "-" {
        Box::new(io::BufReader::new(io::stdin()))
    } else {
        Box::new(io::BufReader::new(fs::File::open(filename)?))
    };

    let mut x = x0;
    let mut line_num = 0;
    for readline in reader.lines() {
        let line = readline?;
        line_num += 1;
        if comment.map_or(true, |ch| (line.is_empty() || line.as_bytes()[0] != ch)) {
            x = f(line_num, x, &line.split('\t').collect())?
        }
    }
    Ok(x)
}

// Iterate `f` over tab-separated lines of the file
pub fn iter_tsv_no_comments<F>(mut f: F, filename: &str, comment: Option<u8>) -> Result<()>
where
    F: FnMut(usize, &Vec<&str>) -> Result<()>,
{
    fold_tsv_no_comments(|line_num, (), tsv| f(line_num, tsv), (), filename, comment)
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
        [],
        |row| row.get(0),
    );
    if let Ok(pg) = pg_result {
        if pg.starts_with("gfabase-v") {
            // for semver comparison, strip front & back matter from the full pg version string
            let parts: Vec<&str> = pg.split('-').collect();
            if let Ok(v) = semver::Version::parse(&parts[1][1..]) {
                return Ok(v);
            }
        }
    }
    Err(Error::NotGfab)
}

pub fn check_gfab_version(gfab_version: &semver::Version) -> Result<()> {
    let req = semver::VersionReq::parse(GFAB_VERSION_REQ).unwrap();
    if req.matches(gfab_version) {
        let my_version_str = env!("CARGO_PKG_VERSION").split('-').next().unwrap();
        let my_version = semver::Version::parse(my_version_str).unwrap();
        if *gfab_version > my_version {
            warn!(
                "input .gfab from a newer version of gfabase ({} > {})",
                gfab_version.to_string(),
                my_version.to_string()
            );
        }
        return Ok(());
    }
    Err(Error::IncompatibleGfab {
        version: gfab_version.to_string(),
        required: req.to_string(),
    })
}

pub fn url_or_extant_file(it: &str) -> Result<()> {
    // not "safe", but usually gives more-helpful error message:
    if !it.starts_with("http:") && !it.starts_with("https:") && !Path::new(it).is_file() {
        return Err(Error::IoError(io::Error::new(
            io::ErrorKind::NotFound,
            "File not found",
        )));
    }
    Ok(())
}

pub fn open_gfab(
    filename: &str,
    flags: rusqlite::OpenFlags,
    genomicsqlite_config: &json::object::Object,
) -> Result<(semver::Version, rusqlite::Connection)> {
    url_or_extant_file(filename)?;
    match genomicsqlite::open(filename, flags, genomicsqlite_config) {
        Ok(db) => {
            let gfab_version = check_gfab_schema(&db, "")?;
            debug!("gfabase v{} created {}", gfab_version, filename);
            check_gfab_version(&gfab_version)?;
            Ok((gfab_version, db))
        }
        Err(err) => {
            debug!("check_gfab_filename_schema: {}", err);
            Err(Error::NotGfab)
        }
    }
}
