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

    #[error("[Bad command] {0}")]
    BadCommand(String),

    #[error("[Invalid GFA input] {0}")]
    InvalidGfa(String),

    #[error("[Invalid .gfab][table = {table:?}, rowid = {rowid:?}] {message:?}")]
    InvalidGfab {
        message: String,
        table: String,
        rowid: i64,
    },

    #[error("Empty .gfab")]
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
