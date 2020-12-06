use genomicsqlite::ConnectionMethods;
use rusqlite::OpenFlags;

pub fn main() -> crate::util::Result<()> {
    let tmp = tempfile::tempdir()?;
    let tmpdb = match tmp.path().join("version.genomicsqlite").to_str() {
        Some(p) => p.to_string(),
        None => panic!("invalid temp directory"),
    };
    let db = genomicsqlite::open(
        &tmpdb,
        OpenFlags::SQLITE_OPEN_CREATE
            | OpenFlags::SQLITE_OPEN_READ_WRITE
            | OpenFlags::SQLITE_OPEN_NO_MUTEX,
        &json::object::Object::new(),
    )?;
    println!("v{}\tSQLite", rusqlite::version());
    println!("{}\tGenomicSQLite", db.genomicsqlite_version());
    println!("v{}\tgfabase", env!("CARGO_PKG_VERSION"));
    Ok(())
}
