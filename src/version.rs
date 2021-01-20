use genomicsqlite::ConnectionMethods;
use rusqlite::OpenFlags;

pub const GFAB_VERSION_REQ: &str = ">= 0.4.0-0";

#[allow(unused)]
mod buildinfo {
    include!(concat!(env!("OUT_DIR"), "/built.rs"));
}

pub fn main() -> crate::util::Result<()> {
    let timestamp = chrono::DateTime::parse_from_rfc2822(buildinfo::BUILT_TIME_UTC).unwrap();
    println!(
        "gfabase {} v{} {}",
        buildinfo::PROFILE,
        env!("CARGO_PKG_VERSION"),
        timestamp.to_rfc3339()
    );
    println!(
        "  compatible .gfab format versions {}",
        semver::VersionReq::parse(GFAB_VERSION_REQ)
            .unwrap()
            .to_string()
    );
    println!(
        "{} RUSTFLAGS='{}'",
        buildinfo::RUSTC_VERSION,
        option_env!("RUSTFLAGS").unwrap_or("")
    );

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
    println!("GenomicSQLite {}", db.genomicsqlite_version());

    println!("SQLite v{} with compile options:", rusqlite::version());
    let mut compile_options_stmt = db.prepare("PRAGMA compile_options")?;
    let mut compile_options_cursor = compile_options_stmt.query(rusqlite::NO_PARAMS)?;
    while let Some(optrow) = compile_options_cursor.next()? {
        let opt: String = optrow.get(0)?;
        println!("  {}", &opt);
    }

    Ok(())
}
