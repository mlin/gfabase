use clap::Clap;
use genomicsqlite::ConnectionMethods;
use rusqlite::{params, Connection, OpenFlags, NO_PARAMS};
use std::path::Path;
#[macro_use]
extern crate log;
extern crate fern;
use fern::colors::{Color, ColoredLevelConfig};

mod load;

#[derive(Clap)]
struct Opts {
    #[clap(subcommand)]
    subcmd: SubCommand,
}

#[derive(Clap)]
enum SubCommand {
    Load(LoadOpts),
}

#[derive(Clap)]
struct LoadOpts {
    /// destination gfab filename
    gfab: String,
    /// input GFA file (omit for standard input)
    gfa: Option<String>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let colors = ColoredLevelConfig::new()
        .debug(Color::Blue)
        .warn(Color::Yellow)
        .error(Color::BrightRed);
    let _logger = fern::Dispatch::new()
        .format(move |out, message, record| {
            out.finish(format_args!(
                "[{}] {}",
                colors.color(record.level()),
                message
            ))
        })
        .level(log::LevelFilter::Debug)
        .chain(std::io::stderr())
        .apply()?;
    let opts = Opts::parse();

    match opts.subcmd {
        SubCommand::Load(subopts) => Ok(main_load(&subopts)?),
    }
}

fn main_load(opts: &LoadOpts) -> Result<(), Box<dyn std::error::Error>> {
    // formulate GenomicSQLite configuration JSON
    let mut dbopts = json::object::Object::new();
    dbopts.insert("unsafe_load", json::JsonValue::from(true));
    dbopts.insert("inner_page_KiB", json::JsonValue::from(64));
    dbopts.insert("outer_page_KiB", json::JsonValue::from(2));

    // delete existing file, if any
    {
        let p = Path::new(&opts.gfab);
        if p.is_file() {
            warn!("delete existing destination file: {}", p.to_str().unwrap());
            std::fs::remove_file(p)?
        }
    }
    // create db
    let mut db = genomicsqlite::open(
        &opts.gfab,
        OpenFlags::SQLITE_OPEN_CREATE | OpenFlags::SQLITE_OPEN_READ_WRITE,
        &dbopts,
    )?;

    // open transaction & apply schema
    let txn = db.transaction()?;
    let mut gfa_sql = include_str!("schema/GFA1.sql").to_string();
    gfa_sql = simple_replace_all(&gfa_sql, "prefix", "");
    println!("{}", gfa_sql);
    txn.execute_batch(&*gfa_sql)?;

    // intake GFA records
    load::insert_gfa1(opts.gfa.as_deref(), &txn, "")?;

    // done
    txn.commit()?;
    Ok(())
}

fn simple_replace_all(template: &str, key: &str, val: &str) -> String {
    let pat = r"\{\{\s*".to_string() + key + r"\s*\}\}";
    regex::Regex::new(&pat)
        .unwrap()
        .replace_all(template, val)
        .to_string()
}
