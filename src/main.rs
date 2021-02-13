use anyhow::Result;
use clap::Clap;
extern crate fern;
extern crate log;
use fern::colors::{Color, ColoredLevelConfig};
use log::error;

mod add_mappings;
mod connectivity;
mod load;
mod sub;
mod util;
mod version;
mod view;

#[derive(Clap)]
#[clap(version = env!("CARGO_PKG_VERSION"))]
struct Opts {
    #[clap(subcommand)]
    subcmd: SubCommand,

    /// log extra progress reports
    #[clap(short, long)]
    pub verbose: bool,

    /// log errors only
    #[clap(short, long)]
    pub quiet: bool,
}

#[derive(Clap)]
enum SubCommand {
    /// display more-detailed version info
    Version,

    /// in.gfa => out.gfab
    Load(load::Opts),

    /// assembly.gfab += mappings.{paf,gaf}
    AddMappings(add_mappings::Opts),

    /// in.gfab => out.gfa
    View(view::Opts),

    /// in.gfab => subgraph.gfab
    Sub(sub::Opts),
}

fn main() -> Result<()> {
    let mut opts = Opts::parse();
    if match &opts.subcmd {
        SubCommand::Version => false,
        SubCommand::Load(subopts) => subopts.verbose,
        SubCommand::AddMappings(subopts) => subopts.verbose,
        SubCommand::View(subopts) => subopts.verbose,
        SubCommand::Sub(subopts) => subopts.verbose,
    } {
        opts.verbose = true;
    }
    if match &opts.subcmd {
        SubCommand::Version => false,
        SubCommand::Load(subopts) => subopts.quiet,
        SubCommand::AddMappings(subopts) => subopts.quiet,
        SubCommand::View(subopts) => subopts.quiet,
        SubCommand::Sub(subopts) => subopts.quiet,
    } {
        opts.quiet = true;
    }

    let t0 = std::time::Instant::now();
    let colors = ColoredLevelConfig::new()
        .info(Color::Blue)
        .warn(Color::Yellow)
        .error(Color::BrightRed);
    let _logger = fern::Dispatch::new()
        .format(move |out, message, record| {
            let dur = t0.elapsed();
            out.finish(format_args!(
                "[{}][{}.{}s] {}",
                colors.color(record.level()),
                dur.as_secs(),
                dur.subsec_millis() / 100,
                message
            ))
        })
        .level(if opts.verbose {
            log::LevelFilter::Debug
        } else if opts.quiet {
            log::LevelFilter::Error
        } else {
            log::LevelFilter::Info
        })
        .chain(std::io::stderr())
        .apply()?;

    let rslt = match &opts.subcmd {
        SubCommand::Version => version::main(),
        SubCommand::Load(subopts) => load::main(subopts),
        SubCommand::AddMappings(subopts) => add_mappings::main(subopts),
        SubCommand::View(subopts) => view::main(subopts),
        SubCommand::Sub(subopts) => sub::main(subopts),
    };

    if let Err(util::Error::EmptyGfab) = rslt {
        error!("empty .gfab; exiting with code 3");
        std::process::exit(3)
    }

    Ok(rslt?)
}
