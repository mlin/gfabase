use anyhow::Result;
use clap::Clap;
extern crate fern;
extern crate log;
use fern::colors::{Color, ColoredLevelConfig};

mod load;
mod util;
mod view;

#[derive(Clap)]
struct Opts {
    #[clap(subcommand)]
    subcmd: SubCommand,
}

#[derive(Clap)]
enum SubCommand {
    /// .gfa => .gfab
    Load(load::Opts),

    /// .gfab => .gfa
    View(view::Opts),
}

fn main() -> Result<()> {
    let t0 = std::time::Instant::now();
    let colors = ColoredLevelConfig::new()
        .debug(Color::Blue)
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
        .level(log::LevelFilter::Debug)
        .chain(std::io::stderr())
        .apply()?;
    let opts = Opts::parse();

    match opts.subcmd {
        SubCommand::Load(subopts) => Ok(load::main(&subopts)?),
        SubCommand::View(subopts) => Ok(view::main(&subopts)?),
    }
}
