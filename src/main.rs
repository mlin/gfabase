use clap::Clap;
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
    Load(load::Opts),
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
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
    }
}
