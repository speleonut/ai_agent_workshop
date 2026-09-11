//! `mytools` — a small reimplementation of a subset of bedtools.
//!
//! `sort`, `merge` and `intersect`. Real `bedtools` is the oracle: where our
//! output differs from it on the same input, we are wrong (`SPEC.md` §8).
//!
//! Layout: `cli` parses the command line, `bed` turns a source into validated
//! records, and one module per subcommand does the work. stdout is data;
//! everything diagnostic goes to stderr.

mod bed;
mod cli;
mod intersect;
mod merge;
mod sort;

use std::env;
use std::process::ExitCode;

/// Exit codes (`SPEC.md` §7): `0` success, `1` bad input data, `2` usage error.
fn main() -> ExitCode {
    let args: Vec<String> = env::args().skip(1).collect();

    let command = match cli::parse(&args) {
        Ok(c) => c,
        Err(usage) => {
            eprint!("{usage}");
            return ExitCode::from(2);
        }
    };

    let result = match command {
        cli::Command::Version => {
            println!("mytools {}", env!("CARGO_PKG_VERSION"));
            return ExitCode::SUCCESS;
        }
        cli::Command::Sort(opts) => sort::run(&opts),
        cli::Command::Merge(opts) => merge::run(&opts),
        cli::Command::Intersect(opts) => intersect::run(&opts),
    };

    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("{e}");
            ExitCode::from(1)
        }
    }
}
