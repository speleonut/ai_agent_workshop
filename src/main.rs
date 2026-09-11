//! `mytools` — a small reimplementation of a subset of bedtools.
//!
//! Only `--version` is wired up so far; subcommands come later.

use std::env;
use std::process::ExitCode;

const USAGE: &str = "usage: mytools --version\n";

fn main() -> ExitCode {
    let args: Vec<String> = env::args().skip(1).collect();

    match args.first().map(String::as_str) {
        Some("--version") => {
            println!("mytools {}", env!("CARGO_PKG_VERSION"));
            ExitCode::SUCCESS
        }
        // No arguments, or anything we don't recognise yet: usage error.
        // Usage goes to stderr because stdout is data and gets piped.
        _ => {
            eprint!("{USAGE}");
            ExitCode::from(2)
        }
    }
}
