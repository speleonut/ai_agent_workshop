//! Command line front end: dispatch, flags, usage errors, warnings, `--quiet`.
//!
//! `SPEC.md` §4, §7. Flag names and meanings match bedtools exactly. Anything
//! not listed here is a usage error — never silently ignored.

// TODO: drop this once sort, merge and intersect are all implemented.
#![allow(dead_code)]

use std::fmt;

pub const USAGE: &str = "\
usage: mytools <subcommand> [options]

subcommands:
  sort      -i <file> [-header]
  merge     -i <file> [-d N] [-s] [-S +|-] [-header]
  intersect -a <file> -b <file> [-u] [-v] [-wa] [-header]

global:
  --version   print version and exit
  --quiet     silence warnings

`-` as a filename means stdin.
";

/// A usage error: exit `2` (`SPEC.md` §7). This is a deliberate, documented
/// deviation from bedtools, which exits `1` for everything — see §8.1. Do not
/// "fix" it to match the oracle.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UsageError(pub Option<String>);

impl UsageError {
    pub fn new(msg: impl Into<String>) -> UsageError {
        UsageError(Some(msg.into()))
    }
    /// Bare `mytools`: usage alone, no message.
    pub fn bare() -> UsageError {
        UsageError(None)
    }
}

impl fmt::Display for UsageError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.0 {
            Some(m) => write!(f, "mytools: {m}\n{USAGE}"),
            None => write!(f, "{USAGE}"),
        }
    }
}

/// Write a warning to stderr. Warnings never change the exit code and never
/// touch stdout, so golden diffs cannot see them (`SPEC.md` §7).
pub fn warn(quiet: bool, msg: &str) {
    if !quiet {
        eprintln!("mytools: warning: {msg}");
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SortOpts {
    /// `-i`: a path, or `-` for stdin.
    pub input: String,
    pub header: bool,
    pub quiet: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MergeOpts {
    pub input: String,
    /// `-d`: maximum distance still merged. Negative requires that much overlap.
    pub d: i64,
    /// `-s`: merge only features sharing a strand.
    pub same_strand: bool,
    /// `-S +` / `-S -`: merge only features on the named strand.
    pub strand: Option<char>,
    pub header: bool,
    pub quiet: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IntersectOpts {
    pub a: String,
    pub b: String,
    /// `-u`: report each overlapping A record once, at original width.
    pub u: bool,
    /// `-v`: report A records with no overlap in B.
    pub v: bool,
    /// `-wa`: report the original A record instead of the clipped overlap.
    pub wa: bool,
    pub header: bool,
    pub quiet: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Command {
    Version,
    Sort(SortOpts),
    Merge(MergeOpts),
    Intersect(IntersectOpts),
}

/// Pull the value that follows a flag, or report the flag as needing one.
fn value(flag: &str, it: &mut std::slice::Iter<'_, String>) -> Result<String, UsageError> {
    match it.next() {
        Some(v) => Ok(v.clone()),
        None => Err(UsageError::new(format!("{flag} requires an argument"))),
    }
}

fn unknown(flag: &str) -> UsageError {
    UsageError::new(format!("unrecognised option '{flag}'"))
}

/// Parse a full argument vector (without argv[0]).
pub fn parse(args: &[String]) -> Result<Command, UsageError> {
    let Some(first) = args.first() else {
        return Err(UsageError::bare());
    };

    if first == "--version" {
        return Ok(Command::Version);
    }

    // `--quiet` is global: strip it wherever it appears, then parse the rest.
    let rest: Vec<String> = args[1..]
        .iter()
        .filter(|a| *a != "--quiet")
        .cloned()
        .collect();
    let quiet = args.iter().any(|a| a == "--quiet");

    match first.as_str() {
        "sort" => parse_sort(&rest, quiet),
        "merge" => parse_merge(&rest, quiet),
        "intersect" => parse_intersect(&rest, quiet),
        other => Err(unknown(other)),
    }
}

fn parse_sort(args: &[String], quiet: bool) -> Result<Command, UsageError> {
    let mut opts = SortOpts {
        input: String::new(),
        header: false,
        quiet,
    };
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "-i" => opts.input = value("-i", &mut it)?,
            "-header" => opts.header = true,
            // `sort -r` does not exist in bedtools; `-g`/`-faidx` and the
            // -size*/-chrThen* family are not in v1. All of them land here.
            other => return Err(unknown(other)),
        }
    }
    if opts.input.is_empty() {
        return Err(UsageError::new("sort requires -i"));
    }
    Ok(Command::Sort(opts))
}

fn parse_merge(args: &[String], quiet: bool) -> Result<Command, UsageError> {
    let mut opts = MergeOpts {
        input: String::new(),
        d: 0,
        same_strand: false,
        strand: None,
        header: false,
        quiet,
    };
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "-i" => opts.input = value("-i", &mut it)?,
            // -d takes its value unconditionally, so `-d -5` is a distance and
            // not a flag. Negative N requires that many bp of overlap.
            "-d" => {
                let raw = value("-d", &mut it)?;
                opts.d = raw
                    .parse::<i64>()
                    .map_err(|_| UsageError::new(format!("-d expects an integer, got '{raw}'")))?;
            }
            "-s" => opts.same_strand = true,
            "-S" => {
                let raw = value("-S", &mut it)?;
                opts.strand = match raw.as_str() {
                    "+" => Some('+'),
                    "-" => Some('-'),
                    _ => return Err(UsageError::new(format!("-S expects + or -, got '{raw}'"))),
                };
            }
            "-header" => opts.header = true,
            // -c/-o column aggregation is not in v1 (SPEC §4) and lands here.
            other => return Err(unknown(other)),
        }
    }
    if opts.input.is_empty() {
        return Err(UsageError::new("merge requires -i"));
    }
    Ok(Command::Merge(opts))
}

fn parse_intersect(args: &[String], quiet: bool) -> Result<Command, UsageError> {
    let mut opts = IntersectOpts {
        a: String::new(),
        b: String::new(),
        u: false,
        v: false,
        wa: false,
        header: false,
        quiet,
    };
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "-a" => opts.a = value("-a", &mut it)?,
            "-b" => opts.b = value("-b", &mut it)?,
            "-u" => opts.u = true,
            "-v" => opts.v = true,
            "-wa" => opts.wa = true,
            "-header" => opts.header = true,
            // -wb/-wo/-wao/-c/-loj/-f/-F/-r/-e/-s/-S are not in v1 (SPEC §4).
            other => return Err(unknown(other)),
        }
    }
    if opts.u && opts.v {
        return Err(UsageError::new("-u and -v are mutually exclusive"));
    }
    if opts.a.is_empty() {
        return Err(UsageError::new("intersect requires -a"));
    }
    if opts.b.is_empty() {
        return Err(UsageError::new("intersect requires -b"));
    }
    Ok(Command::Intersect(opts))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(args: &[&str]) -> Result<Command, UsageError> {
        let owned: Vec<String> = args.iter().map(|s| s.to_string()).collect();
        parse(&owned)
    }

    fn err(args: &[&str]) -> String {
        p(args).unwrap_err().to_string()
    }

    // --- the surface parses ---

    #[test]
    fn version_still_works() {
        assert_eq!(p(&["--version"]).unwrap(), Command::Version);
    }

    #[test]
    fn sort_flags_parse() {
        let Command::Sort(o) = p(&["sort", "-i", "a.bed", "-header"]).unwrap() else {
            panic!("not sort")
        };
        assert_eq!(o.input, "a.bed");
        assert!(o.header);
        assert!(!o.quiet);
    }

    #[test]
    fn stdin_is_a_dash() {
        let Command::Sort(o) = p(&["sort", "-i", "-"]).unwrap() else {
            panic!("not sort")
        };
        assert_eq!(o.input, "-");
    }

    #[test]
    fn merge_flags_parse() {
        let Command::Merge(o) = p(&["merge", "-i", "a.bed", "-d", "10", "-s", "-header"]).unwrap()
        else {
            panic!("not merge")
        };
        assert_eq!(o.d, 10);
        assert!(o.same_strand);
        assert!(o.header);
        assert_eq!(o.strand, None);
    }

    #[test]
    fn merge_d_accepts_a_negative_value() {
        // `-5` must be read as -d's argument, not as an unknown flag.
        let Command::Merge(o) = p(&["merge", "-i", "a.bed", "-d", "-5"]).unwrap() else {
            panic!("not merge")
        };
        assert_eq!(o.d, -5);
    }

    #[test]
    fn merge_strand_flag_parses_both_directions() {
        for (arg, want) in [("+", '+'), ("-", '-')] {
            let Command::Merge(o) = p(&["merge", "-i", "a.bed", "-S", arg]).unwrap() else {
                panic!("not merge")
            };
            assert_eq!(o.strand, Some(want));
        }
    }

    #[test]
    fn merge_defaults_to_d_zero() {
        let Command::Merge(o) = p(&["merge", "-i", "a.bed"]).unwrap() else {
            panic!("not merge")
        };
        assert_eq!(o.d, 0);
    }

    #[test]
    fn intersect_flags_parse() {
        let Command::Intersect(o) =
            p(&["intersect", "-a", "a.bed", "-b", "b.bed", "-wa", "-header"]).unwrap()
        else {
            panic!("not intersect")
        };
        assert_eq!((o.a.as_str(), o.b.as_str()), ("a.bed", "b.bed"));
        assert!(o.wa && o.header && !o.u && !o.v);
    }

    // --- usage errors: one test per row of the SPEC §7 table ---

    #[test]
    fn no_arguments_prints_usage_alone() {
        assert_eq!(err(&[]), USAGE);
    }

    #[test]
    fn unknown_subcommand_is_a_usage_error() {
        assert!(err(&["frobnicate"]).starts_with("mytools: unrecognised option 'frobnicate'\n"));
    }

    #[test]
    fn unknown_flag_is_a_usage_error() {
        // bedtools has no `sort -r`, so neither do we.
        assert!(
            err(&["sort", "-i", "a.bed", "-r"]).starts_with("mytools: unrecognised option '-r'\n")
        );
    }

    #[test]
    fn v1_exclusions_are_usage_errors() {
        for args in [
            vec!["sort", "-i", "a.bed", "-g", "genome.txt"],
            vec!["sort", "-i", "a.bed", "-faidx", "f.fai"],
            vec!["sort", "-i", "a.bed", "-sizeA"],
            vec!["merge", "-i", "a.bed", "-c", "4"],
            vec!["merge", "-i", "a.bed", "-o", "collapse"],
            vec!["intersect", "-a", "a.bed", "-b", "b.bed", "-wb"],
            vec!["intersect", "-a", "a.bed", "-b", "b.bed", "-wo"],
            vec!["intersect", "-a", "a.bed", "-b", "b.bed", "-wao"],
            vec!["intersect", "-a", "a.bed", "-b", "b.bed", "-c"],
            vec!["intersect", "-a", "a.bed", "-b", "b.bed", "-loj"],
            vec!["intersect", "-a", "a.bed", "-b", "b.bed", "-f", "0.5"],
            vec!["intersect", "-a", "a.bed", "-b", "b.bed", "-F", "0.5"],
            vec!["intersect", "-a", "a.bed", "-b", "b.bed", "-r"],
            vec!["intersect", "-a", "a.bed", "-b", "b.bed", "-e"],
            vec!["intersect", "-a", "a.bed", "-b", "b.bed", "-s"],
            vec!["intersect", "-a", "a.bed", "-b", "b.bed", "-S"],
        ] {
            let got = err(&args);
            assert!(
                got.starts_with("mytools: unrecognised option"),
                "{args:?} was accepted or misreported: {got}"
            );
        }
    }

    #[test]
    fn missing_flag_argument_is_a_usage_error() {
        assert!(err(&["sort", "-i"]).starts_with("mytools: -i requires an argument\n"));
        assert!(
            err(&["merge", "-i", "a.bed", "-d"]).starts_with("mytools: -d requires an argument\n")
        );
        assert!(err(&["intersect", "-a"]).starts_with("mytools: -a requires an argument\n"));
    }

    #[test]
    fn u_and_v_together_are_mutually_exclusive() {
        assert!(
            err(&["intersect", "-a", "a.bed", "-b", "b.bed", "-u", "-v"])
                .starts_with("mytools: -u and -v are mutually exclusive\n")
        );
    }

    #[test]
    fn every_usage_error_prints_the_usage_block() {
        for args in [vec!["sort", "-r"], vec!["nope"], vec!["intersect", "-a"]] {
            assert!(err(&args).ends_with(USAGE), "{args:?} did not print usage");
        }
    }

    // --- --quiet ---

    #[test]
    fn quiet_is_global_and_position_independent() {
        for args in [
            vec!["sort", "--quiet", "-i", "a.bed"],
            vec!["sort", "-i", "a.bed", "--quiet"],
        ] {
            let owned: Vec<String> = args.iter().map(|s| s.to_string()).collect();
            let Command::Sort(o) = parse(&owned).unwrap() else {
                panic!("not sort")
            };
            assert!(o.quiet, "{args:?} did not set quiet");
            assert_eq!(o.input, "a.bed");
        }
    }

    #[test]
    fn quiet_silences_nothing_else() {
        // It is not a flag-eater: the rest of the command still parses normally.
        let Command::Merge(o) =
            p(&["merge", "--quiet", "-i", "a.bed", "-d", "-5", "-S", "+"]).unwrap()
        else {
            panic!("not merge")
        };
        assert!(o.quiet);
        assert_eq!(o.d, -5);
        assert_eq!(o.strand, Some('+'));
    }
}
