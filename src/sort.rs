//! `mytools sort` — not implemented yet.

use crate::bed;
use crate::cli::SortOpts;

pub fn run(opts: &SortOpts) -> Result<(), bed::Error> {
    let _ = opts;
    Err(bed::Error::file("sort", "not implemented yet"))
}
