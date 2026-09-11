//! `mytools merge` — not implemented yet.

use crate::bed;
use crate::cli::MergeOpts;

pub fn run(opts: &MergeOpts) -> Result<(), bed::Error> {
    let _ = opts;
    Err(bed::Error::file("merge", "not implemented yet"))
}
