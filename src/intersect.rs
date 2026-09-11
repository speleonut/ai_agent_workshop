//! `mytools intersect` — not implemented yet.

use crate::bed;
use crate::cli::IntersectOpts;

pub fn run(opts: &IntersectOpts) -> Result<(), bed::Error> {
    let _ = opts;
    Err(bed::Error::file("intersect", "not implemented yet"))
}
