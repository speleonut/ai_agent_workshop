//! `mytools intersect` — report features of A against features of B.
//!
//! `SPEC.md` §3 (semantics), §4 (flags), §5 (output), §6 (memory), §7 (warnings).
//!
//! Memory model (§6): **all of B** is held in memory — coordinates only, in a
//! bin index — and **A is streamed** one record at a time, never buffered. The
//! only per-A allocation is a reused hit buffer, bounded by the number of B
//! features one A record overlaps.
//!
//! Two pieces of oracle behaviour are encoded here rather than derived, both
//! verified against bedtools v2.31.1 (see [`inflate`] and [`bin_of`]).

use std::collections::HashMap;
use std::io::{self, BufWriter, Write};

use crate::bed::{self, Reader, Record};
use crate::cli::{self, IntersectOpts};

// ---------------------------------------------------------------- zero-length
//
/// Widen a zero-length feature by one base on each side.
///
/// **Oracle behaviour, not derivable from SPEC §3.** The strict `<` predicate
/// says a zero-length feature overlaps nothing (`500 < 500` is false), but
/// bedtools reports plenty of such pairs. Verified against bedtools v2.31.1
/// with single-line fixtures — every row below is a probe that was actually
/// run, `bare` being the region a flagless `intersect` printed:
///
/// ```text
///   A[500,500) B[400,600) -> chr1 500 500      A[500,500) B[501,600) -> none
///   A[500,500) B[500,600) -> chr1 500 500      A[500,500) B[400,499) -> none
///   A[500,500) B[400,500) -> chr1 500 500      A[500,500) B[502,502) -> none
///   A[500,500) B[500,500) -> chr1 500 500      A[500,500) B[498,498) -> none
///   A[500,500) B[501,501) -> chr1 500 500      A[100,200) B[ 99, 99) -> none
///   A[500,500) B[499,499) -> chr1 500 500      A[100,200) B[201,201) -> none
///   A[  0,  0) B[  0, 10) -> chr1 0 0          A[  0,  0) B[  1, 10) -> none
///   A[  0, 100) B[100,100) -> chr1 99 100
///   A[100,200) B[100,100) -> chr1 100 101
///   A[100,200) B[150,150) -> chr1 149 151
///   A[100,200) B[200,200) -> chr1 199 200
/// ```
///
/// Every one of those is explained by a single rule: a zero-length feature is
/// widened to `[start - 1, end + 1)` before anything is tested — on **both**
/// sides — and the clipped region is then `max(a.start, b_lo) .. min(a.end,
/// b_hi)` using A's *original* coordinates and B's *widened* ones. That is why
/// `A[0,100)` against a zero-length B at 100 prints `99 100`: the region is not
/// a sub-interval of A.
///
/// `crate::bed::overlaps` is untouched — it stays as SPEC §3 defines it, and a
/// unit test in `bed.rs` pins the divergence. The widening lives here, in the
/// one subcommand that has to match the oracle.
///
/// The `saturating_sub` clamp at 0 matters only on the A side (bedtools accepts
/// `A[0,0)`); a zero-length B at position 0 is rejected outright, see
/// [`Db::build`].
fn inflate(start: u64, end: u64) -> (u64, u64) {
    if start == end {
        (start.saturating_sub(1), end + 1)
    } else {
        (start, end)
    }
}

// ---------------------------------------------------------------- bin index
//
// bedtools indexes B in a UCSC-style hierarchy of bins: the finest level is
// 16 kb wide (`>> 14`), each coarser level is 8x wider (`>> 3`), and there are
// 8 levels. A feature lives in the finest bin that contains it whole.
//
// This is not just an optimisation over scanning all of B per A record: the
// order in which a query visits bins IS the row order of a bare `intersect`,
// and row order is part of the answer.
//
// Oracle behaviour, verified against bedtools v2.31.1: with A = `chr1 0 300000`
// and B (in file order) `140000-150000`, `200000-200010`, `0-20000`, `10-20`,
// `0-3000000`, `33000-34000`, bedtools emits them as `10-20`, `33000-34000`,
// `200000-200010`, `0-20000`, `140000-150000`, `0-3000000` — i.e. finest level
// first, ascending bin index within a level, B-file order within a bin.
const FIRST_SHIFT: u32 = 14;
const NEXT_SHIFT: u32 = 3;
const LEVELS: u32 = 8;

/// The single bin holding `[lo, hi)`: the finest level whose bin contains the
/// whole feature. `hi - 1` is the last base covered, BED being half-open.
fn bin_of(lo: u64, hi: u64) -> (u32, u64) {
    let last = hi.saturating_sub(1);
    for level in 0..LEVELS {
        let shift = FIRST_SHIFT + NEXT_SHIFT * level;
        if (lo >> shift) == (last >> shift) {
            return (level, lo >> shift);
        }
    }
    // Wider than the coarsest level (32 Gb) — bedtools drops such a feature in
    // bin 0, which the coarsest level of every query visits.
    (LEVELS - 1, 0)
}

/// Visit every bin a query over `[lo, hi)` must look in, in bedtools' order:
/// finest level first, ascending index. `f` returns `false` to stop early.
fn visit_bins(lo: u64, hi: u64, mut f: impl FnMut(u32, u64) -> bool) {
    let mut start = lo >> FIRST_SHIFT;
    let mut end = hi.saturating_sub(1) >> FIRST_SHIFT;
    for level in 0..LEVELS {
        for idx in start..=end {
            if !f(level, idx) {
                return;
            }
        }
        start >>= NEXT_SHIFT;
        end >>= NEXT_SHIFT;
    }
}

/// All of B: widened coordinates in file order, plus a per-chromosome bin index
/// into them. Lines are **not** kept — no v1 flag reports a B record, so B costs
/// 16 bytes per feature plus the index.
struct Db {
    /// chrom -> (level, bin) -> B ids, in B-file order.
    bins: HashMap<String, HashMap<(u32, u64), Vec<u32>>>,
    /// Widened `[lo, hi)` per B record, in B-file order.
    ivals: Vec<(u64, u64)>,
}

impl Db {
    fn build(recs: Vec<Record>, name: &str) -> Result<Db, bed::Error> {
        let mut db = Db {
            bins: HashMap::new(),
            ivals: Vec::with_capacity(recs.len()),
        };
        for rec in recs {
            // Oracle behaviour: widening a zero-length B feature at position 0
            // gives it a start of -1, which bedtools cannot bin. Verified
            // against bedtools v2.31.1 — with `chr1 0 0` anywhere in B it
            // prints "illegal bin number -1 ... Unable to add record to tree"
            // and exits 1 with zero bytes on stdout, for bare, -u and -v alike.
            // We match the exit code and the empty stdout; the wording is ours.
            if rec.start == rec.end && rec.start == 0 {
                return Err(bed::Error::file(
                    name,
                    "zero-length feature at position 0 cannot be indexed",
                ));
            }
            let (lo, hi) = inflate(rec.start, rec.end);
            let id = db.ivals.len() as u32;
            db.ivals.push((lo, hi));
            db.bins
                .entry(rec.chrom)
                .or_default()
                .entry(bin_of(lo, hi))
                .or_default()
                .push(id);
        }
        Ok(db)
    }

    /// Every B feature `a` overlaps, as the clipped region to print, in the
    /// oracle's row order. `out` is cleared first and reused across A records
    /// so that streaming A costs no allocation per record.
    fn hits_into(&self, a: &Record, out: &mut Vec<(u64, u64)>) {
        out.clear();
        let Some(bins) = self.bins.get(&a.chrom) else {
            return; // a chromosome in A but not in B, and vice versa: no hits.
        };
        let (a_lo, a_hi) = inflate(a.start, a.end);
        visit_bins(a_lo, a_hi, |level, idx| {
            if let Some(ids) = bins.get(&(level, idx)) {
                for &id in ids {
                    let (b_lo, b_hi) = self.ivals[id as usize];
                    // The one predicate (SPEC §3), applied to widened
                    // coordinates. Clipping uses A's original start/end, so a
                    // zero-length A prints itself.
                    if bed::overlaps(a_lo, a_hi, b_lo, b_hi) {
                        out.push((a.start.max(b_lo), a.end.min(b_hi)));
                    }
                }
            }
            true
        });
    }

    /// Does `a` overlap anything in B? Stops at the first hit — all `-u` and
    /// `-v` need.
    fn any_hit(&self, a: &Record) -> bool {
        let Some(bins) = self.bins.get(&a.chrom) else {
            return false;
        };
        let (a_lo, a_hi) = inflate(a.start, a.end);
        let mut found = false;
        visit_bins(a_lo, a_hi, |level, idx| {
            if let Some(ids) = bins.get(&(level, idx)) {
                for &id in ids {
                    let (b_lo, b_hi) = self.ivals[id as usize];
                    if bed::overlaps(a_lo, a_hi, b_lo, b_hi) {
                        found = true;
                        return false;
                    }
                }
            }
            true
        });
        found
    }
}

// ---------------------------------------------------------------- the run
//
/// A data error or a write failure. Two error types meet in the output loop;
/// `run` maps both back to what `main` expects.
enum Fail {
    Data(bed::Error),
    Io(io::Error),
}

impl From<bed::Error> for Fail {
    fn from(e: bed::Error) -> Fail {
        Fail::Data(e)
    }
}

impl From<io::Error> for Fail {
    fn from(e: io::Error) -> Fail {
        Fail::Io(e)
    }
}

pub fn run(opts: &IntersectOpts) -> Result<(), bed::Error> {
    // B first and whole (SPEC §6). collect_records fails on the first bad line.
    let mut b_reader = Reader::open(&opts.b)?;
    let b_recs = b_reader.collect_records()?;
    // Width is taken from B's first record: mixed widths interoperate, but they
    // warn (SPEC §2, §7).
    let b_width = b_recs.first().map(|r| r.ncols());
    let db = Db::build(b_recs, b_reader.name())?;

    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout.lock());

    match emit(opts, &db, b_width, &mut out) {
        Ok(()) => match out.flush() {
            Ok(()) => Ok(()),
            // A closed downstream (`| head`) is not our failure.
            Err(e) if e.kind() == io::ErrorKind::BrokenPipe => Ok(()),
            Err(_) => Err(bed::Error::file("stdout", "cannot be written")),
        },
        Err(Fail::Data(e)) => Err(e),
        Err(Fail::Io(e)) if e.kind() == io::ErrorKind::BrokenPipe => Ok(()),
        Err(Fail::Io(_)) => Err(bed::Error::file("stdout", "cannot be written")),
    }
}

/// Stream A against B and write the results. One A record is live at a time.
fn emit(
    opts: &IntersectOpts,
    db: &Db,
    b_width: Option<usize>,
    out: &mut impl Write,
) -> Result<(), Fail> {
    let mut a_reader = Reader::open(&opts.a)?;

    // Pull the first record, then write the headers seen up to that point.
    //
    // The snapshot point is the whole rule. Oracle behaviour, verified against
    // bedtools v2.31.1: only the header block ABOVE the first data record is
    // replayed. Given `chr2 5 6`, `#mid`, `chr1 1 2`, `track x`, `chr1 0 1`,
    // `-header` prints the records and no header line at all. A file with
    // headers and no data records does still replay them — which is what the
    // `first == None` path gives us, the reader having reached EOF by then.
    // Headers also appear when nothing matches, verified likewise.
    //
    // `Reader::headers()` accumulates mid-file lines too, so taking it at the
    // end instead would silently diverge; data/header.bed has a clean leading
    // block and would not catch it.
    let first = a_reader.next().transpose()?;
    if opts.header {
        let headers: Vec<String> = a_reader.headers().to_vec();
        for h in headers {
            writeln!(out, "{h}")?;
        }
    }

    let mut hits: Vec<(u64, u64)> = Vec::new();
    let mut warned = false;

    for item in first.into_iter().map(Ok).chain(a_reader) {
        let rec = item?;

        // Mixed widths interoperate; they warn once, on stderr, and leave the
        // exit code and stdout alone (SPEC §2, §7, §8.2).
        if !warned
            && let Some(w) = b_width
            && w != rec.ncols()
        {
            cli::warn(
                opts.quiet,
                &format!(
                    "input files have different BED widths: A has {} columns, B has {w}",
                    rec.ncols()
                ),
            );
            warned = true;
        }

        if opts.v {
            if !db.any_hit(&rec) {
                writeln!(out, "{}", rec.line)?;
            }
        } else if opts.u {
            if db.any_hit(&rec) {
                writeln!(out, "{}", rec.line)?;
            }
        } else {
            db.hits_into(&rec, &mut hits);
            if opts.wa {
                // One original A record per overlapping pair, not one per A.
                for _ in 0..hits.len() {
                    writeln!(out, "{}", rec.line)?;
                }
            } else {
                for &(start, end) in &hits {
                    write_clipped(out, &rec, start, end)?;
                }
            }
        }
    }
    Ok(())
}

/// The clipped region, carrying A's columns past the third unchanged.
fn write_clipped(out: &mut impl Write, rec: &Record, start: u64, end: u64) -> io::Result<()> {
    write!(out, "{}\t{}\t{}", rec.chrom, start, end)?;
    for col in rec.line.split('\t').skip(3) {
        write!(out, "\t{col}")?;
    }
    out.write_all(b"\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(chrom: &str, start: u64, end: u64, rest: &str) -> Record {
        let line = if rest.is_empty() {
            format!("{chrom}\t{start}\t{end}")
        } else {
            format!("{chrom}\t{start}\t{end}\t{rest}")
        };
        Record {
            chrom: chrom.to_string(),
            start,
            end,
            line,
        }
    }

    /// Build a Db from `(chrom, start, end)` triples in B-file order.
    fn db(items: &[(&str, u64, u64)]) -> Db {
        let recs: Vec<Record> = items.iter().map(|&(c, s, e)| rec(c, s, e, "")).collect();
        Db::build(recs, "b.bed").expect("B indexed")
    }

    /// The clipped regions a bare `intersect` would print for one A record.
    fn hits(db: &Db, a: &Record) -> Vec<(u64, u64)> {
        let mut out = Vec::new();
        db.hits_into(a, &mut out);
        out
    }

    /// One A interval against one B interval: the region printed, if any.
    fn pair(a: (u64, u64), b: (u64, u64)) -> Option<(u64, u64)> {
        let d = db(&[("chr1", b.0, b.1)]);
        let got = hits(&d, &rec("chr1", a.0, a.1, ""));
        assert!(got.len() <= 1);
        got.first().copied()
    }

    // --- the overlap predicate, through intersect's own path (SPEC §3) ---

    #[test]
    fn bookended_intervals_are_not_a_hit() {
        // a.end == b.start: 100..199 and 200..299 share no base.
        assert_eq!(pair((100, 200), (200, 300)), None);
        assert_eq!(pair((200, 300), (100, 200)), None);
    }

    #[test]
    fn one_shared_base_pair_is_a_hit() {
        // Base 199 is in both, and one bp is the minimum overlap (SPEC §3).
        assert_eq!(pair((100, 200), (199, 300)), Some((199, 200)));
    }

    #[test]
    fn fully_nested_intervals_hit_and_clip_to_the_inner_one() {
        assert_eq!(pair((100, 500), (200, 300)), Some((200, 300)));
        assert_eq!(pair((200, 300), (100, 500)), Some((200, 300)));
    }

    #[test]
    fn identical_coordinates_hit_and_clip_to_themselves() {
        assert_eq!(pair((150, 250), (150, 250)), Some((150, 250)));
    }

    #[test]
    fn interval_at_position_zero() {
        assert_eq!(pair((0, 100), (0, 50)), Some((0, 50)));
        // Bookended at 100, and 0 is not a special case for that.
        assert_eq!(pair((0, 100), (100, 200)), None);
    }

    // --- clipping: the region printed is the intersection, not either input ---

    #[test]
    fn the_emitted_region_is_the_intersection_of_the_two() {
        // Neither (100,200) nor (150,400) is the answer: (150,200) is.
        assert_eq!(pair((100, 200), (150, 400)), Some((150, 200)));
        assert_eq!(pair((150, 400), (100, 200)), Some((150, 200)));
    }

    #[test]
    fn clipping_carries_a_columns_past_the_third_unchanged() {
        let d = db(&[("chr1", 150, 400)]);
        let a = rec("chr1", 100, 200, "a02\t20\t-");
        let mut buf = Vec::new();
        let h = hits(&d, &a);
        write_clipped(&mut buf, &a, h[0].0, h[0].1).unwrap();
        assert_eq!(
            String::from_utf8(buf).unwrap(),
            "chr1\t150\t200\ta02\t20\t-\n"
        );
    }

    #[test]
    fn a_bed3_record_clips_to_three_columns() {
        let a = rec("chr1", 100, 200, "");
        let mut buf = Vec::new();
        write_clipped(&mut buf, &a, 150, 200).unwrap();
        assert_eq!(String::from_utf8(buf).unwrap(), "chr1\t150\t200\n");
    }

    // --- zero-length intervals: oracle behaviour, verified, see `inflate` ---

    #[test]
    fn zero_length_a_is_reported_at_its_own_coordinates() {
        // Oracle behaviour (bedtools v2.31.1): `-a chr1 500 500` against each
        // of these Bs prints `chr1 500 500` — the region is A itself, never the
        // widened interval the test used.
        for b in [(400, 600), (500, 600), (400, 500), (499, 501)] {
            assert_eq!(pair((500, 500), b), Some((500, 500)), "B {b:?}");
        }
    }

    #[test]
    fn zero_length_a_reaches_exactly_one_base_either_side() {
        // Oracle behaviour: bedtools reports neither of these. One base further
        // out in each direction and the widened interval no longer touches.
        assert_eq!(pair((500, 500), (501, 600)), None);
        assert_eq!(pair((500, 500), (400, 499)), None);
    }

    #[test]
    fn two_zero_length_features_meet_where_the_strict_predicate_says_they_cannot() {
        // Oracle behaviour, the divergence CLAUDE.md warns about and the case
        // data/a.bed a07 vs data/b.bed b07 hits: bedtools REPORTS
        // `chr1 500 500` against `chr1 500 500`, while bed::overlaps says no
        // because `500 < 500` is false. Both sides are widened, so they meet.
        assert!(!bed::overlaps(500, 500, 500, 500));
        assert_eq!(pair((500, 500), (500, 500)), Some((500, 500)));
        // ...and one base apart still meets, two bases apart does not.
        assert_eq!(pair((500, 500), (501, 501)), Some((500, 500)));
        assert_eq!(pair((500, 500), (499, 499)), Some((500, 500)));
        assert_eq!(pair((500, 500), (502, 502)), None);
        assert_eq!(pair((500, 500), (498, 498)), None);
    }

    #[test]
    fn zero_length_b_widens_the_region_beyond_a() {
        // Oracle behaviour: the printed region is NOT a sub-interval of A when
        // B is zero-length. Verified against bedtools v2.31.1.
        assert_eq!(pair((0, 100), (100, 100)), Some((99, 100))); // a01 vs b02
        assert_eq!(pair((100, 200), (100, 100)), Some((100, 101))); // a02 vs b02
        assert_eq!(pair((100, 200), (200, 200)), Some((199, 200)));
        assert_eq!(pair((100, 200), (150, 150)), Some((149, 151)));
        // Two bases clear of A on either side: no hit.
        assert_eq!(pair((100, 200), (99, 99)), None);
        assert_eq!(pair((100, 200), (201, 201)), None);
    }

    #[test]
    fn zero_length_a_at_position_zero_does_not_underflow() {
        // Oracle behaviour: bedtools reports `chr1 0 0` against `chr1 0 10` and
        // against `chr1 0 1`, but not against `chr1 1 10`. Widening A's start
        // clamps at 0 rather than going negative, which changes no answer: no
        // legal B can end at or below 0.
        assert_eq!(pair((0, 0), (0, 10)), Some((0, 0)));
        assert_eq!(pair((0, 0), (0, 1)), Some((0, 0)));
        assert_eq!(pair((0, 0), (1, 10)), None);
        assert_eq!(pair((1, 1), (0, 1)), Some((1, 1)));
    }

    #[test]
    fn zero_length_b_at_position_zero_is_an_error() {
        // Oracle behaviour: bedtools cannot bin a feature widened to start -1.
        // With `chr1 0 0` in B it prints "illegal bin number -1" and exits 1
        // with zero bytes on stdout, for bare, -u and -v alike. We match the
        // exit code (a bed::Error is exit 1) and the empty stdout.
        let recs = vec![rec("chr1", 50, 60, ""), rec("chr1", 0, 0, "")];
        let e = match Db::build(recs, "b.bed") {
            Err(e) => e,
            Ok(_) => panic!("indexed a zero-length B feature at position 0"),
        };
        assert_eq!(e.msg, "zero-length feature at position 0 cannot be indexed");
        assert_eq!(
            e.to_string(),
            "mytools: b.bed: zero-length feature at position 0 cannot be indexed"
        );
    }

    // --- chromosomes, and B-only chromosomes ---

    #[test]
    fn a_chromosome_present_only_in_b_is_no_hit_and_no_panic() {
        // data/b.bed has chr3; data/a.bed does not.
        let d = db(&[("chr3", 100, 200), ("chr1", 100, 200)]);
        assert_eq!(hits(&d, &rec("chr1", 150, 250, "")), vec![(150, 200)]);
        assert_eq!(hits(&d, &rec("chrZ", 150, 250, "")), vec![]);
        assert!(!d.any_hit(&rec("chrZ", 150, 250, "")));
    }

    #[test]
    fn identical_coordinates_on_different_chromosomes_do_not_meet() {
        let d = db(&[("chr2", 100, 200)]);
        assert_eq!(hits(&d, &rec("chr1", 100, 200, "")), vec![]);
    }

    // --- row order is part of the answer ---

    #[test]
    fn hits_come_back_in_the_oracles_bin_order() {
        // Oracle behaviour, verified against bedtools v2.31.1 with exactly this
        // input: finest bin level first, ascending bin index within a level,
        // B-file order within a bin. See the `bin_of` comment for the probe.
        let d = db(&[
            ("chr1", 140_000, 150_000), // level 1, index 1
            ("chr1", 200_000, 200_010), // level 0, index 12
            ("chr1", 0, 20_000),        // level 1, index 0
            ("chr1", 10, 20),           // level 0, index 0
            ("chr1", 0, 3_000_000),     // level 2
            ("chr1", 33_000, 34_000),   // level 0, index 2
        ]);
        let got = hits(&d, &rec("chr1", 0, 300_000, ""));
        assert_eq!(
            got,
            vec![
                (10, 20),
                (33_000, 34_000),
                (200_000, 200_010),
                (0, 20_000),
                (140_000, 150_000),
                (0, 300_000),
            ]
        );
    }

    #[test]
    fn hits_within_one_bin_keep_b_file_order() {
        // Everything below 16 kb shares the finest bin, so file order decides —
        // which is why data/a.bed a02 reports b03 (line 1) before b02 (line 3).
        let d = db(&[("chr1", 180, 220), ("chr1", 0, 50), ("chr1", 100, 100)]);
        let got = hits(&d, &rec("chr1", 100, 200, ""));
        assert_eq!(got, vec![(180, 200), (100, 101)]);
    }

    #[test]
    fn every_overlapping_b_feature_produces_its_own_row() {
        // Bare and -wa emit one row per pair, not one per A record.
        let d = db(&[("chr1", 320, 350), ("chr1", 340, 360)]);
        assert_eq!(
            hits(&d, &rec("chr1", 300, 400, "")), // data/a.bed a05
            vec![(320, 350), (340, 360)]
        );
    }

    // --- the bin index itself ---

    #[test]
    fn bin_of_picks_the_finest_level_that_holds_the_feature_whole() {
        assert_eq!(bin_of(10, 20), (0, 0));
        assert_eq!(bin_of(16_384, 16_400), (0, 1));
        // Straddles two 16 kb bins, so it moves up a level.
        assert_eq!(bin_of(16_000, 17_000), (1, 0));
        // Wider than the coarsest level: bin 0 of the top level, which every
        // query visits.
        assert_eq!(bin_of(0, u64::MAX), (LEVELS - 1, 0));
    }

    #[test]
    fn a_query_visits_the_bin_its_target_lives_in() {
        for (lo, hi) in [
            (0u64, 1u64),
            (10, 20),
            (16_000, 17_000),
            (5_000_000, 9_000_000),
        ] {
            let want = bin_of(lo, hi);
            let mut seen = false;
            visit_bins(lo, hi, |level, idx| {
                if (level, idx) == want {
                    seen = true;
                }
                true
            });
            assert!(seen, "query [{lo},{hi}) never visited its own bin {want:?}");
        }
    }

    #[test]
    fn visit_bins_stops_when_the_callback_says_so() {
        // -u and -v lean on this: the first hit ends the search.
        let mut count = 0;
        visit_bins(0, 300_000, |_, _| {
            count += 1;
            false
        });
        assert_eq!(count, 1);
    }

    // --- -u / -v agree with the hit list ---

    #[test]
    fn any_hit_matches_whether_the_hit_list_is_empty() {
        let d = db(&[("chr1", 200, 300), ("chr1", 500, 500)]);
        for a in [
            rec("chr1", 100, 200, ""), // bookended: no
            rec("chr1", 199, 201, ""), // one bp: yes
            rec("chr1", 500, 500, ""), // zero-length pair: yes
            rec("chr1", 900, 950, ""), // nowhere near: no
            rec("chr2", 200, 300, ""), // wrong chromosome: no
        ] {
            assert_eq!(d.any_hit(&a), !hits(&d, &a).is_empty(), "{}", a.line);
        }
    }

    #[test]
    fn empty_b_means_no_hits_for_anything() {
        let d = db(&[]);
        assert_eq!(hits(&d, &rec("chr1", 0, 100, "")), vec![]);
        assert!(!d.any_hit(&rec("chr1", 0, 100, "")));
    }

    // --- the whole run, through emit: flags, output width, headers ---

    /// Run `emit` over an A file written to a temp path, against `b`, and
    /// return stdout as a string.
    fn emit_to_string(
        a_body: &str,
        b: &[(&str, u64, u64)],
        set: impl Fn(&mut IntersectOpts),
    ) -> String {
        use std::io::Write as _;
        use std::sync::atomic::{AtomicU32, Ordering};
        static SEQ: AtomicU32 = AtomicU32::new(0);

        let mut path = std::env::temp_dir();
        path.push(format!(
            "mytools-intersect-test-{}-{}.bed",
            std::process::id(),
            SEQ.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::File::create(&path)
            .unwrap()
            .write_all(a_body.as_bytes())
            .unwrap();

        let mut opts = IntersectOpts {
            a: path.to_str().unwrap().to_string(),
            b: String::new(), // emit never opens B; the Db is already built.
            u: false,
            v: false,
            wa: false,
            header: false,
            quiet: true,
        };
        set(&mut opts);

        let d = db(b);
        let mut buf: Vec<u8> = Vec::new();
        let outcome = emit(&opts, &d, None, &mut buf);
        std::fs::remove_file(&path).ok();
        assert!(matches!(outcome, Ok(())), "emit failed");
        String::from_utf8(buf).unwrap()
    }

    /// Two A records, the first overlapping both B features and the second
    /// overlapping neither.
    const A_BODY: &str = "chr1\t300\t400\ta05\t40\t+\nchr2\t100\t120\ta14\t55\t-\n";
    const B_ITEMS: [(&str, u64, u64); 2] = [("chr1", 320, 350), ("chr1", 340, 360)];

    #[test]
    fn bare_emits_the_clipped_region_once_per_pair() {
        assert_eq!(
            emit_to_string(A_BODY, &B_ITEMS, |_| {}),
            "chr1\t320\t350\ta05\t40\t+\nchr1\t340\t360\ta05\t40\t+\n"
        );
    }

    #[test]
    fn wa_emits_the_original_a_record_once_per_pair() {
        assert_eq!(
            emit_to_string(A_BODY, &B_ITEMS, |o| o.wa = true),
            "chr1\t300\t400\ta05\t40\t+\nchr1\t300\t400\ta05\t40\t+\n"
        );
    }

    #[test]
    fn u_emits_each_overlapping_a_record_once_at_full_width() {
        assert_eq!(
            emit_to_string(A_BODY, &B_ITEMS, |o| o.u = true),
            "chr1\t300\t400\ta05\t40\t+\n"
        );
    }

    #[test]
    fn v_emits_the_non_overlapping_a_records_at_full_width() {
        assert_eq!(
            emit_to_string(A_BODY, &B_ITEMS, |o| o.v = true),
            "chr2\t100\t120\ta14\t55\t-\n"
        );
    }

    #[test]
    fn no_matches_and_empty_input_print_zero_bytes() {
        // SPEC §5: an empty result is zero bytes and exit 0 (emit returning Ok).
        assert_eq!(emit_to_string(A_BODY, &[("chr9", 10, 20)], |_| {}), "");
        assert_eq!(emit_to_string("", &B_ITEMS, |_| {}), "");
        assert_eq!(emit_to_string("", &B_ITEMS, |o| o.v = true), "");
    }

    #[test]
    fn header_replays_only_the_block_above_the_first_record() {
        // Oracle behaviour, verified against bedtools v2.31.1: with headers
        // interleaved after the first data record, `-header` prints no header
        // line at all. Reader::headers() accumulates the mid-file ones, so the
        // rule is the snapshot point, not the list.
        let body = "chr1\t5\t6\n#mid\nchr1\t1\t2\ntrack x\nchr1\t0\t1\n";
        let b = [("chr1", 0, 10)];
        assert_eq!(
            emit_to_string(body, &b, |o| o.header = true),
            "chr1\t5\t6\nchr1\t1\t2\nchr1\t0\t1\n"
        );
    }

    #[test]
    fn header_replays_the_leading_block_even_with_no_results() {
        let body = "#c\ntrack t\nchr1\t500\t600\n";
        let b = [("chr1", 0, 10)];
        assert_eq!(
            emit_to_string(body, &b, |o| o.header = true),
            "#c\ntrack t\n"
        );
        // Without -header they never appear (SPEC §5).
        assert_eq!(emit_to_string(body, &b, |_| {}), "");
    }

    #[test]
    fn a_file_of_headers_alone_still_replays_them() {
        // Oracle behaviour: no data record ever arrives, so the snapshot is the
        // whole file. bedtools prints both lines.
        let body = "#only\ntrack y\n";
        assert_eq!(
            emit_to_string(body, &[("chr1", 0, 10)], |o| o.header = true),
            "#only\ntrack y\n"
        );
    }
}
