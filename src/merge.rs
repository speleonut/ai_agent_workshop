//! `mytools merge` — collapse overlapping and nearby features on the same
//! chromosome into single BED3 records (`SPEC.md` §3, §4, §5).
//!
//! **Memory (`SPEC.md` §6): this streams.** Input is required to be sorted, so
//! the only state is the currently open merge candidates — exactly one cluster
//! per strand key (one without `-s`, at most two with it) plus the set of
//! chromosome names seen, which the sortedness check needs. Nothing scales with
//! the number of records.
//!
//! Several rules here are bedtools' and not anybody's first guess. They are
//! marked `ORACLE:` and were verified against bedtools v2.31.1; the unit tests
//! at the bottom pin each one. Do not "simplify" them back to what seems right.

use std::collections::HashSet;
use std::collections::VecDeque;
use std::io::{self, Write};

use crate::bed::{self, Record};
use crate::cli::{self, MergeOpts};

/// The interval a record occupies *for merging purposes*, as signed
/// coordinates.
///
/// ORACLE: a zero-length record (`start == end`) is inflated by one base in
/// **both** directions — `chr1 500 500` merges as if it were `chr1 499 501`.
/// Verified against bedtools v2.31.1 on throwaway fixtures:
///
/// - `chr1 500 500` + `chr1 500 600`  -> `chr1 499 600` (also the `a.bed` case)
/// - `chr1 500 500` + `chr1 500 500`  -> `chr1 499 501`
/// - `chr1 500 500` + `chr1 501 600`  -> `chr1 499 600` (bookended at 501)
/// - `chr1 500 500` + `chr1 502 600`  -> two records, unmerged
/// - `chr1 100 600` + `chr1 601 601`  -> `chr1 100 602`
///
/// Signed, because the left edge can go below zero: `chr1 0 0` + `chr1 0 100`
/// prints `chr1 -1 100`, negative start and all. bedtools does not clamp it and
/// neither do we.
fn span(rec: &Record) -> (i64, i64) {
    let start = rec.start as i64;
    let end = rec.end as i64;
    if start == end {
        (start - 1, end + 1)
    } else {
        (start, end)
    }
}

/// Column 6 narrowed to the two strands `-s`/`-S` recognise.
///
/// ORACLE: under a strand flag, a record whose strand is neither `+` nor `-`
/// (`.`, say) is dropped entirely rather than forming its own group —
/// `bedtools merge -s` over two overlapping `.` records prints nothing.
/// `None` here therefore means both "no column 6" and "not a strand".
fn strand_of(rec: &Record) -> Option<char> {
    match rec.strand() {
        Some("+") => Some('+'),
        Some("-") => Some('-'),
        _ => None,
    }
}

/// One run of merge candidates: the records collapsed so far into one output
/// record, plus enough to print it.
struct Cluster {
    chrom: String,
    /// Strand group this cluster belongs to: `None` without `-s`/`-S`.
    key: Option<char>,
    /// Union of the [`span`]s added so far — what gets printed for a cluster
    /// of two or more.
    start: i64,
    end: i64,
    /// The first record's *original* coordinates, printed when the cluster
    /// turns out to hold only that one record. See [`Cluster::bounds`].
    solo: (i64, i64),
    n: u64,
    /// A closed cluster takes no more records but may still be waiting behind
    /// an older one for its turn to be printed.
    open: bool,
}

impl Cluster {
    fn new(key: Option<char>, rec: &Record) -> Cluster {
        let (start, end) = span(rec);
        Cluster {
            chrom: rec.chrom.clone(),
            key,
            start,
            end,
            solo: (rec.start as i64, rec.end as i64),
            n: 1,
            open: true,
        }
    }

    /// Does `rec` join this cluster? The gap between them is
    /// `rec.start - cluster.end`, negative when they overlap, so one signed
    /// comparison covers every case of `-d`: `0` merges overlapping and
    /// bookended features (gap `0`), a positive `N` merges across a gap of up
    /// to `N`, and a negative `N` demands `N` bp of overlap. Boundaries are
    /// inclusive both ways — verified against the oracle: gap `10` merges under
    /// `-d 10` and gap `11` does not; 5 bp of overlap merges under `-d -5` and
    /// 4 bp does not.
    fn accepts(&self, rec: &Record, d: i64) -> bool {
        span(rec).0 - self.end <= d
    }

    fn add(&mut self, rec: &Record) {
        let (start, end) = span(rec);
        // `start` can still fall: input is sorted on the *original* start, and
        // a zero-length record's span reaches one base further left than that.
        self.start = self.start.min(start);
        self.end = self.end.max(end);
        self.n += 1;
    }

    /// ORACLE: a cluster of exactly one record prints that record's original
    /// coordinates, *not* its inflated [`span`] — `bedtools merge` on a lone
    /// `chr1 500 500` prints `chr1 500 500`, while two of them print
    /// `chr1 499 501`. The inflation is only ever visible once something
    /// actually merged. For every record that is not zero-length the two are
    /// identical, so this only bites on the fixtures' zero-length features.
    fn bounds(&self) -> (i64, i64) {
        if self.n == 1 {
            self.solo
        } else {
            (self.start, self.end)
        }
    }
}

/// stdout, with the `-header` lines held back until we know what precedes the
/// first data record.
struct Sink<'a, W: Write> {
    out: &'a mut W,
    header: bool,
    /// Header lines as they stood when the first data record was read. Taken
    /// and written once, by [`Sink::begin`].
    pending: Option<Vec<String>>,
    started: bool,
}

impl<'a, W: Write> Sink<'a, W> {
    fn new(out: &'a mut W, header: bool) -> Sink<'a, W> {
        Sink {
            out,
            header,
            pending: None,
            started: false,
        }
    }

    /// Record the header lines to reproduce. ORACLE: only the ones ahead of the
    /// first data record count — `bedtools merge -header` on a file with a `#`
    /// line in the middle prints the leading one and swallows the other. So
    /// this is called once, the moment the first record arrives (or at end of
    /// input, if there were no records at all).
    fn set_headers(&mut self, lines: Vec<String>) {
        if !self.started && self.pending.is_none() {
            self.pending = Some(lines);
        }
    }

    fn begin(&mut self) -> Result<(), bed::Error> {
        if !self.started {
            self.started = true;
            let lines = self.pending.take().unwrap_or_default();
            if self.header {
                for line in lines {
                    wr(writeln!(self.out, "{line}"))?;
                }
            }
        }
        Ok(())
    }

    /// Headers come out even when the result is empty: `bedtools merge -header`
    /// over input that merges to nothing still prints them.
    fn finish(&mut self) -> Result<(), bed::Error> {
        self.begin()?;
        wr(self.out.flush())
    }

    fn emit(&mut self, cluster: &Cluster) -> Result<(), bed::Error> {
        self.begin()?;
        let (start, end) = cluster.bounds();
        wr(writeln!(self.out, "{}\t{}\t{}", cluster.chrom, start, end))
    }
}

/// stdout is data; a failure to write it is not something we can report *on*
/// stdout, so it becomes a data error on stderr like any other (`SPEC.md` §7).
fn wr(r: io::Result<()>) -> Result<(), bed::Error> {
    r.map_err(|_| bed::Error::file("<stdout>", "cannot be written"))
}

// High-water mark of the cluster queue, so the streaming claim in `SPEC.md` §6
// is something a test can assert rather than something this file merely says.
// Thread-local, because `cargo test` runs each test in its own thread.
#[cfg(test)]
thread_local! {
    static PEAK_QUEUE: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

#[cfg(test)]
fn note_queue_len(len: usize) {
    PEAK_QUEUE.with(|p| p.set(p.get().max(len)));
}

#[cfg(not(test))]
fn note_queue_len(_len: usize) {}

pub fn run(opts: &MergeOpts) -> Result<(), bed::Error> {
    let mut reader = bed::Reader::open(&opts.input)?;
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());
    merge_into(&mut reader, opts, &mut out)
}

/// The whole of `merge`, with the output sink injected so unit tests can run it
/// without a subprocess.
fn merge_into<W: Write>(
    reader: &mut bed::Reader,
    opts: &MergeOpts,
    out: &mut W,
) -> Result<(), bed::Error> {
    let name = reader.name().to_string();
    let stranded = opts.same_strand || opts.strand.is_some();

    let mut sink = Sink::new(out, opts.header);
    // At most one open cluster per strand key, so at most two entries under
    // `-s`. Ordered by when each cluster was *opened*, which is the order
    // bedtools prints them in — see the loop below.
    let mut queue: VecDeque<Cluster> = VecDeque::new();
    let mut seen_chroms: HashSet<String> = HashSet::new();
    let mut prev: Option<(String, u64)> = None;
    let mut records: u64 = 0;
    let mut warned = false;

    while let Some(item) = reader.next() {
        let rec = item?;
        records += 1;
        if records == 1 {
            sink.set_headers(reader.headers().to_vec());
        }

        // The true file line of the record just taken, blank and header lines
        // included, so the message points where `sed -n '<n>p'` would.
        let lineno = reader.lineno();
        let unsorted = || bed::Error::at(&name, lineno, "input is not sorted");

        match &prev {
            // Within a chromosome, starts must not go backwards. Equal starts
            // are fine and `end` is never compared.
            Some((chrom, start)) if *chrom == rec.chrom => {
                if rec.start < *start {
                    return Err(unsorted());
                }
            }
            // A new chromosome. ORACLE: bedtools does *not* require
            // chromosomes in lexicographic order — `chr2` before `chr1` is
            // accepted — only that each chromosome's records are contiguous.
            // Coming back to one already finished is the error.
            _ => {
                if !seen_chroms.insert(rec.chrom.clone()) {
                    return Err(unsorted());
                }
                // Nothing merges across a chromosome boundary.
                for cluster in queue.iter_mut() {
                    cluster.open = false;
                }
                while let Some(cluster) = queue.pop_front() {
                    sink.emit(&cluster)?;
                }
            }
        }
        prev = Some((rec.chrom.clone(), rec.start));

        // The sortedness check above runs on every record, before this filter:
        // bedtools rejects an out-of-order record even when a strand flag would
        // have discarded it.
        if stranded && rec.strand().is_none() && !warned {
            warned = true;
            cli::warn(
                opts.quiet,
                &format!("{name}: no strand column; -s/-S can merge nothing"),
            );
        }
        let key = match (opts.strand, opts.same_strand) {
            (Some(want), _) => {
                if strand_of(&rec) != Some(want) {
                    continue;
                }
                Some(want)
            }
            (None, true) => match strand_of(&rec) {
                Some(s) => Some(s),
                None => continue,
            },
            (None, false) => None,
        };

        match queue.iter().position(|c| c.open && c.key == key) {
            Some(i) if queue[i].accepts(&rec, opts.d) => queue[i].add(&rec),
            Some(i) => {
                queue[i].open = false;
                queue.push_back(Cluster::new(key, &rec));
            }
            None => queue.push_back(Cluster::new(key, &rec)),
        }
        note_queue_len(queue.len());

        // ORACLE: output order is the order clusters were *opened*, not the
        // order they finished and not their merged start. Under `-s`, a `+`
        // cluster opened at 10 and running to 900 is printed before a `-`
        // cluster 30..40 that both opened and closed while it was still
        // growing. So a finished cluster waits behind any older one.
        while queue.front().is_some_and(|c| !c.open) {
            let cluster = queue.pop_front().expect("front checked");
            sink.emit(&cluster)?;
        }
    }

    if records == 0 {
        sink.set_headers(reader.headers().to_vec());
    }
    while let Some(cluster) = queue.pop_front() {
        sink.emit(&cluster)?;
    }
    sink.finish()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::sync::atomic::{AtomicU64, Ordering};

    static SEQ: AtomicU64 = AtomicU64::new(0);

    fn opts(input: &str) -> MergeOpts {
        MergeOpts {
            input: input.to_string(),
            d: 0,
            same_strand: false,
            strand: None,
            header: false,
            quiet: true,
        }
    }

    /// Run `merge` over `body` with `opts` already set up, returning stdout.
    fn merge(body: &str, tweak: impl FnOnce(&mut MergeOpts)) -> Result<String, bed::Error> {
        let mut path = std::env::temp_dir();
        path.push(format!(
            "mytools-merge-test-{}-{}",
            std::process::id(),
            SEQ.fetch_add(1, Ordering::Relaxed)
        ));
        File::create(&path)
            .unwrap()
            .write_all(body.as_bytes())
            .unwrap();

        let mut o = opts(path.to_str().unwrap());
        tweak(&mut o);
        let mut reader = bed::Reader::open(&o.input).unwrap();
        let mut sink: Vec<u8> = Vec::new();
        let result = merge_into(&mut reader, &o, &mut sink);
        std::fs::remove_file(&path).ok();
        result.map(|()| String::from_utf8(sink).unwrap())
    }

    /// The common case: defaults, expected to succeed.
    fn plain(body: &str) -> String {
        merge(body, |_| {}).unwrap()
    }

    // --- -d 0: overlapping and bookended ---

    #[test]
    fn bookended_features_merge_at_d_zero_without_overlapping() {
        // The two halves of the rule that trips everyone: `chr1 100 200` and
        // `chr1 200 300` share no base, yet `-d 0` merges them.
        assert!(!bed::overlaps(100, 200, 200, 300));
        assert_eq!(
            plain("chr1\t100\t200\nchr1\t200\t300\n"),
            "chr1\t100\t300\n"
        );
    }

    #[test]
    fn overlapping_features_merge_at_d_zero() {
        assert_eq!(
            plain("chr1\t100\t200\nchr1\t150\t250\n"),
            "chr1\t100\t250\n"
        );
    }

    #[test]
    fn features_with_a_gap_do_not_merge_at_d_zero() {
        assert_eq!(
            plain("chr1\t100\t200\nchr1\t201\t300\n"),
            "chr1\t100\t200\nchr1\t201\t300\n"
        );
    }

    #[test]
    fn nothing_merges_across_chromosomes() {
        assert_eq!(
            plain("chr1\t100\t200\nchr2\t100\t200\n"),
            "chr1\t100\t200\nchr2\t100\t200\n"
        );
    }

    // --- nested intervals and position 0 ---

    #[test]
    fn a_nested_interval_does_not_shrink_the_cluster() {
        // The cluster's end is the running maximum. If `200 300` were allowed
        // to set it, `450 600` would not reach and we would print two records.
        assert_eq!(
            plain("chr1\t100\t500\nchr1\t200\t300\nchr1\t450\t600\n"),
            "chr1\t100\t600\n"
        );
    }

    #[test]
    fn intervals_at_position_zero_merge_like_any_other() {
        assert_eq!(plain("chr1\t0\t10\nchr1\t10\t20\n"), "chr1\t0\t20\n");
        assert_eq!(
            plain("chr1\t0\t10\nchr1\t11\t20\n"),
            "chr1\t0\t10\nchr1\t11\t20\n"
        );
    }

    // --- -d N ---

    #[test]
    fn positive_d_merges_up_to_exactly_n_apart() {
        let d10 = |body: &str| merge(body, |o| o.d = 10).unwrap();
        assert_eq!(d10("chr1\t100\t200\nchr1\t210\t300\n"), "chr1\t100\t300\n");
        assert_eq!(
            d10("chr1\t100\t200\nchr1\t211\t300\n"),
            "chr1\t100\t200\nchr1\t211\t300\n"
        );
    }

    #[test]
    fn negative_d_requires_that_many_base_pairs_of_overlap() {
        let d5 = |body: &str| merge(body, |o| o.d = -5).unwrap();
        // 5 bp shared (195..199): merges. 4 bp shared: does not.
        assert_eq!(d5("chr1\t100\t200\nchr1\t195\t300\n"), "chr1\t100\t300\n");
        assert_eq!(
            d5("chr1\t100\t200\nchr1\t196\t300\n"),
            "chr1\t100\t200\nchr1\t196\t300\n"
        );
        // Bookended is zero overlap, so a negative -d splits what -d 0 joined.
        assert_eq!(
            d5("chr1\t100\t200\nchr1\t200\t300\n"),
            "chr1\t100\t200\nchr1\t200\t300\n"
        );
    }

    // --- zero-length features: all ORACLE behaviour, bedtools v2.31.1 ---

    #[test]
    fn a_lone_zero_length_feature_prints_its_original_coordinates() {
        // ORACLE: `bedtools merge` on a file holding only `chr1 500 500`
        // prints `chr1 500 500` — the inflation in `span` never shows.
        assert_eq!(plain("chr1\t500\t500\n"), "chr1\t500\t500\n");
        assert_eq!(
            plain("chr1\t500\t500\nchr1\t700\t800\n"),
            "chr1\t500\t500\nchr1\t700\t800\n"
        );
    }

    #[test]
    fn a_zero_length_feature_merges_one_base_leftwards() {
        // ORACLE: this is the `a.bed` case from the issue — `chr1 500 500` and
        // `chr1 500 600` merge to `chr1 499 600`, expanding LEFTWARDS past
        // anything either record covers.
        assert_eq!(
            plain("chr1\t500\t500\nchr1\t500\t600\n"),
            "chr1\t499\t600\n"
        );
    }

    #[test]
    fn a_zero_length_feature_merges_one_base_rightwards_too() {
        // ORACLE: the inflation is symmetric. `chr1 500 500` reaches 501, so a
        // feature starting at 501 is bookended with it and merges, one at 502
        // is not; and two copies of the same zero-length feature print as
        // `chr1 499 501`.
        assert_eq!(
            plain("chr1\t500\t500\nchr1\t501\t600\n"),
            "chr1\t499\t600\n"
        );
        assert_eq!(
            plain("chr1\t500\t500\nchr1\t502\t600\n"),
            "chr1\t500\t500\nchr1\t502\t600\n"
        );
        assert_eq!(
            plain("chr1\t500\t500\nchr1\t500\t500\n"),
            "chr1\t499\t501\n"
        );
        assert_eq!(
            plain("chr1\t100\t600\nchr1\t601\t601\n"),
            "chr1\t100\t602\n"
        );
    }

    #[test]
    fn a_zero_length_feature_at_zero_can_produce_a_negative_start() {
        // ORACLE: bedtools prints `chr1 -1 100` here. It does not clamp at 0,
        // so neither do we — hence signed coordinates throughout.
        assert_eq!(plain("chr1\t0\t0\nchr1\t0\t100\n"), "chr1\t-1\t100\n");
        // Alone, though, it is still printed verbatim.
        assert_eq!(plain("chr1\t0\t0\n"), "chr1\t0\t0\n");
    }

    #[test]
    fn zero_length_features_obey_negative_d_through_their_inflated_span() {
        // ORACLE: `chr1 500 500` against `chr1 500 600` shares exactly 1 bp
        // once inflated, so -d -1 merges and -d -2 does not.
        assert_eq!(
            merge("chr1\t500\t500\nchr1\t500\t600\n", |o| o.d = -1).unwrap(),
            "chr1\t499\t600\n"
        );
        assert_eq!(
            merge("chr1\t500\t500\nchr1\t500\t600\n", |o| o.d = -2).unwrap(),
            "chr1\t500\t500\nchr1\t500\t600\n"
        );
    }

    // --- strand ---

    #[test]
    fn same_strand_merges_only_within_a_strand() {
        let body = "chr1\t100\t200\ta\t0\t+\n\
                    chr1\t150\t250\tb\t0\t-\n\
                    chr1\t180\t400\tc\t0\t+\n\
                    chr1\t500\t600\td\t0\t-\n";
        assert_eq!(
            merge(body, |o| o.same_strand = true).unwrap(),
            "chr1\t100\t400\nchr1\t150\t250\nchr1\t500\t600\n"
        );
    }

    #[test]
    fn strand_clusters_are_printed_in_the_order_they_opened() {
        // ORACLE: not by merged start, and not by the order they finished.
        // The `+` run 10..900 opens second and finishes last, yet prints
        // second; the `-` run 30..40 opens and closes inside it and still
        // waits its turn.
        let body = "chr1\t10\t20\ta\t0\t-\n\
                    chr1\t10\t900\tb\t0\t+\n\
                    chr1\t30\t40\tc\t0\t-\n\
                    chr1\t800\t850\td\t0\t-\n";
        assert_eq!(
            merge(body, |o| o.same_strand = true).unwrap(),
            "chr1\t10\t20\nchr1\t10\t900\nchr1\t30\t40\nchr1\t800\t850\n"
        );
    }

    #[test]
    fn named_strand_keeps_only_that_strand() {
        let body = "chr1\t100\t200\ta\t0\t+\n\
                    chr1\t150\t250\tb\t0\t-\n\
                    chr1\t180\t400\tc\t0\t+\n";
        assert_eq!(
            merge(body, |o| o.strand = Some('+')).unwrap(),
            "chr1\t100\t400\n"
        );
        assert_eq!(
            merge(body, |o| o.strand = Some('-')).unwrap(),
            "chr1\t150\t250\n"
        );
    }

    #[test]
    fn records_whose_strand_is_neither_plus_nor_minus_are_dropped() {
        // ORACLE: `bedtools merge -s` over two overlapping `.` records prints
        // nothing — they do not form a group of their own. Without a strand
        // flag they merge like anything else.
        let body = "chr1\t100\t200\ta\t0\t.\nchr1\t150\t250\tb\t0\t.\n";
        assert_eq!(merge(body, |o| o.same_strand = true).unwrap(), "");
        assert_eq!(plain(body), "chr1\t100\t250\n");
    }

    #[test]
    fn strand_flags_on_strandless_input_produce_nothing_and_still_succeed() {
        // ORACLE: zero bytes on stdout, exit 0 — not an error (SPEC §3). The
        // warning goes to stderr, which stdout-diffing golden tests never see.
        let bed3 = "chr1\t100\t200\nchr1\t150\t250\n";
        assert_eq!(merge(bed3, |o| o.same_strand = true).unwrap(), "");
        assert_eq!(merge(bed3, |o| o.strand = Some('+')).unwrap(), "");
        // Without the flags the same input merges normally.
        assert_eq!(plain(bed3), "chr1\t100\t250\n");
    }

    // --- the unsorted-input detector ---

    #[test]
    fn a_start_going_backwards_is_a_data_error() {
        let e = merge("chr1\t100\t200\nchr1\t50\t60\n", |_| {}).unwrap_err();
        assert_eq!(e.msg, "input is not sorted");
        assert_eq!(e.line, Some(2));
        assert!(e.to_string().ends_with(":2: input is not sorted"));
    }

    #[test]
    fn the_unsorted_line_number_counts_header_lines() {
        let e = merge("#h\ntrack name=x\nchr1\t100\t200\nchr1\t50\t60\n", |_| {}).unwrap_err();
        assert_eq!(e.line, Some(4));
    }

    #[test]
    fn returning_to_a_finished_chromosome_is_a_data_error() {
        let e = merge("chr1\t100\t200\nchr2\t10\t20\nchr1\t300\t400\n", |_| {}).unwrap_err();
        assert_eq!(e.msg, "input is not sorted");
        assert_eq!(e.line, Some(3));
    }

    #[test]
    fn chromosomes_out_of_lexicographic_order_are_accepted() {
        // ORACLE: bedtools merge requires each chromosome's records to be
        // contiguous and ascending, not that the chromosomes themselves are
        // sorted. `chr2` before `chr1` runs clean and keeps that order.
        assert_eq!(
            plain("chr2\t100\t200\nchr1\t10\t20\n"),
            "chr2\t100\t200\nchr1\t10\t20\n"
        );
    }

    #[test]
    fn equal_starts_are_sorted_and_end_is_never_compared() {
        assert_eq!(
            plain("chr1\t100\t900\nchr1\t100\t200\n"),
            "chr1\t100\t900\n"
        );
    }

    #[test]
    fn an_out_of_order_record_is_rejected_even_when_a_strand_flag_would_drop_it() {
        // ORACLE: the check runs at read time, ahead of the strand filter.
        let body = "chr1\t100\t200\ta\t0\t+\n\
                    chr1\t50\t60\tb\t0\t-\n\
                    chr1\t300\t400\tc\t0\t+\n";
        let e = merge(body, |o| o.strand = Some('+')).unwrap_err();
        assert_eq!(e.msg, "input is not sorted");
        assert_eq!(e.line, Some(2));
    }

    // --- output shape ---

    #[test]
    fn output_is_bed3_whatever_the_input_width() {
        let out = plain("chr1\t100\t200\tname\t60\t+\textra\n");
        assert_eq!(out, "chr1\t100\t200\n");
        assert_eq!(out.trim_end().split('\t').count(), 3);
    }

    #[test]
    fn empty_input_prints_nothing() {
        assert_eq!(plain(""), "");
    }

    // --- -header ---

    #[test]
    fn header_flag_reproduces_the_leading_header_lines() {
        let body = "#comment\ntrack name=x\nchr1\t100\t200\nchr1\t150\t250\n";
        assert_eq!(
            merge(body, |o| o.header = true).unwrap(),
            "#comment\ntrack name=x\nchr1\t100\t250\n"
        );
        assert_eq!(plain(body), "chr1\t100\t250\n");
    }

    #[test]
    fn header_lines_after_the_first_record_are_not_reproduced() {
        // ORACLE: only the header block *above the first data record* is
        // replayed. bedtools prints `#top` here and swallows `#mid`...
        let body = "#top\nchr1\t100\t200\n#mid\nchr1\t150\t250\n";
        assert_eq!(
            merge(body, |o| o.header = true).unwrap(),
            "#top\nchr1\t100\t250\n"
        );
        // ...and where the data starts first, it prints no header at all,
        // however many header lines follow. `Reader::headers()` accumulates
        // every one it has seen, so this only holds because the list is
        // snapshotted the moment the first record arrives. Verified against
        // bedtools v2.31.1; the `header.bed` fixture has a clean leading block
        // and would pass either way, so the divergence would be silent.
        let body = "chr1\t100\t200\n#mid\ntrack x\nchr1\t150\t250\n";
        assert_eq!(
            merge(body, |o| o.header = true).unwrap(),
            "chr1\t100\t250\n"
        );
    }

    #[test]
    fn headers_are_printed_even_when_the_result_is_empty() {
        // ORACLE: `-header -s` on strandless input prints the headers and no
        // data. Empty *input* under `-header` prints the headers too.
        let body = "#comment\nchr1\t100\t200\nchr1\t150\t250\n";
        assert_eq!(
            merge(body, |o| {
                o.header = true;
                o.same_strand = true;
            })
            .unwrap(),
            "#comment\n"
        );
        assert_eq!(merge("#only\n", |o| o.header = true).unwrap(), "#only\n");
    }

    // --- streaming (SPEC §6) ---

    #[test]
    fn only_the_open_clusters_are_ever_held() {
        // The memory claim of SPEC §6, asserted rather than asserted-in-prose:
        // however many records go past, the queue holds at most one open
        // cluster per strand group — plus, for the instant between opening a
        // cluster and flushing the queue, the one it displaced. So two without
        // `-s` and three with it, whatever the input size. A buffering
        // implementation would end up holding all 4000 of these.
        let mut body = String::new();
        for i in 0..2000u64 {
            let s = i * 10;
            body.push_str(&format!("chr1\t{}\t{}\ta\t0\t+\n", s, s + 5));
            body.push_str(&format!("chr1\t{}\t{}\tb\t0\t-\n", s + 1, s + 6));
        }

        PEAK_QUEUE.with(|p| p.set(0));
        let out = plain(&body);
        assert_eq!(out.lines().count(), 2000, "every pair merges without -s");
        assert_eq!(PEAK_QUEUE.with(|p| p.get()), 2);

        PEAK_QUEUE.with(|p| p.set(0));
        let out = merge(&body, |o| o.same_strand = true).unwrap();
        assert_eq!(
            out.lines().count(),
            4000,
            "no two share both strand and gap"
        );
        assert_eq!(PEAK_QUEUE.with(|p| p.get()), 3);
    }
}
