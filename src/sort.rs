//! `mytools sort` — `SPEC.md` §4, §5, §6.
//!
//! # Ordering
//!
//! Chromosome **lexicographically**, then `start` ascending. Nothing else.
//! The chromosome comparison is a raw byte comparison, so `chr10` sorts before
//! `chr2` and `chr17` before `chr7` (and `Chr2` before `chr2`). Natural/numeric
//! ordering was considered and rejected: it breaks oracle parity on any genome
//! with `chr10`-style names. Verified against bedtools v2.31.1.
//!
//! `end` is **never** compared and the sort is **stable**: records tying on
//! `(chrom, start)` come out in input order. Verified against bedtools v2.31.1
//! with `chr1 100 900`, `chr1 100 200`, `chr1 100 500` in that input order —
//! bedtools emits them in that order, not by ascending end.
//!
//! # Memory: spill and collect
//!
//! `sort` is the one command allowed not to stream, and it may hold **one
//! chromosome** at a time — never the whole file (`CLAUDE.md`, `SPEC.md` §6).
//! A chromosome's records are not contiguous in the input, and lexicographic
//! output order means nothing can be emitted before EOF, so a running buffer
//! cannot work. Hence two phases:
//!
//! 1. **Spill.** Stream the input once, appending each record to a temp file
//!    for its chromosome. Only the chromosome names, one line, and a bounded
//!    set of open writers are in memory.
//! 2. **Collect.** Walk the chromosome names in order; for each, read *that
//!    chromosome's* spill file into a vector, stable-sort by start, write it
//!    out, and drop it before touching the next.
//!
//! The high-water mark is therefore one chromosome, whatever the file size.

use std::cmp::Ordering;
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process;
use std::time::{SystemTime, UNIX_EPOCH};

use crate::bed;
use crate::cli::SortOpts;

/// How many spill files may be open at once. Beyond this they are closed and
/// reopened in append mode on demand, so a reference with thousands of contigs
/// cannot exhaust the process' file descriptors.
const MAX_OPEN_SPILLS: usize = 64;

/// Chromosome order: a raw byte comparison, which is what bedtools does.
/// `chr10` < `chr2`, `chr17` < `chr7`, `Chr2` < `chr2`. Oracle behaviour,
/// verified against bedtools v2.31.1.
fn chrom_order(a: &str, b: &str) -> Ordering {
    a.as_bytes().cmp(b.as_bytes())
}

/// The full record ordering: chromosome, then `start`. **`end` is not a
/// tiebreaker** — records tying here compare `Equal` and a stable sort then
/// leaves them in input order. Oracle behaviour, verified against bedtools
/// v2.31.1; no fixture exposes it, so it is pinned by unit test only.
fn record_order(a: &Spilled, b: &Spilled) -> Ordering {
    chrom_order(&a.chrom, &b.chrom).then(a.start.cmp(&b.start))
}

/// One record on its way through the spill: the sort key plus the original
/// line, replayed verbatim so trailing columns survive untouched (`SPEC.md` §5).
#[derive(Debug, Clone, PartialEq, Eq)]
struct Spilled {
    chrom: String,
    start: u64,
    line: String,
}

pub fn run(opts: &SortOpts) -> Result<(), bed::Error> {
    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    sort_to(opts, &mut out)?;
    out.flush().map_err(|_| write_failed("stdout"))
}

/// The whole of `sort`, writing to any sink so the unit tests can read it back.
fn sort_to(opts: &SortOpts, out: &mut dyn Write) -> Result<(), bed::Error> {
    let mut reader = bed::Reader::open(&opts.input)?;
    let mut spill = Spill::new()?;

    // Headers are captured as they stood when the first data record arrived.
    // bedtools only replays the header block at the top of the file: given
    // `chr2 5 6`, `#mid`, `chr1 1 2`, `-header` emits no header at all.
    // Oracle behaviour, verified against bedtools v2.31.1. The shared reader
    // collects every `#`/`track`/`browser` line it sees, mid-file ones
    // included, so sort must take the snapshot itself rather than ask at EOF.
    let mut headers: Option<Vec<String>> = None;

    while let Some(rec) = reader.next() {
        let rec = rec?;
        if headers.is_none() {
            headers = Some(reader.headers().to_vec());
        }
        spill.push(&rec)?;
    }
    // No data records at all: every header line in the file is still a header,
    // and bedtools replays them (a header-only file under `-header` prints the
    // headers and exits 0). Verified against bedtools v2.31.1.
    let headers = headers.unwrap_or_else(|| reader.headers().to_vec());

    if opts.header {
        for h in &headers {
            writeln!(out, "{h}").map_err(|_| write_failed("stdout"))?;
        }
    }

    // Collect: one chromosome in memory at a time, in lexicographic order.
    for chrom in spill.chroms()? {
        let mut batch = spill.load(&chrom)?;
        // Stable, so a `(chrom, start)` tie keeps input order.
        batch.sort_by(record_order);
        for rec in &batch {
            // Always terminated, even when the input's last line was not.
            writeln!(out, "{}", rec.line).map_err(|_| write_failed("stdout"))?;
        }
        // Explicit: the batch dies here, before the next chromosome is read.
        drop(batch);
    }

    Ok(())
}

fn write_failed(what: &str) -> bed::Error {
    bed::Error::file(what, "cannot be written")
}

fn read_failed(path: &Path) -> bed::Error {
    bed::Error::file(&path.display().to_string(), "cannot be read")
}

/// The spill area: one temp file per chromosome, plus the mapping from
/// chromosome name to file. Files are named by index, not by chromosome, so a
/// contig called `../etc` or `a/b` cannot escape the directory.
///
/// Memory here is O(number of chromosomes), not O(number of records).
struct Spill {
    dir: PathBuf,
    index: HashMap<String, usize>,
    writers: HashMap<usize, BufWriter<File>>,
}

impl Spill {
    fn new() -> Result<Spill, bed::Error> {
        let base = std::env::temp_dir();
        let stamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        for attempt in 0..64 {
            let dir = base.join(format!("mytools-sort-{}-{stamp}-{attempt}", process::id()));
            match fs::create_dir(&dir) {
                Ok(()) => {
                    return Ok(Spill {
                        dir,
                        index: HashMap::new(),
                        writers: HashMap::new(),
                    });
                }
                Err(e) if e.kind() == io::ErrorKind::AlreadyExists => continue,
                Err(_) => break,
            }
        }
        Err(bed::Error::file(
            &base.display().to_string(),
            "cannot create a temporary directory for sort",
        ))
    }

    fn path(&self, idx: usize) -> PathBuf {
        self.dir.join(format!("{idx}.spill"))
    }

    /// Append one record to its chromosome's file. The spill line is
    /// `<start>\t<original line>`: `start` is re-read on the way back rather
    /// than re-parsed out of the record, and BED lines never contain a newline
    /// (the reader strips terminators), so one spilled record is one line.
    fn push(&mut self, rec: &bed::Record) -> Result<(), bed::Error> {
        let next = self.index.len();
        let idx = *self.index.entry(rec.chrom.clone()).or_insert(next);

        if !self.writers.contains_key(&idx) && self.writers.len() >= MAX_OPEN_SPILLS {
            self.close_all()?;
        }
        let path = self.path(idx);
        let writer = match self.writers.entry(idx) {
            Entry::Occupied(e) => e.into_mut(),
            Entry::Vacant(e) => {
                // create+append, so reopening after a close just continues.
                let f = OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(&path)
                    .map_err(|_| write_failed(&path.display().to_string()))?;
                e.insert(BufWriter::new(f))
            }
        };
        writeln!(writer, "{}\t{}", rec.start, rec.line)
            .map_err(|_| write_failed(&path.display().to_string()))
    }

    /// Flush and close every open spill file.
    fn close_all(&mut self) -> Result<(), bed::Error> {
        for (idx, mut w) in self.writers.drain() {
            w.flush().map_err(|_| {
                write_failed(&self.dir.join(format!("{idx}.spill")).display().to_string())
            })?;
        }
        Ok(())
    }

    /// Every chromosome seen, in output order. Names only — no records.
    fn chroms(&mut self) -> Result<Vec<String>, bed::Error> {
        self.close_all()?;
        let mut names: Vec<String> = self.index.keys().cloned().collect();
        names.sort_by(|a, b| chrom_order(a, b));
        Ok(names)
    }

    /// Read back exactly one chromosome, in input order. This is the only
    /// place records live in memory, and it holds one chromosome's worth.
    fn load(&mut self, chrom: &str) -> Result<Vec<Spilled>, bed::Error> {
        self.close_all()?;
        let Some(&idx) = self.index.get(chrom) else {
            return Ok(Vec::new());
        };
        let path = self.path(idx);
        let f = File::open(&path).map_err(|_| read_failed(&path))?;
        let mut out = Vec::new();
        for line in BufReader::new(f).lines() {
            let line = line.map_err(|_| read_failed(&path))?;
            // Our own format, written by push: `<start>\t<line>`.
            let (start, rest) = line.split_once('\t').ok_or_else(|| read_failed(&path))?;
            let start: u64 = start.parse().map_err(|_| read_failed(&path))?;
            out.push(Spilled {
                chrom: chrom.to_string(),
                start,
                line: rest.to_string(),
            });
        }
        Ok(out)
    }
}

impl Drop for Spill {
    fn drop(&mut self) {
        self.writers.clear();
        // Best effort: a leftover temp directory is not worth an error path.
        let _ = fs::remove_dir_all(&self.dir);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};

    static COUNTER: AtomicU64 = AtomicU64::new(0);

    fn tmpfile(body: &[u8]) -> PathBuf {
        let n = COUNTER.fetch_add(1, AtomicOrdering::SeqCst);
        let mut p = std::env::temp_dir();
        p.push(format!("mytools-sort-test-{}-{n}.bed", process::id()));
        let mut f = File::create(&p).unwrap();
        f.write_all(body).unwrap();
        p
    }

    fn opts(input: &str, header: bool) -> SortOpts {
        SortOpts {
            input: input.to_string(),
            header,
            quiet: false,
        }
    }

    /// Sort `body` and return stdout as a string.
    fn sorted(body: &str) -> String {
        sorted_opt(body, false).unwrap()
    }

    fn sorted_opt(body: &str, header: bool) -> Result<String, bed::Error> {
        let p = tmpfile(body.as_bytes());
        let mut out: Vec<u8> = Vec::new();
        let r = sort_to(&opts(p.to_str().unwrap(), header), &mut out);
        fs::remove_file(&p).ok();
        r.map(|()| String::from_utf8(out).unwrap())
    }

    fn spilled(chrom: &str, start: u64) -> Spilled {
        Spilled {
            chrom: chrom.to_string(),
            start,
            line: format!("{chrom}\t{start}\t{}", start + 1),
        }
    }

    // --- the comparator itself ---

    #[test]
    fn chromosomes_compare_lexicographically_not_numerically() {
        // The case a.bed cannot show: chr10 before chr2, chr17 before chr7.
        // Oracle behaviour, verified against bedtools v2.31.1.
        assert_eq!(chrom_order("chr10", "chr2"), Ordering::Less);
        assert_eq!(chrom_order("chr17", "chr7"), Ordering::Less);
        assert_eq!(chrom_order("chr2", "chrX"), Ordering::Less);
        assert_eq!(chrom_order("chr2", "chr2"), Ordering::Equal);
        // Raw bytes, so uppercase sorts first: bedtools puts Chr2 before chr2.
        assert_eq!(chrom_order("Chr2", "chr2"), Ordering::Less);
    }

    #[test]
    fn record_order_is_chrom_then_start_only() {
        assert_eq!(
            record_order(&spilled("chr1", 500), &spilled("chr2", 0)),
            Ordering::Less
        );
        assert_eq!(
            record_order(&spilled("chr1", 100), &spilled("chr1", 200)),
            Ordering::Less
        );
    }

    #[test]
    fn record_order_never_compares_end() {
        // Same chrom and start, wildly different ends: Equal, so a stable sort
        // leaves them in input order. Oracle behaviour (bedtools v2.31.1).
        let long = Spilled {
            chrom: "chr1".into(),
            start: 100,
            line: "chr1\t100\t900".into(),
        };
        let short = Spilled {
            chrom: "chr1".into(),
            start: 100,
            line: "chr1\t100\t200".into(),
        };
        assert_eq!(record_order(&long, &short), Ordering::Equal);
        assert_eq!(record_order(&short, &long), Ordering::Equal);
    }

    // --- ordering, end to end ---

    #[test]
    fn chr10_before_chr2_end_to_end() {
        let got = sorted("chr2\t50\t60\nchr10\t5\t6\nchr2\t10\t20\nchr1\t7\t8\nchr10\t1\t2\n");
        assert_eq!(
            got, "chr1\t7\t8\nchr10\t1\t2\nchr10\t5\t6\nchr2\t10\t20\nchr2\t50\t60\n",
            "chromosomes must sort lexicographically"
        );
    }

    #[test]
    fn chromosome_records_need_not_be_contiguous_in_the_input() {
        // The spill-and-collect case: chr1's records are split by a chr2 line.
        let got = sorted("chr1\t9\t10\nchr2\t0\t1\nchr1\t1\t2\n");
        assert_eq!(got, "chr1\t1\t2\nchr1\t9\t10\nchr2\t0\t1\n");
    }

    #[test]
    fn start_tie_keeps_input_order_and_ignores_end() {
        // The case no golden test can catch: a.bed, b.bed and genes.bed all
        // sort identically under "then end ascending". bedtools v2.31.1 emits
        // these three in input order, 900 first. Oracle behaviour.
        let got = sorted("chr1\t100\t900\tx1\nchr1\t100\t200\tx2\nchr1\t100\t500\tx3\n");
        assert_eq!(
            got,
            "chr1\t100\t900\tx1\nchr1\t100\t200\tx2\nchr1\t100\t500\tx3\n"
        );
    }

    #[test]
    fn full_ties_are_stable() {
        // a03/a04 and a09/a10 from data/a.bed: identical coordinates, distinct
        // names. bedtools keeps them in input order.
        let got = sorted(
            "chr1\t700\t800\ta09\t60\t+\n\
             chr1\t150\t250\ta03\t30\t+\n\
             chr1\t150\t250\ta04\t30\t-\n\
             chr1\t700\t800\ta10\t60\t+\n",
        );
        assert_eq!(
            got,
            "chr1\t150\t250\ta03\t30\t+\n\
             chr1\t150\t250\ta04\t30\t-\n\
             chr1\t700\t800\ta09\t60\t+\n\
             chr1\t700\t800\ta10\t60\t+\n"
        );
    }

    #[test]
    fn zero_length_and_position_zero_records_sort_by_start_like_the_oracle() {
        // `chr2 0 0` leads chr2 and `chr1 0 100` leads chr1; `chr2 300 300`
        // precedes `chr2 300 400` only because it came first in the input, not
        // because its end is smaller. Verified against bedtools v2.31.1 on
        // data/a.bed, which contains exactly these records (a12, a01, a16, a17).
        let got = sorted(
            "chr2\t300\t300\ta16\t0\t+\n\
             chr2\t300\t400\ta17\t65\t+\n\
             chr2\t0\t0\ta12\t0\t+\n\
             chr1\t100\t200\ta02\t20\t-\n\
             chr1\t0\t100\ta01\t10\t+\n",
        );
        assert_eq!(
            got,
            "chr1\t0\t100\ta01\t10\t+\n\
             chr1\t100\t200\ta02\t20\t-\n\
             chr2\t0\t0\ta12\t0\t+\n\
             chr2\t300\t300\ta16\t0\t+\n\
             chr2\t300\t400\ta17\t65\t+\n"
        );
    }

    // --- output shape ---

    #[test]
    fn bed3_and_bed6_are_echoed_at_their_original_width() {
        // Trailing columns are untouched, and a wider record is not padded or
        // trimmed to match a narrower one.
        let got = sorted(
            "chr1\t200\t300\tname\t5\t+\textra\tcolumns\n\
             chr1\t100\t200\n\
             chr1\t150\t250\tn4\n",
        );
        assert_eq!(
            got,
            "chr1\t100\t200\n\
             chr1\t150\t250\tn4\n\
             chr1\t200\t300\tname\t5\t+\textra\tcolumns\n"
        );
    }

    #[test]
    fn the_last_line_is_terminated_even_when_the_input_is_not() {
        assert_eq!(sorted("chr1\t5\t6\tz"), "chr1\t5\t6\tz\n");
    }

    #[test]
    fn empty_input_produces_zero_bytes() {
        assert_eq!(sorted(""), "");
        // And with -header, since there are no headers to replay either.
        assert_eq!(sorted_opt("", true).unwrap(), "");
    }

    // --- headers ---

    #[test]
    fn headers_are_replayed_only_with_the_header_flag() {
        let body = "#comment line\n\
                    track name=fixture\n\
                    browser position chr1:1-1000\n\
                    chr2\t0\t1\n\
                    chr1\t0\t1\n";
        assert_eq!(
            sorted_opt(body, true).unwrap(),
            "#comment line\n\
             track name=fixture\n\
             browser position chr1:1-1000\n\
             chr1\t0\t1\n\
             chr2\t0\t1\n"
        );
        assert_eq!(sorted_opt(body, false).unwrap(), "chr1\t0\t1\nchr2\t0\t1\n");
    }

    #[test]
    fn header_lines_after_the_first_record_are_dropped() {
        // Oracle behaviour: bedtools v2.31.1 on `chr2 5 6`, `#mid`, `chr1 1 2`,
        // `track x`, `chr1 0 1` with -header prints the three records and no
        // header line at all. Do not "fix" this to replay every header.
        let body = "chr2\t5\t6\n#mid\nchr1\t1\t2\ntrack x\nchr1\t0\t1\n";
        assert_eq!(
            sorted_opt(body, true).unwrap(),
            "chr1\t0\t1\nchr1\t1\t2\nchr2\t5\t6\n"
        );
    }

    #[test]
    fn a_header_only_file_still_replays_its_headers() {
        // No data records, so every header line is part of the leading block.
        // Verified against bedtools v2.31.1: exit 0, headers on stdout.
        let body = "#comment line\ntrack name=fixture\n";
        assert_eq!(sorted_opt(body, true).unwrap(), body);
        assert_eq!(sorted_opt(body, false).unwrap(), "");
    }

    // --- data errors (SPEC §7): message exact, exit 1 via main ---

    #[test]
    fn fewer_than_three_columns_is_a_data_error() {
        let e = sorted_opt("chr1\t100\n", false).unwrap_err();
        assert_eq!(e.msg, "fewer than 3 columns");
        assert_eq!(e.line, Some(1));
    }

    #[test]
    fn non_integer_coordinates_are_a_data_error() {
        let e = sorted_opt("chr1\t0\t100\nchr1\tten\t100\n", false).unwrap_err();
        assert_eq!(e.msg, "start/end must be integers");
        assert_eq!(e.line, Some(2));
    }

    #[test]
    fn start_greater_than_end_is_a_data_error() {
        let e = sorted_opt("chr1\t300\t200\n", false).unwrap_err();
        assert_eq!(e.msg, "start greater than end");
    }

    #[test]
    fn a_missing_input_file_is_a_data_error() {
        let e = sort_to(&opts("/nonexistent/nope.bed", false), &mut Vec::new()).unwrap_err();
        assert_eq!(e.msg, "no such file");
        assert_eq!(
            e.to_string(),
            "mytools: /nonexistent/nope.bed: no such file"
        );
    }

    #[test]
    fn a_data_error_writes_nothing_to_stdout() {
        // Nothing is emitted before EOF, so a bad line at the end still leaves
        // stdout untouched — errors never half-write a result.
        let p = tmpfile(b"chr1\t0\t100\nchr1\tbad\t1\n");
        let mut out: Vec<u8> = Vec::new();
        let r = sort_to(&opts(p.to_str().unwrap(), false), &mut out);
        fs::remove_file(&p).ok();
        assert!(r.is_err());
        assert!(out.is_empty(), "stdout was written before the error");
    }

    // --- the memory contract ---

    #[test]
    fn the_spill_hands_back_one_chromosome_at_a_time() {
        // The acceptance criterion made testable: load() returns the records of
        // exactly one chromosome, in input order, never the whole file.
        let mut spill = Spill::new().unwrap();
        for line in [
            "chr2\t50\t60",
            "chr10\t5\t6",
            "chr2\t10\t20",
            "chr1\t7\t8",
            "chr10\t1\t2",
        ] {
            let cols: Vec<&str> = line.split('\t').collect();
            spill
                .push(&bed::Record {
                    chrom: cols[0].to_string(),
                    start: cols[1].parse().unwrap(),
                    end: cols[2].parse().unwrap(),
                    line: line.to_string(),
                })
                .unwrap();
        }
        assert_eq!(spill.chroms().unwrap(), ["chr1", "chr10", "chr2"]);

        let chr2 = spill.load("chr2").unwrap();
        assert_eq!(chr2.len(), 2, "load() returned more than one chromosome");
        // Input order preserved on the way back, which is what makes the
        // stable sort produce input order for ties.
        assert_eq!(chr2[0].start, 50);
        assert_eq!(chr2[1].start, 10);
        assert_eq!(spill.load("chr1").unwrap().len(), 1);
        assert_eq!(spill.load("chrZ").unwrap().len(), 0);
    }

    #[test]
    fn the_spill_survives_more_chromosomes_than_it_may_keep_open() {
        // Past MAX_OPEN_SPILLS the writers are closed and reopened in append
        // mode; nothing may be lost or truncated when that happens.
        let n = MAX_OPEN_SPILLS * 2 + 3;
        let mut body = String::new();
        // Two passes, so every chromosome is written to after a reopen.
        for _ in 0..2 {
            for i in 0..n {
                body.push_str(&format!("c{i:04}\t{i}\t{}\n", i + 1));
            }
        }
        let got = sorted(&body);
        assert_eq!(got.lines().count(), 2 * n);
        // c0000 twice, then c0001 twice, ... — lexicographic, then stable.
        let first: Vec<&str> = got.lines().take(4).collect();
        assert_eq!(
            first,
            ["c0000\t0\t1", "c0000\t0\t1", "c0001\t1\t2", "c0001\t1\t2"]
        );
    }

    #[test]
    fn gzipped_input_sorts_identically_to_plain_input() {
        use flate2::Compression;
        use flate2::write::GzEncoder;

        let body = "chr17\t7\t8\tb\nchr7\t1\t2\ta\nchr17\t1\t2\tc\n";
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(body.as_bytes()).unwrap();
        let p = tmpfile(&enc.finish().unwrap());
        let mut out: Vec<u8> = Vec::new();
        sort_to(&opts(p.to_str().unwrap(), false), &mut out).unwrap();
        fs::remove_file(&p).ok();
        assert_eq!(String::from_utf8(out).unwrap(), sorted(body));
        // And chr17 leads chr7, as the oracle has it.
        assert_eq!(
            sorted(body),
            "chr17\t1\t2\tc\nchr17\t7\t8\tb\nchr7\t1\t2\ta\n"
        );
    }
}
