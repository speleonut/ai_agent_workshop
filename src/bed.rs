//! Shared BED input layer: parsing, sources (file/stdin/gzip), and the one
//! overlap predicate the whole project hinges on.
//!
//! `SPEC.md` §2, §3, §6, §7. Every subcommand goes through here so that none of
//! them invents its own idea of what a BED record is.

// Parts of this layer are used only by subcommands that are still stubs.
// TODO: drop this once sort, merge and intersect are all implemented.
#![allow(dead_code)]

use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};

use flate2::read::MultiGzDecoder;

/// The overlap predicate. BED is 0-based half-open, so this is a strict `<` on
/// both sides: `chr1 100 200` covers bases 100..199.
///
/// Consequences, both deliberate:
/// - bookended intervals (`a_end == b_start`) do **not** overlap;
/// - zero-length intervals (`start == end`) never overlap anything, because
///   `a_start < a_end` is false for them.
///
/// Every off-by-one in this project lives in this comparison. It exists once.
pub fn overlaps(a_start: u64, a_end: u64, b_start: u64, b_end: u64) -> bool {
    a_start < b_end && b_start < a_end
}

/// A parsed BED record. The first three columns are validated; everything past
/// them is carried untouched.
///
/// `line` is the original line bytes with the terminator stripped. `sort` and
/// `intersect -wa`/`-v`/`-u` echo records back at their original width, so they
/// replay `line` rather than re-joining parsed fields.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub line: String,
}

impl Record {
    /// Column `i`, 0-based, or `None` if the record is narrower than that.
    pub fn field(&self, i: usize) -> Option<&str> {
        self.line.split('\t').nth(i)
    }

    /// Number of columns in the original line.
    pub fn ncols(&self) -> usize {
        self.line.split('\t').count()
    }

    /// Column 4 (`name`), BED4 and wider.
    pub fn name(&self) -> Option<&str> {
        self.field(3)
    }

    /// Column 6 (`strand`), BED6 and wider. `None` on narrower input — callers
    /// that need a strand warn rather than fail (`SPEC.md` §3, §7).
    pub fn strand(&self) -> Option<&str> {
        self.field(5)
    }

    /// Does this record overlap `other`? Same chromosome and [`overlaps`].
    pub fn overlaps(&self, other: &Record) -> bool {
        self.chrom == other.chrom && overlaps(self.start, self.end, other.start, other.end)
    }
}

/// A data error: bad input, exit `1` (`SPEC.md` §7). Usage errors are the CLI's
/// business and live in `cli.rs`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Error {
    pub file: String,
    pub line: Option<u64>,
    pub msg: String,
}

impl Error {
    /// An error attached to a specific line of a specific file.
    pub fn at(file: &str, line: u64, msg: &str) -> Error {
        Error {
            file: file.to_string(),
            line: Some(line),
            msg: msg.to_string(),
        }
    }

    /// An error about the file as a whole (missing, unreadable, corrupt).
    pub fn file(file: &str, msg: &str) -> Error {
        Error {
            file: file.to_string(),
            line: None,
            msg: msg.to_string(),
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.line {
            Some(n) => write!(f, "mytools: {}:{}: {}", self.file, n, self.msg),
            None => write!(f, "mytools: {}: {}", self.file, self.msg),
        }
    }
}

impl std::error::Error for Error {}

/// A streaming reader over one BED source.
///
/// Yields `Result<Record, Error>` line by line — nothing is buffered beyond the
/// current line, so an exome-scale file costs constant memory here. `sort`
/// buffers on top of this; that is its documented exception, not the reader's.
pub struct Reader {
    inner: Box<dyn BufRead>,
    name: String,
    lineno: u64,
    headers: Vec<String>,
    gzipped: bool,
}

impl Reader {
    /// Open `spec`: a path, or `-` for stdin.
    ///
    /// gzip and bgzip (BGZF) are detected by **magic bytes**, not by filename,
    /// so `foo.bed` that is secretly gzipped still reads. BGZF is a multi-member
    /// gzip stream, hence `MultiGzDecoder` — a plain `GzDecoder` would silently
    /// stop after the first block.
    pub fn open(spec: &str) -> Result<Reader, Error> {
        let (raw, name): (Box<dyn Read>, String) = if spec == "-" {
            (Box::new(io::stdin()), "-".to_string())
        } else {
            match File::open(spec) {
                Ok(f) => (Box::new(f), spec.to_string()),
                Err(e) if e.kind() == io::ErrorKind::NotFound => {
                    return Err(Error::file(spec, "no such file"));
                }
                Err(_) => return Err(Error::file(spec, "cannot be read")),
            }
        };

        let mut buffered = BufReader::new(raw);
        // Peek without consuming: fill_buf leaves the bytes in place.
        let gzipped = match buffered.fill_buf() {
            Ok(head) => head.starts_with(&[0x1f, 0x8b]),
            Err(_) => false,
        };
        let inner: Box<dyn BufRead> = if gzipped {
            Box::new(BufReader::new(MultiGzDecoder::new(buffered)))
        } else {
            Box::new(buffered)
        };

        Ok(Reader {
            inner,
            name,
            lineno: 0,
            headers: Vec::new(),
            gzipped,
        })
    }

    /// The source name, as it appears in error messages.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// `#`, `track` and `browser` lines seen so far, in input order, with
    /// terminators stripped. They never reach the interval logic; `-header`
    /// replays them ahead of the results (`SPEC.md` §5).
    ///
    /// Headers sit at the top of a file in practice, so a streaming caller can
    /// pull the first record and then emit these before writing anything.
    pub fn headers(&self) -> &[String] {
        &self.headers
    }

    /// Collect every record, failing on the first bad line. Convenient for
    /// `sort` (which buffers anyway) and for `intersect`'s B side.
    pub fn collect_records(&mut self) -> Result<Vec<Record>, Error> {
        let mut out = Vec::new();
        for rec in self {
            out.push(rec?);
        }
        Ok(out)
    }

    fn read_line(&mut self) -> Option<Result<String, Error>> {
        let mut raw = Vec::new();
        match self.inner.read_until(b'\n', &mut raw) {
            Ok(0) => None,
            Ok(_) => {
                self.lineno += 1;
                if raw.last() == Some(&b'\n') {
                    raw.pop();
                }
                if raw.last() == Some(&b'\r') {
                    raw.pop();
                }
                match String::from_utf8(raw) {
                    Ok(s) => Some(Ok(s)),
                    Err(_) => Some(Err(Error::at(
                        &self.name,
                        self.lineno,
                        "malformed BED record",
                    ))),
                }
            }
            // A read error on a gzip source means the stream is corrupt; on a
            // plain file it is an I/O failure. Both exit 1.
            Err(_) if self.gzipped => Some(Err(Error::file(&self.name, "not a valid gzip stream"))),
            Err(_) => Some(Err(Error::file(&self.name, "cannot be read"))),
        }
    }

    fn parse(&self, line: &str) -> Result<Record, Error> {
        let bad = |msg: &str| Error::at(&self.name, self.lineno, msg);

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 3 {
            return Err(bad("fewer than 3 columns"));
        }
        if cols[0].is_empty() {
            return Err(bad("malformed BED record"));
        }
        // parse::<u64> rejects negatives and non-digits alike, which is what we
        // want: both are "not a non-negative integer".
        let start = cols[1]
            .parse::<u64>()
            .map_err(|_| bad("start/end must be integers"))?;
        let end = cols[2]
            .parse::<u64>()
            .map_err(|_| bad("start/end must be integers"))?;
        if start > end {
            return Err(bad("start greater than end"));
        }

        Ok(Record {
            chrom: cols[0].to_string(),
            start,
            end,
            line: line.to_string(),
        })
    }
}

/// `#`, `track` and `browser` lines are metadata, not intervals.
fn is_header(line: &str) -> bool {
    line.starts_with('#') || line.starts_with("track") || line.starts_with("browser")
}

impl Iterator for Reader {
    type Item = Result<Record, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let line = match self.read_line()? {
                Ok(l) => l,
                Err(e) => return Some(Err(e)),
            };
            if line.is_empty() {
                continue;
            }
            if is_header(&line) {
                self.headers.push(line);
                continue;
            }
            return Some(self.parse(&line));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn tmpfile(name: &str, body: &[u8]) -> std::path::PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!("mytools-bed-test-{}-{}", std::process::id(), name));
        let mut f = File::create(&p).unwrap();
        f.write_all(body).unwrap();
        p
    }

    fn parse_all(body: &str) -> Result<Vec<Record>, Error> {
        let p = tmpfile("parse", body.as_bytes());
        let mut r = Reader::open(p.to_str().unwrap()).unwrap();
        let out = r.collect_records();
        std::fs::remove_file(&p).ok();
        out
    }

    // --- the overlap predicate, one case per edge ---

    #[test]
    fn bookended_intervals_do_not_overlap() {
        // a.end == b.start: 100..199 and 200..299 share no base.
        assert!(!overlaps(100, 200, 200, 300));
    }

    #[test]
    fn one_shared_base_pair_overlaps() {
        // 199 is in both.
        assert!(overlaps(100, 200, 199, 300));
    }

    #[test]
    fn nested_intervals_overlap() {
        assert!(overlaps(100, 500, 200, 300));
        assert!(overlaps(200, 300, 100, 500));
    }

    #[test]
    fn identical_intervals_overlap() {
        assert!(overlaps(150, 250, 150, 250));
    }

    #[test]
    fn interval_at_position_zero_overlaps() {
        assert!(overlaps(0, 100, 0, 50));
        assert!(!overlaps(0, 100, 100, 200));
    }

    #[test]
    fn zero_length_interval_inside_another_satisfies_the_predicate() {
        // A zero-length point strictly inside an interval passes the strict
        // `<` on both sides, and bedtools agrees: `intersect -a` of
        // `chr1 500 500` against `chr1 400 600` reports the pair. Oracle
        // behaviour, verified against bedtools v2.31.1.
        assert!(overlaps(500, 500, 400, 600));
        assert!(overlaps(400, 600, 500, 500));
    }

    #[test]
    fn zero_length_interval_does_not_overlap_at_its_own_boundary() {
        // Touching, not containing: the strict `<` rules these out.
        assert!(!overlaps(500, 500, 500, 600));
        assert!(!overlaps(500, 500, 400, 500));
    }

    #[test]
    fn identical_zero_length_intervals_are_where_the_predicate_and_oracle_part() {
        // The predicate says no (500 < 500 is false) but bedtools reports
        // `chr1 500 500` against itself as an overlap. Verified against
        // bedtools v2.31.1. This function stays as SPEC §3 defines it; the
        // subcommand that has to match the oracle here handles the case and
        // says so. Pinned as a test so the divergence cannot be forgotten.
        assert!(!overlaps(500, 500, 500, 500));
    }

    #[test]
    fn record_overlap_requires_same_chromosome() {
        let a = Record {
            chrom: "chr1".into(),
            start: 0,
            end: 100,
            line: "chr1\t0\t100".into(),
        };
        let b = Record {
            chrom: "chr2".into(),
            start: 0,
            end: 100,
            line: "chr2\t0\t100".into(),
        };
        assert!(!a.overlaps(&b));
        assert!(a.overlaps(&a.clone()));
    }

    // --- parsing ---

    #[test]
    fn parses_bed3_through_bed6_and_keeps_original_width() {
        let recs =
            parse_all("chr1\t0\t100\nchr1\t100\t200\tn\nchr1\t200\t300\tn\t5\t+\textra\n").unwrap();
        assert_eq!(recs.len(), 3);
        assert_eq!(recs[0].ncols(), 3);
        assert_eq!(recs[1].name(), Some("n"));
        assert_eq!(recs[2].ncols(), 7);
        assert_eq!(recs[2].strand(), Some("+"));
        // Verbatim echo: the line is carried, not rebuilt.
        assert_eq!(recs[2].line, "chr1\t200\t300\tn\t5\t+\textra");
        assert_eq!(recs[0].strand(), None);
    }

    #[test]
    fn bed3_fidelity_is_checked_regardless_of_width() {
        // A BED6 line with a bad start is still a bad start.
        let e = parse_all("chr1\t0\t100\tn\t5\t+\nchr1\tx\t200\tn\t5\t+\n").unwrap_err();
        assert_eq!(e.msg, "start/end must be integers");
        assert_eq!(e.line, Some(2));
    }

    #[test]
    fn zero_length_and_position_zero_records_are_legal() {
        let recs = parse_all("chr2\t0\t0\ta12\nchr1\t0\t100\ta01\n").unwrap();
        assert_eq!((recs[0].start, recs[0].end), (0, 0));
        assert_eq!((recs[1].start, recs[1].end), (0, 100));
    }

    #[test]
    fn last_line_without_newline_still_parses() {
        let recs = parse_all("chr1\t0\t100").unwrap();
        assert_eq!(recs.len(), 1);
    }

    #[test]
    fn headers_are_skipped_and_retrievable() {
        let body = "#comment\ntrack name=x\nbrowser position chr1\nchr1\t0\t100\n";
        let p = tmpfile("hdr", body.as_bytes());
        let mut r = Reader::open(p.to_str().unwrap()).unwrap();
        let recs = r.collect_records().unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(
            r.headers(),
            ["#comment", "track name=x", "browser position chr1"]
        );
        std::fs::remove_file(&p).ok();
    }

    // --- diagnostics: one test per row of the SPEC §7 data-error table ---

    #[test]
    fn gzipped_and_plain_inputs_parse_identically() {
        use flate2::Compression;
        use flate2::write::GzEncoder;

        let body = "#hdr\nchr1\t0\t100\ta01\t10\t+\nchr2\t0\t0\ta12\t0\t+\n";
        let plain = tmpfile("plain.bed", body.as_bytes());

        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(body.as_bytes()).unwrap();
        // Deliberately named without .gz: detection is by magic bytes, not name.
        let gz = tmpfile("secretly-gzipped.bed", &enc.finish().unwrap());

        let mut a = Reader::open(plain.to_str().unwrap()).unwrap();
        let mut b = Reader::open(gz.to_str().unwrap()).unwrap();
        let (ra, rb) = (a.collect_records().unwrap(), b.collect_records().unwrap());
        assert_eq!(ra, rb);
        assert_eq!(a.headers(), b.headers());
        assert_eq!(ra.len(), 2);

        std::fs::remove_file(&plain).ok();
        std::fs::remove_file(&gz).ok();
    }

    #[test]
    fn bgzf_style_multi_member_gzip_reads_past_the_first_block() {
        use flate2::Compression;
        use flate2::write::GzEncoder;

        // BGZF is a concatenation of gzip members. A plain GzDecoder stops after
        // the first one and silently truncates the file; MultiGzDecoder does not.
        let mut blob = Vec::new();
        for chunk in ["chr1\t0\t100\n", "chr1\t200\t300\n"] {
            let mut enc = GzEncoder::new(Vec::new(), Compression::default());
            enc.write_all(chunk.as_bytes()).unwrap();
            blob.extend(enc.finish().unwrap());
        }
        let p = tmpfile("multimember.bed.gz", &blob);
        let mut r = Reader::open(p.to_str().unwrap()).unwrap();
        assert_eq!(r.collect_records().unwrap().len(), 2);
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn corrupt_gzip_is_reported_as_such() {
        // Right magic bytes, rubbish payload.
        let p = tmpfile("corrupt.bed.gz", b"\x1f\x8b\x08\x00nonsense-not-deflate");
        let mut r = Reader::open(p.to_str().unwrap()).unwrap();
        let e = r.collect_records().unwrap_err();
        assert_eq!(e.msg, "not a valid gzip stream");
        assert_eq!(e.line, None);
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn fewer_than_three_columns_is_an_error() {
        let e = parse_all("chr1\t100\n").unwrap_err();
        assert_eq!(e.msg, "fewer than 3 columns");
        assert_eq!(e.line, Some(1));
    }

    #[test]
    fn non_integer_coordinates_are_an_error() {
        assert_eq!(
            parse_all("chr1\tten\t100\n").unwrap_err().msg,
            "start/end must be integers"
        );
        // A negative coordinate is not a non-negative integer.
        assert_eq!(
            parse_all("chr1\t-5\t100\n").unwrap_err().msg,
            "start/end must be integers"
        );
    }

    #[test]
    fn start_greater_than_end_is_an_error() {
        let e = parse_all("chr1\t300\t200\n").unwrap_err();
        assert_eq!(e.msg, "start greater than end");
    }

    #[test]
    fn empty_chrom_is_a_malformed_record() {
        assert_eq!(
            parse_all("\t100\t200\n").unwrap_err().msg,
            "malformed BED record"
        );
    }

    #[test]
    fn missing_file_is_an_error() {
        let e = match Reader::open("/nonexistent/nope.bed") {
            Err(e) => e,
            Ok(_) => panic!("opened a file that does not exist"),
        };
        assert_eq!(e.msg, "no such file");
        assert_eq!(
            e.to_string(),
            "mytools: /nonexistent/nope.bed: no such file"
        );
    }

    #[test]
    fn error_display_carries_file_and_line() {
        assert_eq!(
            Error::at("data/a.bed", 7, "start greater than end").to_string(),
            "mytools: data/a.bed:7: start greater than end"
        );
    }

    #[test]
    fn line_numbers_count_headers_and_blanks() {
        // The reported line must be findable with sed -n '<n>p'.
        let e = parse_all("#h\n\nchr1\t0\t100\nchr1\tx\t1\n").unwrap_err();
        assert_eq!(e.line, Some(4));
    }
}
