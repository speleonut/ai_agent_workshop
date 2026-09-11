# SPEC.md — mytools

A small reimplementation of a subset of bedtools. This file is the contract: it is
what the golden tests encode and what `CLAUDE.md` defers to on design questions.

**Tiebreaker rule.** Where this spec is silent, ambiguous, or in conflict with
observed behaviour, **match real `bedtools`**. The oracle wins. Deviations are legal
only where §8 lists them explicitly, with a reason.

---

## 1. Scope

**Subcommands in v1:** `sort`, `merge`, `intersect`.

**Explicitly NOT in v1:** `subtract`, `closest`. Also out: every bedtools subcommand
not named above (`window`, `flank`, `slop`, `genomecov`, `map`, …).

Rationale: three subcommands finished and oracle-clean beats five half-built. Adding a
fourth is a new decision, not a natural next step.

## 2. Input formats

- **Formats accepted:** BED3 and wider — BED3, BED4, BED5, BED6, and beyond. Columns
  past the ones a subcommand needs are carried or ignored per that subcommand's rules
  (§5), never rejected for existing.
- **BED3 fidelity is always checked.** The first three columns must be present and
  well-formed on every data line: `chrom` non-empty, `start` and `end` non-negative
  integers, `start <= end`. This check runs regardless of how wide the record is.
- **Mixed widths may interoperate.** A BED3 file and a BED6 file can be used together
  in one `intersect`. Where the operation needs a column the narrower file lacks, a
  warning is emitted (§7) and the oracle's behaviour is reproduced exactly — including
  when that behaviour is to return nothing (see §3, `merge -s` on BED3).
- **Input sources:** a file argument or stdin. `-` means stdin. (From `CLAUDE.md`.)
- **Compressed input:** `.gz` and bgzipped (BGZF) files are read transparently,
  detected by magic bytes rather than filename. BGZF is a valid gzip stream, so
  sequential reading needs one code path, not two. Indexed random access (`.tbi`,
  `.csi`) is **not** supported — irrelevant to streaming.
  Real bedtools reads `.gz` natively, so this is required for oracle parity, not a
  convenience. Implemented with `flate2`; see §9.
- **`track`, `browser`, and `#` lines are skipped** and do not reach the interval
  logic. They may be reproduced on output under `-header` (§5).

## 3. Interval semantics

- **Coordinate system: 0-based, half-open.** `chr1 100 200` covers bases 100..199.
  Not a decision — BED says so.
- **Overlap predicate:** `a.start < b.end && b.start < a.end`. Note the strict `<`.
  Every off-by-one in this project lives in that comparison.
- **Bookended intervals** (`a.end == b.start`) do **not** overlap. They **do** merge
  under `merge -d 0`.
- **Zero-length intervals** (`start == end`) are legal and appear in the fixtures.
  Their overlap behaviour is whatever `bedtools` does — do not reason from first
  principles, and do not "fix" them. Encode the oracle, comment that it is oracle
  behaviour, move on.
- **Minimum overlap:** one base pair. No fractional-overlap flags in v1, so there is
  no threshold to configure.
- **Strand on strandless input:** `merge -s` / `-S` against BED3 produces **zero bytes
  on stdout and exit 0** — verified against bedtools v2.31.1. It is not an error. A
  warning is emitted (§7); stdout stays byte-identical to the oracle.

## 4. Flags per subcommand

Flag names and meanings match bedtools exactly. Any flag not listed here is a usage
error (§7) — never silently ignored.

| Subcommand  | Flags in v1              | Notes |
|-------------|--------------------------|-------|
| `sort`      | `-i`, `-header`          | Ordering is fixed and ties keep input order; see below |
| `merge`     | `-i`, `-d`, `-s`, `-S`, `-header` | Requires sorted input |
| `intersect` | `-u`, `-v`, `-wa`, `-header` | `-a` and `-b` name the two inputs |

Global: `--version` (exit 0), `--quiet` (silence warnings, §7).

**`sort`**
- Default and only ordering: **chromosome lexicographically, then start ascending.**
  `chr10` sorts before `chr2`. This is bedtools' actual behaviour, verified —
  natural/numeric ordering was considered and **rejected** because it breaks oracle
  parity on any genome with `chr10`-style names. `a.bed` has only `chr1`, `chr2`,
  `chrX` and does not expose this, but `genes.bed` does: `chr17` sorts before `chr7`.
- **Ties on `(chrom, start)` keep input order. The sort is stable and `end` is never
  compared.** Verified against bedtools v2.31.1: given `chr1 100 900`, `chr1 100 200`,
  `chr1 100 500` in that input order, bedtools emits them in that same order, not by
  ascending end. An earlier draft of this spec said "then end ascending" — that was
  wrong, and the tiebreaker rule at the top of this file applies: the oracle wins.
  No fixture exposes it (`a.bed`, `b.bed` and `genes.bed` sort identically under
  either comparator), so this is pinned by a **unit** test, not a golden one.
- **Input is named by `-i`:** `mytools sort -i <file>`, and `-i -` for stdin. bedtools
  rejects a positional filename (`Unrecognized parameter`), and golden tests pass
  identical arguments to both binaries, so `-i` is required for parity rather than a
  style choice. `merge` takes `-i` likewise; `intersect` takes `-a` and `-b`.
- No `-r`. `bedtools sort` has no such flag; passing it is a usage error.
- `-g` / `-faidx` (explicit chromosome order from a file) are **not** in v1.
- Not in v1: `-sizeA`, `-sizeD`, `-chrThenSizeA/D`, `-chrThenScoreA/D`.

**`merge`**
- `-d N`: maximum distance between features still merged. Default `0` — overlapping
  and bookended features merge. **Negative values require that many bp of overlap.**
- `-s`: merge only features on the same strand.
- `-S +` or `-S -`: merge only features on the named strand.
- Both `-s` and `-S` read column 6; see §3 for strandless input.
- **Requires pre-sorted input.** `merge` does not sort for you. Unsorted input is a
  data error (§7), matching bedtools.
- Not in v1: `-c` / `-o` column aggregation (18 operations plus broadcasting rules —
  deliberately cut as the single largest chunk of surface in `merge`).

**`intersect`**
- Default, no flag: report the **overlapping portion** of each A/B pair.
- `-u`: report each A feature **once** if it overlaps any B feature.
- `-v`: report A features with **no** overlap in B.
- `-wa`: report the **original A record**, not the overlapping portion.
- `-u` and `-v` are mutually exclusive; together they are a usage error.
- Does **not** require sorted input.
- Not in v1: `-wb`, `-wo`, `-wao`, `-c`, `-loj`, `-f`, `-F`, `-r`, `-e`, `-s`, `-S`.

## 5. Output

- **Separator:** a single tab. **Line ending:** `\n`, including on the final line.
- **Empty results print nothing** — zero bytes on stdout — and **exit 0**. This
  covers no-match results and empty input files alike. Verified against the oracle.
  A literal "no output" marker was considered and **rejected**: it would break the
  stdout diff on every empty case and corrupt downstream pipes.
- **stdout is data.** Warnings, errors, and diagnostics go to stderr, always.
- **Column width per subcommand:**
  - `sort` — echoes input records unchanged, full width preserved. Records are
    replayed as the original line bytes, so trailing columns survive untouched.
  - `merge` — BED3 (`chrom`, `start`, `end`). Without `-c`/`-o` (not in v1), merged
    records carry no name, score, or strand.
  - `intersect` — `-wa` and `-v` emit the full original A record; `-u` emits the full
    original A record once; the default emits the clipped overlap region.
- **`-header`:** when given, the header lines skipped on input (`#`, `track`,
  `browser`) are reproduced ahead of the results, taken from the A file. Without it
  they never appear.

## 6. Memory model

- **Streaming by default.** Read line by line, write as you go.
- **Target scale: exome, ~200k intervals.** This is the size `mytools` promises to
  handle comfortably. Whole-genome inputs in the tens of millions are out of scope —
  use real bedtools.
- **Per subcommand:**
  - `sort` — cannot stream. Holds **one chromosome** in memory at a time, the single
    exception `CLAUDE.md` allows.
  - `merge` — streams. Needs only the current run of merge-candidates, because input
    is required to be sorted already.
  - `intersect` — holds **all of B** in memory and streams A against it. B's size, not
    A's, sets the ceiling. At the ~200k target this is a few tens of MB.

## 7. Errors, warnings, and exit codes

**Exit codes** (from `CLAUDE.md`; see §8 for the deviation this creates):

| Situation                        | stderr message                                          | exit |
|----------------------------------|---------------------------------------------------------|------|
| Success                          | —                                                        | `0`  |
| Success, empty result            | —                                                        | `0`  |
| Malformed BED line               | `mytools: <file>:<line>: malformed BED record`           | `1`  |
| Non-integer start/end            | `mytools: <file>:<line>: start/end must be integers`     | `1`  |
| `start > end`                    | `mytools: <file>:<line>: start greater than end`         | `1`  |
| Fewer than 3 columns             | `mytools: <file>:<line>: fewer than 3 columns`           | `1`  |
| Missing input file               | `mytools: <file>: no such file`                          | `1`  |
| Unreadable / corrupt gzip        | `mytools: <file>: not a valid gzip stream`               | `1`  |
| Unsorted input to `merge`        | `mytools: <file>:<line>: input is not sorted`            | `1`  |
| Unknown flag                     | `mytools: unrecognised option '<flag>'` + usage          | `2`  |
| Missing required argument        | `mytools: <flag> requires an argument` + usage           | `2`  |
| Mutually exclusive flags         | `mytools: -u and -v are mutually exclusive` + usage      | `2`  |
| No arguments                     | usage                                                    | `2`  |

All messages are prefixed `mytools:` and carry file and line number where one exists.

**Warnings** are diagnostics that do not stop the run:

- Written to **stderr**, prefixed `mytools: warning:`.
- **Do not change the exit code** — a run that only warns still exits `0`.
- **Do not touch stdout**, so they are invisible to golden-test stdout diffs.
- Silenced by `--quiet`.
- Emitted when: a strand flag is used against input with no strand column; files of
  differing BED widths are combined in one operation.

## 8. Correctness

**Oracle:** real `bedtools` on the files in `data/`. Non-negotiable. If our output
differs from bedtools on the same input, we are wrong.

**Golden coverage in v1 — "each thing once."** Every flag is exercised at least once
and every input source at least once, without testing every pairing of the two. The
residual risk of untested pairings is accepted deliberately.

| # | Case                                        | Input source |
|---|---------------------------------------------|--------------|
| 1 | `sort` `a.bed`                              | file         |
| 2 | `sort` `b.bed`                              | stdin (`-`)  |
| 3 | `sort` `genes.bed`                          | `.gz`        |
| 4 | `sort -header` on a file with header lines  | file         |
| 5 | `merge -d 0` on sorted `a.bed`              | file         |
| 6 | `merge -d 10` on sorted `a.bed`             | stdin (`-`)  |
| 7 | `merge -d -5` (negative: requires overlap)  | file         |
| 8 | `merge -s` on sorted `a.bed` (BED6)         | file         |
| 9 | `merge -S +` on sorted `a.bed`              | file         |
|10 | `merge -S -` on sorted `a.bed`              | `.gz`        |
|11 | `merge -s` on sorted BED3 (strandless)      | file         |
|12 | `intersect` bare, `a.bed` vs `b.bed`        | file         |
|13 | `intersect -u`, `a.bed` vs `b.bed`          | file         |
|14 | `intersect -v`, `a.bed` vs `b.bed`          | stdin (`-`)  |
|15 | `intersect -wa`, `a.bed` vs `b.bed`         | `.gz`        |
|16 | `intersect -u`, `a.bed` vs `genes.bed`      | file         |
|17 | `intersect -v` with no matches (empty out)  | file         |
|18 | any subcommand on an empty input file       | file         |

Each case compares **stdout and exit code** against bedtools. `merge` and `sort` cases
requiring sorted input sort into a temp file first — `a.bed` and `b.bed` are
deliberately unsorted.

**Error-path cases are ours alone.** Because exit codes deviate (below), error cases
assert `mytools`' own documented codes from §7 rather than diffing against bedtools.

**Known deviations from bedtools, accepted:**

1. **Usage errors exit `2`, not `1`.** bedtools exits `1` for everything, including
   unrecognised flags. `CLAUDE.md` reserves `2` for usage errors, `mytools --version`
   already ships that behaviour, and `2` is the conventional CLI choice. Cost: error
   paths cannot be exit-code-compared against the oracle.
2. **Warnings on stderr that bedtools does not emit** (strand flags on strandless
   input, mixed BED widths). stdout is unaffected, so golden diffs are unaffected.
3. **`--quiet` is an extension.** No bedtools equivalent; covered by unit tests only.

## 9. Language and layout

- **Implementation language: Rust**, edition 2024. Every subcommand and every test.
- **Dependencies:** standard library only, with one exception — `flate2` for
  gzip/bgzip input (§2). Any further dependency needs an explicit decision.
- **Entry point:** `src/main.rs`, binary `mytools`, invoked as
  `mytools <subcommand> [flags]`.
- **Tests:** unit tests in `#[cfg(test)]` modules run by `cargo test`, requiring no
  bedtools. Golden tests in `tests/run_golden.sh`, which shells out to the built
  binary and to real bedtools.
- **Linter:** `cargo clippy -- -D warnings`, plus `cargo fmt --check`.

---

## Open assumptions

Recorded rather than silently decided. Correct any of these and the spec follows.

- `-header` is offered on all three subcommands, as bedtools does.
- "Strand aware" for `merge` was read as both `-s` and `-S`, the pair bedtools ships.
- `sort` takes no flags beyond `-header`; `-g`/`-faidx` were not requested.
- Warnings go to stderr, leave the exit code at `0`, and `--quiet` is the silencing
  flag — the mechanism was not specified beyond "can be silenced".
