# intersect cases: SPEC §8 rows 12-17 (issue #8).
#
# bedtools is the oracle: identical arguments to both binaries, stdout and exit
# code compared. Row order is part of the answer, so nothing is sorted before
# the diff. stderr is not part of the contract — our warnings live there.

# Row 12: bare intersect emits the clipped overlap region per overlapping pair.
# Also the chr3-only-in-B case (b.bed has chr3, a.bed does not) and the
# zero-length pairs a07/b07, a12, a16.
check "intersect bare, a.bed vs b.bed"        -- intersect -a "$DATA/a.bed" -b "$DATA/b.bed"

# Row 13: -u emits each overlapping A record once, at full original width.
check "intersect -u, a.bed vs b.bed"          -- intersect -a "$DATA/a.bed" -b "$DATA/b.bed" -u

# Row 14: -v emits the non-overlapping A records, A read from stdin.
check_stdin "intersect -v, a.bed on stdin" "$DATA/a.bed" -- \
  intersect -a - -b "$DATA/b.bed" -v

# Row 15: -wa emits the original A record per pair, B read from a .gz.
check "intersect -wa, b.bed gzipped"          -- intersect -a "$DATA/a.bed" -b "$tmp/b.bed.gz" -wa

# Row 16: -u against genes.bed, whose coordinates are megabase-scale.
check "intersect -u, a.bed vs genes.bed"      -- intersect -a "$DATA/a.bed" -b "$DATA/genes.bed" -u

# Row 17: an empty result is zero bytes and exit 0 (SPEC §5) — both when nothing
# matches and when A is empty. Diffed against the oracle, then asserted directly
# so the "zero bytes" half cannot pass by matching a wrong-but-equal oracle.
check  "intersect -u, no matches at all"      -- intersect -a "$DATA/a.bed" -b "$tmp/nomatch.bed" -u
expect "intersect -u, no matches: empty, exit 0" 0 empty -- \
  intersect -a "$DATA/a.bed" -b "$tmp/nomatch.bed" -u
check  "intersect -v, empty A file"           -- intersect -a "$tmp/empty.bed" -b "$DATA/b.bed" -v
expect "intersect -v, empty A: empty, exit 0" 0 empty -- \
  intersect -a "$tmp/empty.bed" -b "$DATA/b.bed" -v

# Mixed widths interoperate: BED3 A against BED6 B. stdout must stay
# byte-identical to the oracle and the exit code 0; the warning is ours (SPEC
# §8.2) and goes to stderr only.
check "intersect bare, BED3 A vs BED6 B"      -- intersect -a "$tmp/a3.sorted.bed" -b "$DATA/b.bed"
expect "mixed widths still exit 0"             0 nonempty -- \
  intersect -a "$tmp/a3.sorted.bed" -b "$DATA/b.bed"
expect_stderr "mixed widths warn on stderr" "different BED widths" -- \
  intersect -a "$tmp/a3.sorted.bed" -b "$DATA/b.bed"
