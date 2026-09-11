# `sort` cases: SPEC §8 rows 1-4 and 18. Diffed against real bedtools — stdout
# bytes and exit code, in row order, never re-sorted before the diff.
#
# The tie rule (records equal on chrom and start keep input order, `end` is
# never compared) has no golden case on purpose: a.bed, b.bed and genes.bed all
# sort identically under either comparator, so it is pinned by unit test in
# src/sort.rs instead.
#
# Fixtures in $tmp are built by run_golden.sh: genes.bed.gz, empty.bed and
# header.bed (a.bed behind #/track/browser lines).

check "sort a.bed"                -- sort -i "$DATA/a.bed"
check_stdin "sort b.bed on stdin" "$DATA/b.bed" -- sort -i -
# genes.bed is the only fixture with chr10-style names: it proves chr17 sorts
# before chr7, which a.bed (chr1/chr2/chrX) cannot show.
check "sort genes.bed.gz"         -- sort -i "$tmp/genes.bed.gz"
check "sort -header header.bed"   -- sort -header -i "$tmp/header.bed"
check "sort empty file"           -- sort -i "$tmp/empty.bed"
