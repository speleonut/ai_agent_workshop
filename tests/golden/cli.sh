# CLI-level cases: --version and bare invocation must not regress (issue #5).
# Not diffed against bedtools — our exit codes are our own (SPEC §8.1).

expect "--version exits 0"            0 nonempty -- --version
expect "bare mytools is a usage error" 2 empty    -- 
expect "unknown subcommand exits 2"    2 empty    -- frobnicate
expect "sort -r exits 2"               2 empty    -- sort -i "$DATA/a.bed" -r
expect "intersect -u -v exits 2"       2 empty    -- intersect -a "$DATA/a.bed" -b "$DATA/b.bed" -u -v

expect_stderr "usage error names the option" "unrecognised option '-r'" -- sort -i "$DATA/a.bed" -r
expect_stderr "mutually exclusive message"   "-u and -v are mutually exclusive" -- \
  intersect -a "$DATA/a.bed" -b "$DATA/b.bed" -u -v
