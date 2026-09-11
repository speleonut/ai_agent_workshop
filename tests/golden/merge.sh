# merge cases (SPEC §8, rows 5-11), plus this subcommand's own error paths.
#
# merge requires sorted input and data/a.bed is deliberately unsorted, so every
# oracle case runs on the pre-sorted copies the harness built in $tmp.

# A sorted file with a leading header block, for the -header case: $tmp/header.bed
# is headers over the *unsorted* a.bed, which merge rejects.
{ echo "#merge header case"
  echo "track name=merge description=\"header case\""
  cat "$tmp/a.sorted.bed"
} > "$tmp/merge.header.bed"

# ---- oracle: stdout and exit code diffed against real bedtools ----------
check "merge -d 0 on sorted a.bed"        -- merge -i "$tmp/a.sorted.bed"
check_stdin "merge -d 10 from stdin" "$tmp/a.sorted.bed" -- merge -d 10 -i -
check "merge -d -5 requires overlap"      -- merge -d -5 -i "$tmp/a.sorted.bed"
check "merge -s on sorted a.bed"          -- merge -s -i "$tmp/a.sorted.bed"
check "merge -S + on sorted a.bed"        -- merge -S + -i "$tmp/a.sorted.bed"
check "merge -S - on sorted a.bed.gz"     -- merge -S - -i "$tmp/a.sorted.bed.gz"
check "merge -s on sorted BED3"           -- merge -s -i "$tmp/a3.sorted.bed"
check "merge on an empty file"            -- merge -i "$tmp/empty.bed"
check "merge -header reproduces headers"  -- merge -header -i "$tmp/merge.header.bed"

# ---- ours alone: exit codes deviate (SPEC §8.1), so assert against SPEC ----
expect "merge -s on BED3 is empty and exits 0" 0 empty -- merge -s -i "$tmp/a3.sorted.bed"
expect_stderr "merge -s on BED3 warns" "no strand column" -- merge -s -i "$tmp/a3.sorted.bed"

# --quiet silences that warning without touching stdout or the exit code.
"$MYTOOLS" merge --quiet -s -i "$tmp/a3.sorted.bed" > "$tmp/got" 2>"$tmp/got.err"
if [[ -s $tmp/got.err ]]; then
  echo "FAIL merge --quiet silences the strandless warning"
  sed 's/^/      /' "$tmp/got.err" | head -3
  (( fail++ ))
else
  echo "ok   merge --quiet silences the strandless warning"; (( pass++ ))
fi

expect "merge on unsorted input exits 1"  1 any   -- merge -i "$DATA/a.bed"
expect_stderr "unsorted input is named as such" "input is not sorted" -- merge -i "$DATA/a.bed"
expect "merge -c is a usage error"        2 empty -- merge -i "$tmp/a.sorted.bed" -c 4
expect "merge -o is a usage error"        2 empty -- merge -i "$tmp/a.sorted.bed" -o collapse
