#!/usr/bin/env bash
# Golden tests: diff mytools against real bedtools.
#
#   ./tests/run_golden.sh            run every case
#   ./tests/run_golden.sh sort       run only cases from tests/golden/sort.sh
#
# bedtools is the oracle. If we differ from it on the same input, we are wrong.
# Error paths are the exception: our exit codes deviate deliberately (SPEC §8.1),
# so those are asserted against SPEC, not diffed against bedtools.
set -uo pipefail

here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/.." && pwd)
DATA=$root/data

# Build unless the caller supplied a binary to test.
if [[ -z ${MYTOOLS:-} ]]; then
  cargo build --quiet --manifest-path "$root/Cargo.toml" || exit 1
  MYTOOLS=$root/target/debug/mytools
fi
export MYTOOLS DATA

command -v bedtools >/dev/null || { echo "bedtools not installed; golden tests cannot run"; exit 1; }

tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
export tmp
pass=0; fail=0

# ---------------------------------------------------------------- fixtures
# Built once here rather than in each case file. a.bed and b.bed are
# deliberately unsorted, and merge needs sorted input.
bedtools sort -i "$DATA/a.bed" > "$tmp/a.sorted.bed"
bedtools sort -i "$DATA/b.bed" > "$tmp/b.sorted.bed"
cut -f1-3 "$tmp/a.sorted.bed" > "$tmp/a3.sorted.bed"          # BED3: strandless
gzip -c "$tmp/a.sorted.bed" > "$tmp/a.sorted.bed.gz"
gzip -c "$DATA/genes.bed"   > "$tmp/genes.bed.gz"
gzip -c "$DATA/b.bed"       > "$tmp/b.bed.gz"
: > "$tmp/empty.bed"
{ echo "#comment line"
  echo "track name=fixture description=\"header case\""
  echo "browser position chr1:1-1000"
  cat "$DATA/a.bed"
} > "$tmp/header.bed"
# A file whose intervals cannot overlap a.bed: for the empty-result cases.
printf 'chr9\t10\t20\tno1\t0\t+\nchr9\t30\t40\tno2\t0\t-\n' > "$tmp/nomatch.bed"

# ---------------------------------------------------------------- helpers
#
# check <name> -- <args...>
#   Run mytools and bedtools with identical arguments; diff stdout and exit code.
check() {
  local name=$1; shift; shift                 # drop the literal --
  "$MYTOOLS" "$@" > "$tmp/got"  2>"$tmp/got.err"; local got_rc=$?
  bedtools    "$@" > "$tmp/want" 2>/dev/null;     local want_rc=$?
  _verdict "$name" "$got_rc" "$want_rc"
}

# check_stdin <name> <file> -- <args...>
#   Same, but <file> is fed to both binaries on stdin. Use with `-i -`, `-a -`.
check_stdin() {
  local name=$1 infile=$2; shift 2; shift
  "$MYTOOLS" "$@" < "$infile" > "$tmp/got"  2>"$tmp/got.err"; local got_rc=$?
  bedtools    "$@" < "$infile" > "$tmp/want" 2>/dev/null;     local want_rc=$?
  _verdict "$name" "$got_rc" "$want_rc"
}

_verdict() {
  local name=$1 got_rc=$2 want_rc=$3
  if [[ $got_rc -ne $want_rc ]]; then
    echo "FAIL $name (exit $got_rc, bedtools gave $want_rc)"
    sed 's/^/      /' "$tmp/got.err" | head -3
    (( fail++ )); return
  fi
  if diff -q "$tmp/want" "$tmp/got" >/dev/null; then
    echo "ok   $name"; (( pass++ ))
  else
    echo "FAIL $name"
    diff -u "$tmp/want" "$tmp/got" | sed 's/^/      /' | head -20
    (( fail++ ))
  fi
}

# expect <name> <exit-code> <empty|nonempty|any> -- <args...>
#   mytools alone, against SPEC rather than the oracle. For error paths,
#   --quiet, and anything else where bedtools' behaviour is its own.
expect() {
  local name=$1 want_rc=$2 shape=$3; shift 3; shift
  "$MYTOOLS" "$@" > "$tmp/got" 2>"$tmp/got.err"; local got_rc=$?
  if [[ $got_rc -ne $want_rc ]]; then
    echo "FAIL $name (exit $got_rc, expected $want_rc)"
    sed 's/^/      /' "$tmp/got.err" | head -3
    (( fail++ )); return
  fi
  case $shape in
    empty)    [[ -s $tmp/got ]] && { echo "FAIL $name (stdout not empty)"; (( fail++ )); return; } ;;
    nonempty) [[ -s $tmp/got ]] || { echo "FAIL $name (stdout empty)";     (( fail++ )); return; } ;;
  esac
  echo "ok   $name"; (( pass++ ))
}

# expect_stderr <name> <substring> -- <args...>
#   Assert our own diagnostic wording (SPEC §7). Never used against bedtools.
expect_stderr() {
  local name=$1 want=$2; shift 2; shift
  "$MYTOOLS" "$@" > "$tmp/got" 2>"$tmp/got.err"
  if grep -qF -- "$want" "$tmp/got.err"; then
    echo "ok   $name"; (( pass++ ))
  else
    echo "FAIL $name (stderr did not contain: $want)"
    sed 's/^/      /' "$tmp/got.err" | head -3
    (( fail++ ))
  fi
}

# ---------------------------------------------------------------- cases
# One file per subcommand in tests/golden/, so they can be written independently.
only=${1:-}
for casefile in "$here"/golden/*.sh; do
  [[ -e $casefile ]] || continue
  [[ -n $only && $(basename "$casefile" .sh) != "$only" ]] && continue
  echo "== $(basename "$casefile" .sh)"
  # shellcheck source=/dev/null
  source "$casefile"
done

echo "---"
echo "$pass passed, $fail failed"
[[ $fail -eq 0 ]]
