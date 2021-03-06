#!/bin/bash

set -o pipefail

REPO="$(dirname "$0")/.."
cd "$REPO"
export BASH_TAP_ROOT=test/bash-tap
source test/bash-tap/bash-tap-bootstrap

plan tests 7

aria2c -c -d /tmp -s 10 -x 10 --retry-wait 2 \
    ftp://ftp.dfci.harvard.edu/pub/hli/hifiasm/HG002-trio-v0.11/HG002-v0.11.dip.r_utg.gfa.gz
is "$?" "0" "download HG002-v0.11.dip.r_utg.gfa.gz"

gfabase="cargo run --release -- --verbose"

export TMPDIR=$(mktemp -d --tmpdir gfabase_hifiasm_test_XXXXXX)

bgzip -dc -@ 4 /tmp/HG002-v0.11.dip.r_utg.gfa.gz \
    | time $gfabase load -o "${TMPDIR}/HG002-v0.11.dip.r_utg.gfab" --compress 1
is "$?" "0" "gfabase load"

ls -lh /tmp/HG002-v0.11.dip.r_utg.gfa.gz "${TMPDIR}/HG002-v0.11.dip.r_utg.gfab"

gsql() {
    genomicsqlite "${TMPDIR}/HG002-v0.11.dip.r_utg.gfab" "$1"| tail -n 1
}
export -f gsql

is "$(gsql 'select count(1) from gfa1_segment_meta')"  "53405" "gfab segment count"
is "$(gsql 'select count(1) from gfa1_link')" "145558" "gfab link count"
is "$(gsql 'select sum(sequence_length) from gfa1_segment_meta')" "6506642704" "gfab base metacount"
is "$(gsql 'select sum(twobit_length(sequence_twobit)) from gfa1_segment_sequence')" "6506642704" "gfab base count"

# sub by segment name
$gfabase sub --view "${TMPDIR}/HG002-v0.11.dip.r_utg.gfab" utg000042l utg050830l utg021888l > "${TMPDIR}/sub.gfa"
is "$(cat "${TMPDIR}/sub.gfa" | wc -l)" "8" "sub by segment name"

rm -rf "$TMPDIR"
