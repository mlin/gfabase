#!/bin/bash

set -o pipefail

REPO="$(dirname "$0")/.."
cd "$REPO"
export BASH_TAP_ROOT=test/bash-tap
source test/bash-tap/bash-tap-bootstrap

plan tests 6

aria2c -c -d /tmp -s 10 -x 10 --retry-wait 2 \
    ftp://ftp.dfci.harvard.edu/pub/hli/hifiasm/HG002-trio-v0.11/HG002-v0.11.dip.r_utg.gfa.gz
is "$?" "0" "download HG002-v0.11.dip.r_utg.gfa.gz"

gfabase="cargo run --release --"

export TMPDIR=$(mktemp -d --tmpdir gfabase_hifiasm_test_XXXXXX)

bgzip -dc -@ 4 /tmp/HG002-v0.11.dip.r_utg.gfa.gz \
    | time $gfabase load - "${TMPDIR}/HG002-v0.11.dip.r_utg.gfab"
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

rm -rf "$TMPDIR"
