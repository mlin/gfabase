#!/bin/bash

set -o pipefail

REPO="$(dirname "$0")/.."
cd "$REPO"
export BASH_TAP_ROOT=test/bash-tap
source test/bash-tap/bash-tap-bootstrap
export LC_ALL=C

plan tests 7

aria2c -c -d /tmp -o shasta-HG002-Guppy-3.6.0-run4-UL.gfa -s 16 -x 16 --retry-wait 2 \
    https://s3-us-west-2.amazonaws.com/czi.paolo-public/HG002-Guppy-3.6.0-run4-UL/Assembly.gfa
is "$?" "0" "download gfa"

gfabase="cargo run --release -- --verbose"

export TMPDIR=$(mktemp -d --tmpdir gfabase_shasta_test_XXXXXX)

time $gfabase load /tmp/shasta-HG002-Guppy-3.6.0-run4-UL.gfa "${TMPDIR}/Assembly.gfab" --compress 1
is "$?" "0" "gfabase load"

time $gfabase view "${TMPDIR}/Assembly.gfab" | wc -c
is "$?" "0" "gfabase view"

ls -lh "${TMPDIR}/Assembly.gfab"

gsql() {
    genomicsqlite "${TMPDIR}/Assembly.gfab" "$1"| tail -n 1
}
export -f gsql

is "$(gsql 'select count(1) from gfa1_segment_meta')"  "12520" "gfab segment count"
is "$(gsql 'select count(1) from gfa1_link')" "18456" "gfab link count"
is "$(gsql 'select sum(sequence_length) from gfa1_segment_meta')" "5100857444" "gfab base metacount"
is "$(gsql 'select sum(twobit_length(sequence_twobit)) from gfa1_segment_sequence')" "5100857444" "gfab base count"

rm -rf "$TMPDIR"
