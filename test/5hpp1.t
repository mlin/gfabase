#!/bin/bash

set -o pipefail

REPO="$(dirname "$0")/.."
cd "$REPO"
export BASH_TAP_ROOT=test/bash-tap
source test/bash-tap/bash-tap-bootstrap
export LC_ALL=C

plan tests 14

aria2c -c -d /tmp -s 16 -x 16 --retry-wait 2 ftp://ftp.dfci.harvard.edu/pub/hli/minigraph/HPP/GRCh38-freeze1.gfa.gz
is "$?" "0" "download gfa"

gfabase="cargo run --release -- --verbose"

export TMPDIR=$(mktemp -d --tmpdir gfabase_hpp1_test_XXXXXX)

bgzip -@ 4 -dc /tmp/GRCh38-freeze1.gfa.gz | time $gfabase load -o "${TMPDIR}/GRCh38-freeze1.gfab" --compress 1
is "$?" "0" "gfabase load"

time $gfabase view "${TMPDIR}/GRCh38-freeze1.gfab" | wc -c
is "$?" "0" "gfabase view"

ls -lh "${TMPDIR}/GRCh38-freeze1.gfab"

gsql() {
    genomicsqlite "$1" "$2"| tail -n 1
}
export -f gsql

is "$(gsql "${TMPDIR}/GRCh38-freeze1.gfab" 'pragma application_id')"  "1734762850" "gfab application_id"
is "$(gsql "${TMPDIR}/GRCh38-freeze1.gfab" 'select count(1) from gfa1_segment_meta')"  "236294" "gfab segment count"
is "$(gsql "${TMPDIR}/GRCh38-freeze1.gfab" 'select count(1) from gfa1_link')" "348614" "gfab link count"
is "$(gsql "${TMPDIR}/GRCh38-freeze1.gfab" 'select sum(sequence_length) from gfa1_segment_meta')" "3216020942" "gfab base metacount"
is "$(gsql "${TMPDIR}/GRCh38-freeze1.gfab" 'select sum(twobit_length(sequence_twobit)) from gfa1_segment_sequence')" "3216020942" "gfab base count"

$gfabase sub "${TMPDIR}/GRCh38-freeze1.gfab" -o "${TMPDIR}/CR1.gfab" --biconnected 2 215577
is "$(gsql "${TMPDIR}/CR1.gfab" 'select count(1) from gfa1_segment_meta')" "16" "CR1 --cutpoints 2 segment count"
is "$(gsql "${TMPDIR}/CR1.gfab" 'select count(1) from gfa1_link')" "22" "CR1 --cutpoints 3 link count"
$gfabase sub "${TMPDIR}/GRCh38-freeze1.gfab" -o "${TMPDIR}/CR1.3.gfab" --biconnected 4 215577
is "$(gsql "${TMPDIR}/CR1.3.gfab" 'select count(1) from gfa1_segment_meta')" "29" "CR1 --cutpoints 4 segment count"
is "$(gsql "${TMPDIR}/CR1.3.gfab" 'select count(1) from gfa1_link')" "41" "CR1 --cutpoints 3 link count"

$gfabase sub "${TMPDIR}/GRCh38-freeze1.gfab" -o "${TMPDIR}/RHD.gfab" --biconnected 2 --range chr1:25,272,509-25,330,445
is "$(gsql "${TMPDIR}/RHD.gfab" 'select count(1) from gfa1_segment_meta')" "17" "RHD --cutpoints 2 segment count"
is "$(gsql "${TMPDIR}/RHD.gfab" 'select count(1) from gfa1_link')" "25" "RHD --cutpoints 3 link count"

rm -rf "$TMPDIR"
