#!/bin/bash

set -o pipefail

REPO="$(dirname "$0")/.."
cd "$REPO"
export BASH_TAP_ROOT=test/bash-tap
source test/bash-tap/bash-tap-bootstrap

plan tests 18

cargo build --release
is "$?" "0" "cargo build"
PATH="$(pwd)/target/release:${PATH}"
gfabase version
is "$?" "0" "gfabase version"

export TMPDIR=$(mktemp -d --tmpdir gfabase_metaspades_test_XXXXXX)

# extract metaspades GFA
zstd -dc test/data/atcc_staggered.assembly_graph_with_scaffolds.gfa.zst > "${TMPDIR}/atcc_staggered.assembly_graph_with_scaffolds.gfa"

# roundtrip it and check fidelity
gfabase load --compress 1 "${TMPDIR}/atcc_staggered.assembly_graph_with_scaffolds.gfa" "${TMPDIR}/atcc_staggered.assembly_graph_with_scaffolds.gfab"
is "$?" "0" "gfabase load"

time gfabase view "${TMPDIR}/atcc_staggered.assembly_graph_with_scaffolds.gfab" > "${TMPDIR}/atcc_staggered.assembly_graph_with_scaffolds.roundtrip.gfa"
is "$?" "0" "gfabase view"

export LC_ALL=C
seq_digest=$(cat "${TMPDIR}/atcc_staggered.assembly_graph_with_scaffolds.gfa" | grep "^S" | cut -f1-3 | sort -V | sha256sum)
roundtrip_seq_digest=$(cat "${TMPDIR}/atcc_staggered.assembly_graph_with_scaffolds.roundtrip.gfa" | grep "^S" | cut -f1-3 | sort -V | sha256sum)
is "$roundtrip_seq_digest" "$seq_digest" "sequences roundtrip identical"

link_digest=$(cat "${TMPDIR}/atcc_staggered.assembly_graph_with_scaffolds.gfa" | grep "^L" | cut -f1-6 | sort -V | sha256sum)
roundtrip_link_digest=$(cat "${TMPDIR}/atcc_staggered.assembly_graph_with_scaffolds.roundtrip.gfa" | grep "^L" | cut -f1-6 | sort -V | sha256sum)
is "$roundtrip_link_digest" "$link_digest" "links roundtrip identical"

path_digest=$(cat "${TMPDIR}/atcc_staggered.assembly_graph_with_scaffolds.gfa" | grep "^P" | cut -f1-4 | sort -V | sha256sum)
roundtrip_path_digest=$(cat "${TMPDIR}/atcc_staggered.assembly_graph_with_scaffolds.roundtrip.gfa" | grep "^P" | cut -f1-4 | sort -V | sha256sum)
is "$roundtrip_path_digest" "$path_digest" "paths roundtrip identical"

# additional roundtrip back through load
gfabase load "${TMPDIR}/atcc_staggered.assembly_graph_with_scaffolds.roundtrip.gfa" "${TMPDIR}/roundtrip2.gfab"
is "$?" "0" "roundtrip2"

gfabase view "${TMPDIR}/roundtrip2.gfab" --no-sequences | grep GATTACA
is "$?" "1" "view --no-sequences"
gfabase view "${TMPDIR}/roundtrip2.gfab" | grep GATTACA > /dev/null
is "$?" "0" "view --no-sequences (control)"

# sub two scaffolds and make sure we get those Paths
time gfabase sub "${TMPDIR}/atcc_staggered.assembly_graph_with_scaffolds.gfab" "${TMPDIR}/sub.gfab" \
    20412 106423 14364 133587 133589 17280
is "$?" "0" "sub scaffolds"
gfabase view "${TMPDIR}/sub.gfab" | grep NODE_2_length_747618_cov_15.708553_3
is "$?" "0" "sub scaffold NODE_2_length_747618_cov_15.708553_3"
gfabase view "${TMPDIR}/sub.gfab" | grep NODE_2_length_747618_cov_15.708553_4
is "$?" "0" "sub scaffold grep NODE_2_length_747618_cov_15.708553_4"

time gfabase sub --view "${TMPDIR}/atcc_staggered.assembly_graph_with_scaffolds.gfab" - \
    20412 106423 14364 133587 133589 17280 > "${TMPDIR}/sub.gfa"
is "$?" "0" "sub --view scaffolds"
grep NODE_2_length_747618_cov_15.708553_3 "${TMPDIR}/sub.gfa"
is "$?" "0" "sub --view scaffold NODE_2_length_747618_cov_15.708553_3"

# sub by path
time gfabase sub --view "${TMPDIR}/atcc_staggered.assembly_graph_with_scaffolds.gfab" "${TMPDIR}/sub_by_path.gfa" --path \
    NODE_2_length_747618_cov_15.708553_3 NODE_2_length_747618_cov_15.708553_4
is "$?" "0" "sub --view by path"
is "$(cat "${TMPDIR}/sub_by_path.gfa" | wc -l)" "14" "sub --view by path line count"

# test behavior w/ empty input
gfabase load /dev/null "${TMPDIR}/empty.gfab"
is "$?" "3" "gfabase load empty"

rm -rf "$TMPDIR"
