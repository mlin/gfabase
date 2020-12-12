#!/bin/bash

set -o pipefail

REPO="$(dirname "$0")/.."
cd "$REPO"
export BASH_TAP_ROOT=test/bash-tap
source test/bash-tap/bash-tap-bootstrap

plan tests 7

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

rm -rf "$TMPDIR"
