#!/bin/bash

# bash-tap tests for gfabase CLI ops
# requires cargo toolchain and:
#   apt install -y zstd && pip3 install --upgrade genomicsqlite

set -o pipefail

REPO="$(dirname "$0")/.."
cd "$REPO"
export BASH_TAP_ROOT=test/bash-tap
source test/bash-tap/bash-tap-bootstrap

plan tests 15

cargo build --release
is "$?" "0" "cargo build"
PATH="$(pwd)/target/release:${PATH}"
gfabase version
is "$?" "0" "gfabase version"

export TMPDIR=$(mktemp -d --tmpdir gfabase_cli_test_XXXXXX)

# extract test rGFA
zstd -dc test/data/GRCh38-20-0.10b.chr22_chrY.gfa.zst > "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfa"

# roundtrip it and check fidelity
gfabase load --compress 1 "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfa" "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfab"
is "$?" "0" "gfabase load"
cat "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfa" \
    | gfabase --verbose load --compress 1 --topology - "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfab"
is "$?" "0" "gfabase load stdin"

gfabase view "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfab" > "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.roundtrip.gfa"
is "$?" "0" "gfabase view"

seq_digest=$(cat "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfa" | grep "^S" | cut -f3 | LC_ALL=C sort | sha256sum)
roundtrip_seq_digest=$(cat "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.roundtrip.gfa" | grep "^S" | cut -f3 | LC_ALL=C sort | sha256sum)
is "$roundtrip_seq_digest" "$seq_digest" "sequences roundtrip identical"

is `cat "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfa" | grep "^L" | wc -l` \
   `cat "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.roundtrip.gfa" | grep "^L" | wc -l` \
   "links count roundtrip identical"

# extract chr22 segments only
gfabase sub \
    "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfab" "${TMPDIR}/GRCh38-20-0.10b.chr22only.gfab" \
    chr22:1-999999999 --reference
chr22_digest=$(cat "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfa" | grep chr22 | cut -f3 | LC_ALL=C sort | sha256sum)
sub_chr22_digest=$(gfabase view "${TMPDIR}/GRCh38-20-0.10b.chr22only.gfab" | grep "^S" | cut -f3 | LC_ALL=C sort | sha256sum)
is "$sub_chr22_digest" "$chr22_digest" "chr22 reference segments"

# extract chr22 connected component
gfabase sub \
    "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfab" "${TMPDIR}/GRCh38-20-0.10b.chr22.gfab" \
    chr22:1-999999999 --reference --connected
is "$(gfabase view "${TMPDIR}/GRCh38-20-0.10b.chr22.gfab" | grep "^S" | wc -l)" "3319" "chr22-connected segments"
is "$(gfabase view "${TMPDIR}/GRCh38-20-0.10b.chr22.gfab" | grep "^L" | wc -l)" "4795" "chr22-connected links"

# sub --view to stream GFA directly
gfabase sub --view \
    "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfab" - \
    chr22:1-999999999 --reference --connected \
    > "${TMPDIR}/GRCh38-20-0.10b.chr22.gfa"
is "$(grep "^S" "${TMPDIR}/GRCh38-20-0.10b.chr22.gfa" | wc -l)" "3319" "chr22-connected segments --view"
is "$(grep "^L" "${TMPDIR}/GRCh38-20-0.10b.chr22.gfa" | wc -l)" "4795" "chr22-connected links --view"

# sub --radius-bp
gfabase sub "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfab" - \
    chr22:11,000,000-12,000,000 --reference --radius-bp 100000 --view \
    > "${TMPDIR}/radius.gfa"
is $(cat "${TMPDIR}/radius.gfa" | grep "^S" | cut -f3 | LC_ALL=C sort | sha256sum | cut -f1 -d ' ') \
   "fcbc6fd76fd6b4220aeeaad1772044ce923cc8e3d488bbcdf292c13eb668b40f" "sub --radius-bp segments"
is $(cat "${TMPDIR}/radius.gfa" | grep "^L" | wc -l) "65" "sub --radius-bp link count"

is "$(gfabase sub "${TMPDIR}/GRCh38-20-0.10b.chr22_chrY.gfab" - \
      chr22:15,000,000-15,000,001 --reference --radius-bp 999999999 --view | grep "^S" | wc -l)" \
   "3319" "sub --radius-bp 999999999"

rm -rf "$TMPDIR"
