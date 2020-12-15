# gfabase

### WIP for show &amp; tell, don't use yet

`gfabase` is a command-line tool for random-access storage of [Graphical Fragment Assembly (GFA)](https://github.com/GFA-spec/GFA-spec) data. It imports a .gfa file into a compressed **.gfab** file, from which it can later access subgraphs quickly (reading only the necessary parts), producing .gfa or .gfab. Beyond ID-based access, .gfab can index segments by their mappings to linear reference coordinates. This facilitates navigation within *de novo* assemblies or [rGFA reference graphs](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md), with quick access to subgraphs connected to linear coordinate ranges.

Effectively, .gfab is a new GFA-superset format with built-in compression and indexing. It is in fact a SQLite (+ [Genomics Extension](https://github.com/mlin/GenomicSQLite)) database populated with a [GFA-like schema](src/schema/GFA1.sql). Programmers have the option to access .gfab files directly using SQLite (+ Genomics Extension), without requiring `gfabase` nor even a low-level parser for .gfa/.gfab.

### Quick start

![CI](https://github.com/mlin/gfabase/workflows/CI/badge.svg?branch=main)

[Install rust toolchain](https://rustup.rs/) and:

```bash
git clone https://github.com/mlin/gfabase.git
cd gfabase
cargo build --release
```

Then, a few things to try:

```bash
alias gfabase=$(pwd)/target/release/gfabase

# import a rGFA file into .gfab
gfabase load <(zstd -dc test/data/GRCh38-20-0.10b.chr22_chrY.gfa.zst) /tmp/chr22_chrY.gfab

# extract subgraph of chrY segments (only) and links between them
gfabase sub /tmp/chr22_chrY.gfab /tmp/chrYonly.gfab --reference chrY:1-999,999,999
gfabase view /tmp/chrYonly.gfab > /tmp/chrYonly.gfa

# extract the entire connected component of chr22
gfabase sub /tmp/chr22_chrY.gfab /tmp/chr22.gfa --view --reference --connected chr22:1-999,999,999
```

The test rGFA file for chr22 &amp; chrY are originally from [the minigraph publication (Li, Feng & Chu 2020)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02168-z#availability-of-data-and-materials).
