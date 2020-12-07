# gfabase

### WIP for show &amp; tell, don't use yet

`gfabase` is a command-line utility for random-access storage of [Graphical Fragment Assembly (GFA)](https://github.com/GFA-spec/GFA-spec) data. It imports a .gfa file into a new **.gfab** file, from which it can later access subgraphs quickly (without reading the whole .gfab). It also has specialized features for [reference GFA (rGFA)](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md), e.g. accessing the subgraph connected to a linear coordinate range.

Effectively, .gfab is a new GFA-equivalent format with built-in compression and indexing. It is in fact a [GenomicSQLite](https://github.com/mlin/GenomicSQLite) database that `gfabase` populates with a [GFA-like schema](src/schema/GFA1.sql) to [query with SQL](src/query). Programmers have the option to access .gfab files directly using SQLite (+ Genomics Extension), without requiring `gfabase` nor even a low-level parser for .gfa/.gfab.

### Quick start

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
gfabase load <(zstd -dc test/data/GRCh38-20-0.10b.chr22_chrY.gfa.zst) /tmp/chr22_chrY.gfab --rgfa

# extract subgraph of chrY segments (only) and links between them
gfabase sub /tmp/chr22_chrY.gfab /tmp/chrYonly.gfab --reference chrY:1-999,999,999
gfabase view /tmp/chrYonly.gfab > /tmp/chrYonly.gfa

# extract the entire connected component of chr22
gfabase sub /tmp/chr22_chrY.gfab /tmp/chr22.gfab --reference --connected chr22:1-999,999,999
gfabase view /tmp/chr22.gfab > /tmp/chr22.gfa
```

The test rGFA file for chr22 &amp; chrY are originally from [the minigraph publication (Li, Feng & Chu 2020)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02168-z#availability-of-data-and-materials).
