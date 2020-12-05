# gfabase

## WIP for show &amp; tell; do not use yet

`gfabase` is a command-line utility for random-access storage of [Graphical Fragment Assembly (GFA)](https://github.com/GFA-spec/GFA-spec) data. It imports a .gfa file into a new **.gfab** file, which it can later search and operate on subgraphs without first loading everything into memory. It also has specialized features for [reference GFA (rGFA)](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md), e.g. accessing the subgraph connected to a linear coordinate range.

Effectively, .gfab is a new GFA-equivalent format with built-in compression and indexing. It is in fact a [GenomicSQLite](https://github.com/mlin/GenomicSQLite) database that `gfabase` populates with a GFA-like schema to query with SQL. Programmers have the option to access .gfab files directly using SQLite (+ Genomics Extension), without requiring `gfabase` nor even a low-level parser for .gfa/.gfab.
