-- addendum to GFA1.sql defining an extra table for rGFA (equivalent of SN+SO+SR tags), enabling
-- range queries with GenomicSQLite's Genomic Range Index (GRI)
CREATE TABLE {{prefix}}gfa1_reference(
    segment_id INTEGER NOT NULL
        REFERENCES {{prefix}}gfa1_segment_meta(segment_id),
    rid INTEGER NOT NULL,       -- reference sequence ID; typically REFERENCES _gri_refsq(_gri_rid)
    position INTEGER NOT NULL,  -- zero-based offset
    length INTEGER NOT NULL,    -- segment sequence length (needed here for GRI)
    rank INTEGER NOT NULL
);
