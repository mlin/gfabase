-- Segment metadata
CREATE TABLE gfa1_segment_meta(
    segment_id INTEGER NOT NULL PRIMARY KEY,
    name TEXT COLLATE UINT,   -- if distinct from ID, otherwise NULL
    sequence_length INTEGER,  -- if no sequence (*): length taken from LN:i tag if present, otherwise NULL
    tags_json TEXT
);

-- Segment sequences (stored separately & with two-bit encoding)
CREATE TABLE gfa1_segment_sequence(
    segment_id INTEGER NOT NULL PRIMARY KEY
        REFERENCES gfa1_segment_meta(segment_id),
    sequence_twobit BLOB NOT NULL          -- not null: omit row if no sequence available
);

-- Convenience: view joining segment_meta & segment_sequence
CREATE VIEW gfa1_segment AS
    SELECT
        segment_id, name, sequence_length, tags_json,
        twobit_dna(sequence_twobit) AS sequence
    FROM
        gfa1_segment_meta LEFT JOIN gfa1_segment_sequence
        USING (segment_id);

-- gfabase extension: "mappings" of segments to linear reference sequence coordinates, suitable for
-- Genomic Range Indexing
CREATE TABLE gfa1_segment_mapping(
    segment_id INTEGER NOT NULL              -- nb: one segment may have multiple mappings
        REFERENCES gfa1_segment_meta(segment_id),
    refseq_name TEXT NOT NULL COLLATE UINT,  -- associated reference sequence (e.g. chromosome name)
    refseq_begin INTEGER NOT NULL,           -- zero-based begin of associated range
    refseq_end INTEGER NOT NULL,             -- end (exclusive) of associated range
    cigar TEXT,                              -- alignment of segment sequence to associated refseq range, if known
    tags_json TEXT                           -- currently unused
);

-- Link
CREATE TABLE gfa1_link(
    from_segment INTEGER NOT NULL
        REFERENCES gfa1_segment_meta(segment_id),
    from_reverse INTEGER NOT NULL,                  -- {0,1}
    to_segment INTEGER NOT NULL
        REFERENCES gfa1_segment_meta(segment_id),
    to_reverse INTEGER NOT NULL,                    -- {0,1}
    cigar TEXT,
    tags_json TEXT
);

-- Containment
CREATE TABLE gfa1_containment(
    container_segment INTEGER NOT NULL
        REFERENCES gfa1_segment_meta(segment_id),
    container_reverse INTEGER NOT NULL,             -- {0,1}
    contained_segment INTEGER NOT NULL
        REFERENCES gfa1_segment_meta(segment_id),
    contained_reverse INTEGER NOT NULL,             -- {0,1}
    position INTEGER NOT NULL,
    cigar TEXT,
    tags_json TEXT
);

-- Path
CREATE TABLE gfa1_path(
    path_id INTEGER NOT NULL PRIMARY KEY,
    name TEXT COLLATE UINT,  -- if distinct from ID, otherwise NULL
    tags_json TEXT
);

-- Path elements (derived from Path SegmentNames and Overlaps)
CREATE TABLE gfa1_path_element(
    path_id INTEGER NOT NULL
        REFERENCES gfa1_path(path_id),
    ordinal INTEGER NOT NULL,
    segment_id INTEGER NOT NULL
        REFERENCES gfa1_segment_meta(segment_id),
    reverse INTEGER NOT NULL,                       -- {0,1}
    cigar_vs_previous TEXT,                         -- NULL for the first element
    PRIMARY KEY (path_id,ordinal)
) WITHOUT ROWID;

-- Header
CREATE TABLE gfa1_header(
    tags_json TEXT
);
