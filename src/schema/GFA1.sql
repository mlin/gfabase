-- Segment
CREATE TABLE {{prefix}}gfa1_segment(
    _rowid_ INTEGER NOT NULL PRIMARY KEY,
    name TEXT,                             -- if distinct from _rowid_, otherwise NULL
    tags_json TEXT NOT NULL DEFAULT '{}',
    sequence_twobit BLOB                   -- nucleotides_twobit() of the sequence text
);

-- Link
CREATE TABLE {{prefix}}gfa1_link(
    from_segment INTEGER NOT NULL
        REFERENCES {{prefix}}gfa1_segment(_rowid_),
    from_reverse INTEGER NOT NULL,                  -- {0,1}
    to_segment INTEGER NOT NULL
        REFERENCES {{prefix}}gfa1_segment(_rowid_),
    to_reverse INTEGER NOT NULL,                    -- {0,1}
    cigar TEXT,
    tags_json TEXT NOT NULL DEFAULT '{}'
);

-- Containment
CREATE TABLE {{prefix}}gfa1_containment(
    container_segment INTEGER NOT NULL
        REFERENCES {{prefix}}gfa1_segment(_rowid_),
    container_reverse INTEGER NOT NULL,             -- {0,1}
    contained_segment INTEGER NOT NULL
        REFERENCES {{prefix}}gfa1_segment(_rowid_),
    contained_reverse INTEGER NOT NULL,             -- {0,1}
    position INTEGER NOT NULL,
    cigar TEXT,
    tags_json TEXT NOT NULL DEFAULT '{}'
);

-- Path
CREATE TABLE {{prefix}}gfa1_path(
    _rowid_ INTEGER NOT NULL PRIMARY KEY,
    name TEXT UNIQUE,                      -- if distinct from _rowid_, otherwise NULL
    tags_json TEXT NOT NULL DEFAULT '{}'
);

-- Path elements (derived from Path SegmentNames and Overlaps)
CREATE TABLE {{prefix}}gfa1_path_element(
    path INTEGER NOT NULL
        REFERENCES {{prefix}}gfa1_path(_rowid_),
    ordinal INTEGER NOT NULL,
    segment INTEGER NOT NULL
        REFERENCES {{prefix}}gfa1_segment(_rowid_),
    reverse INTEGER NOT NULL,                       -- {0,1}
    cigar_vs_prev TEXT,                             -- NULL for the first element
    PRIMARY KEY (path,ordinal)
) WITHOUT ROWID;
