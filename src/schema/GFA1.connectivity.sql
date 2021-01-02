-- Connected components of segment graph (treating links as undirected). Lone disconnected segments
-- are omitted from this table.
CREATE TABLE gfa1_connectivity(
    segment_id INTEGER PRIMARY KEY
        REFERENCES gfa1_segment_meta(segment_id),
    component_id INTEGER NOT NULL,
    is_cutpoint INTEGER NOT NULL  -- {0,1}: 1 iff deleting this segment (& links touching it) would
                                  -- cut its connected component into two or more
);

-- created in code:
-- CREATE INDEX gfa1_connectivity_component ON gfa1_connectivity(component_id);
