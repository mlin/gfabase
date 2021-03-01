-- Connected components of segment graph (treating links as undirected). Omits disconnected
-- segments.
CREATE TABLE gfa1_connectivity(
    segment_id INTEGER PRIMARY KEY
        REFERENCES gfa1_segment_meta(segment_id),
    component_id INTEGER NOT NULL,  -- smallest segment_id undirectedly linked (possibly self)
    is_cutpoint INTEGER NOT NULL    -- {0,1}, 1 iff deletion would increase # connected components
);

-- index Walks to connected components they touch
CREATE TABLE gfa1_walk_connectivity(
    walk_id INTEGER NOT NULL REFERENCES gfa1_walk(walk_id),
    component_id INTEGER NOT NULL REFERENCES gfa1_connectivity(component_id),
    PRIMARY KEY (walk_id, component_id)
) WITHOUT ROWID;

-- created in code:
-- CREATE INDEX gfa1_connectivity_component ON gfa1_connectivity(component_id);
-- CREATE INDEX gfa1_walk_connectivity_component ON gfa1_walk_connectivity(component_id);
