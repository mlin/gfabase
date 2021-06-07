-- Connected components of segment graph (treating links as undirected). Omits disconnected
-- segments.
CREATE TABLE gfa1_connectivity(
    segment_id INTEGER PRIMARY KEY
        REFERENCES gfa1_segment_meta(segment_id),
    component_id INTEGER NOT NULL,  -- smallest segment_id undirectedly linked (possibly self)
    is_cutpoint INTEGER NOT NULL,    -- {0,1}, 1 iff deletion would increase # connected components
    CHECK (segment_id >= component_id)
);

-- index Walks to connected components they touch
CREATE TABLE gfa1_walk_connectivity(
    walk_id INTEGER NOT NULL
        REFERENCES gfa1_walk(walk_id),
    component_id INTEGER NOT NULL, 
        -- can't REFERENCES gfa1_connectivity(component_id) because that's not unique
    PRIMARY KEY (walk_id, component_id)
) WITHOUT ROWID;

-- Biconnected components: sets of >=3 segments that remain undirectedly connected if any one is
-- deleted. Cutpoint segments may be part of zero, one, or multiple biconnected components.
-- Non-cutpoint segments, zero or one. The ID of a biconnected component is the tuple of its min
-- and max segment ID.
CREATE TABLE gfa1_biconnectivity(
    segment_id INTEGER REFERENCES gfa1_segment_meta(segment_id),
    bicomponent_min INTEGER NOT NULL,
    bicomponent_max INTEGER NOT NULL,
    PRIMARY KEY (segment_id, bicomponent_min, bicomponent_max),
    CHECK (bicomponent_min < bicomponent_max),
    CHECK (bicomponent_min <= segment_id AND segment_id <= bicomponent_max)
) WITHOUT ROWID;

-- created in code:
-- CREATE INDEX gfa1_connectivity_component ON gfa1_connectivity(component_id);
-- CREATE INDEX gfa1_walk_connectivity_component ON gfa1_walk_connectivity(component_id);
-- CREATE INDEX gfa1_biconnectivity_component ON gfa1_biconnectivity(bicomponent_min,bicomponent_max,segment_id);
