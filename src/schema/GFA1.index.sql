CREATE UNIQUE INDEX gfa1_segment_name ON
    gfa1_segment_meta(name) WHERE name IS NOT NULL;
CREATE INDEX gfa1_connected_component ON
    gfa1_segment_meta(connected_component);

CREATE INDEX gfa1_segment_mapping_segment ON
    gfa1_segment_mapping(segment_id);

CREATE INDEX gfa1_containment_container_contained ON
    gfa1_containment(container_segment,contained_segment);
CREATE INDEX gfa1_containment_contained_container ON
    gfa1_containment(contained_segment,container_segment);

CREATE UNIQUE INDEX gfa1_path_name ON
    gfa1_path(name) WHERE name IS NOT NULL;
CREATE INDEX gfa1_path_segment ON
    gfa1_path_element(segment_id)

-- Not shown:
-- 1. Not shown: gfa1_link indexes (created in code because they're needed for connectivity analysis)
-- 2. GenomicSQLite Genomic Range Index on gfa1_segment_mapping (added in code)
