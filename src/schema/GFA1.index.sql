CREATE UNIQUE INDEX gfa1_segment_name ON
    gfa1_segment_meta(name) WHERE name IS NOT NULL;

CREATE INDEX gfa1_segment_mapping_segment ON
    gfa1_segment_mapping(segment_id);

CREATE INDEX gfa1_link_from ON
    gfa1_link(from_segment);
CREATE INDEX gfa1_link_to ON
    gfa1_link(to_segment);

CREATE INDEX gfa1_containment_container ON
    gfa1_containment(container_segment);
CREATE INDEX gfa1_containment_contained ON
    gfa1_containment(contained_segment);

CREATE UNIQUE INDEX gfa1_path_name ON
    gfa1_path(name) WHERE name IS NOT NULL;
CREATE INDEX gfa1_path_segment ON
    gfa1_path_element(segment_id)

-- Not shown: GenomicSQLite Genomic Range Index on gfa1_segment_mapping (added in code)
