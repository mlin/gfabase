CREATE UNIQUE INDEX {{prefix}}gfa1_segment_name ON
    {{prefix}}gfa1_segment_meta(name) WHERE name IS NOT NULL;

CREATE INDEX {{prefix}}gfa1_link_from_to ON
    {{prefix}}gfa1_link(from_segment,to_segment);
CREATE INDEX {{prefix}}gfa1_link_to_from ON
    {{prefix}}gfa1_link(to_segment,from_segment);

CREATE INDEX {{prefix}}gfa1_containment_container_contained ON
    {{prefix}}gfa1_containment(container_segment,contained_segment);
CREATE INDEX {{prefix}}gfa1_containment_contained_container ON
    {{prefix}}gfa1_containment(contained_segment,container_segment);

CREATE UNIQUE INDEX {{prefix}}gfa1_path_name ON
    {{prefix}}gfa1_path(name) WHERE name IS NOT NULL;
CREATE INDEX {{prefix}}gfa1_path_segment ON
    {{prefix}}gfa1_path_element(segment_id);
