-- This batch script expects a one-column temp table temp.sub_segments populated with segment_id's
-- of a desired subgraph. From attached input.{{prefix}}* it extracts the segments and the links
-- between them (only the links with both source and sink in the subgraph) into the main db.

INSERT INTO {{prefix}}gfa1_segment_sequence(segment_id, sequence_twobit)
    SELECT segment_id, sequence_twobit FROM input.{{prefix}}gfa1_segment_sequence
    WHERE segment_id IN temp.sub_segments
    ORDER BY segment_id;

INSERT INTO {{prefix}}gfa1_segment_meta(segment_id, name, tags_json)
    SELECT segment_id, name, tags_json FROM input.{{prefix}}gfa1_segment_meta
    WHERE segment_id IN temp.sub_segments
    ORDER BY segment_id;

INSERT INTO {{prefix}}gfa1_segment_mapping(segment_id, refseq_name, refseq_begin, refseq_end, cigar, tags_json)
    SELECT segment_id, refseq_name, refseq_begin, refseq_end, cigar, tags_json
    FROM input.{{prefix}}gfa1_segment_mapping
    WHERE segment_id IN temp.sub_segments
    ORDER BY segment_id;

INSERT INTO {{prefix}}gfa1_link(from_segment, from_reverse, to_segment, to_reverse, cigar, tags_json)
    SELECT from_segment, from_reverse, to_segment, to_reverse, cigar, tags_json FROM input.{{prefix}}gfa1_link
    WHERE from_segment IN temp.sub_segments AND to_segment IN temp.sub_segments
    ORDER BY from_segment, to_segment;
