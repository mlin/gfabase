-- This batch script expects a one-column temp table temp.sub_segments populated with segment_id's
-- of a desired subgraph. From attached input.* it copies the segments and the links
-- between them (only the links with both source and sink in the subgraph) into the main db.

INSERT INTO gfa1_segment_meta(segment_id, name, sequence_length, tags_json)
    SELECT segment_id, name, sequence_length, tags_json FROM input.gfa1_segment_meta
    WHERE segment_id IN temp.sub_segments;

-- gfa1_segment_sequences copied in code (if not --no-sequences)

INSERT INTO gfa1_segment_mapping(segment_id, refseq_name, refseq_begin, refseq_end, tags_json)
    SELECT segment_id, refseq_name, refseq_begin, refseq_end, tags_json
    FROM input.gfa1_segment_mapping
    WHERE segment_id IN temp.sub_segments
    ORDER BY segment_id;

INSERT INTO gfa1_link(from_segment, from_reverse, to_segment, to_reverse, tags_json)
    SELECT from_segment, from_reverse, to_segment, to_reverse, tags_json FROM input.gfa1_link
    WHERE +from_segment IN temp.sub_segments AND to_segment IN temp.sub_segments
    -- FIXME: the unary plus hint +from_segment is a temporary workaround for a SQLite problem:
    --        https://sqlite.org/forum/forumpost/b4fcb8a598?t=h
    ORDER BY from_segment, to_segment;

-- Identify and copy the paths with no segments missing from temp.sub_segments
CREATE TABLE temp.sub_paths(path_id INTEGER PRIMARY KEY);

INSERT INTO temp.sub_paths(path_id)
    SELECT path_id FROM input.gfa1_path
    WHERE path_id NOT IN
        -- paths with missing segments:
        (SELECT DISTINCT path_id FROM input.gfa1_path_element
         WHERE segment_id NOT IN temp.sub_segments);

INSERT INTO gfa1_path(path_id, name, tags_json)
    SELECT path_id, name, tags_json FROM input.gfa1_path
    WHERE path_id IN temp.sub_paths;

INSERT INTO gfa1_path_element(path_id, ordinal, segment_id, reverse, cigar_vs_previous)
    SELECT path_id, ordinal, segment_id, reverse, cigar_vs_previous FROM input.gfa1_path_element
    WHERE path_id IN temp.sub_paths;

--DROP TABLE temp.sub_paths;

INSERT INTO gfa1_header(_rowid_, tags_json)
    SELECT _rowid_, tags_json FROM input.gfa1_header
    WHERE _rowid_ = 1;
