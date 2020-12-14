-- Compute all segments in the connected component(s) including a given set of segments.
--  IN: segment IDs in temp.start_segments
-- OUT: segment IDs in temp.connected_segments
--
-- Requires SQLite 3.34.0 (2020-12-01) for the compound recursive select feature added therein. See
-- connected_old.sql for a less-elegant equivalent for older versions.

DROP TABLE IF EXISTS temp.connected_segments;
CREATE TABLE temp.connected_segments(segment_id INTEGER PRIMARY KEY);

WITH RECURSIVE
    connected(segment_id) AS (
        SELECT segment_id FROM temp.start_segments
        UNION
        SELECT from_segment FROM gfa1_link WHERE to_segment = connected.segment_id
        UNION
        SELECT to_segment FROM gfa1_link WHERE from_segment = connected.segment_id
    )
INSERT INTO temp.connected_segments(segment_id) SELECT segment_id FROM connected ORDER BY segment_id;
