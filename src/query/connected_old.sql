-- Compute all segments in the connected component(s) including a given set of segments.
--  IN: segment IDs in temp.start_segments
-- OUT: segment IDs in temp.connected_segments

DROP TABLE IF EXISTS temp.connected_upstream;
DROP TABLE IF EXISTS temp.connected_segments;
CREATE TABLE temp.connected_upstream(segment_id INTEGER PRIMARY KEY);
CREATE TABLE temp.connected_segments(segment_id INTEGER PRIMARY KEY);

WITH RECURSIVE
    upstream(segment_id) AS (
        SELECT segment_id FROM temp.start_segments
        UNION
        SELECT from_segment FROM {{prefix}}gfa1_link, upstream WHERE to_segment = upstream.segment_id
    )
INSERT INTO temp.connected_upstream(segment_id) SELECT segment_id FROM upstream;

WITH RECURSIVE
    connected(segment_id) AS (
        SELECT segment_id FROM temp.connected_upstream
        UNION
        SELECT to_segment FROM {{prefix}}gfa1_link, connected WHERE from_segment = connected.segment_id
    )
INSERT INTO temp.connected_segments(segment_id) SELECT segment_id FROM connected ORDER BY segment_id;

DROP TABLE temp.connected_upstream;
