-- Compute all segments in the connected component(s) including a given set of segments.
--  IN: segment IDs in temp.start_segments
-- OUT: segment IDs in temp.connected_segments

DROP TABLE IF EXISTS temp.connected_segments;
CREATE TABLE temp.connected_segments(segment_id INTEGER PRIMARY KEY);

WITH RECURSIVE
    connected(segment_id) AS (
        SELECT segment_id FROM temp.start_segments
        UNION
        SELECT
            IIF(to_segment = connected.segment_id, from_segment, to_segment)
            FROM gfa1_link, connected
            WHERE to_segment = connected.segment_id OR from_segment = connected.segment_id
    )
    INSERT INTO temp.connected_segments(segment_id)
        SELECT segment_id FROM connected ORDER BY segment_id;
