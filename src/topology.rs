use rusqlite::{params, OptionalExtension, NO_PARAMS};
use std::cmp;
use std::collections::{BTreeMap, HashSet};

use crate::util::Result;

pub fn index(db: &rusqlite::Connection) -> Result<()> {
    db.execute_batch(include_str!("schema/GFA1.topology.sql"))?;

    let mut visited_query = db.prepare("SELECT 1 from gfa1_topology WHERE segment_id = ?")?;
    let mut neighbors = db.prepare(
        // remove directionality from links
        "  SELECT from_segment FROM gfa1_link WHERE to_segment = ?1
         UNION
           SELECT to_segment FROM gfa1_link WHERE from_segment = ?1",
    )?;
    let mut insert = db.prepare(
        "INSERT INTO gfa1_topology(segment_id,component_id,cuts_component) VALUES(?,?,?)",
    )?;

    // traverse DFS forest to discover connected components
    let mut component_id: i64 = 1;
    let mut prev_component = HashSet::new();
    let mut all_segments = db.prepare("SELECT segment_id FROM gfa1_segment_meta")?;
    let mut all_segments_cursor = all_segments.query(NO_PARAMS)?;
    while let Some(segrow) = all_segments_cursor.next()? {
        let segment_id: i64 = segrow.get(0)?;
        let mut visited = prev_component.contains(&segment_id);
        if !visited  // we used prev_component like a cache to elide many visited_query ops
            && visited_query
                .query_row(params!(segment_id), |_| Ok(()))
                .optional()?
                .is_some()
        {
            visited = true;
        }
        if !visited {
            let component = component_dfs(component_id, segment_id, &mut neighbors, &mut insert)?;
            if !component.is_empty() {
                component_id += 1;
                prev_component = component;
            }
        }
    }

    db.execute_batch("CREATE INDEX gfa1_topology_component ON gfa1_topology(component_id)")?;
    Ok(())
}

// DFS traversal from given start segment; populate gfa1_topology with the discovered connected
// component, also marking the cut segments therein. https://cp-algorithms.com/graph/cutpoints.html
// Return nonempty segment ID set iff the connected component contains at least two segments.

// cutpoint algo state for each discovered segment
struct DfsSegmentState {
    // timestamp when we first reach segment
    t_in: u64,
    // lowest timestamp reachable via segment (excluding the edge with the predecessor via which we
    // first reach segment)
    t_low: u64,
    // marked true once this segment is a proven cut point
    cuts_component: bool,
}
// stack frame for iterative DFS
struct DfsStackFrame {
    // segment to visit
    segment: i64,
    // the predecessor segment; undefined when segment == start_segment_id
    pred: i64,
    // None on first reaching segment, Some(neighbor) when we return after searching a neighbor
    maybe_neighbor: Option<i64>,
}
fn component_dfs(
    component_id: i64,
    start_segment_id: i64,
    neighbors: &mut rusqlite::Statement,
    insert: &mut rusqlite::Statement,
) -> Result<HashSet<i64>> {
    let mut timestamp: u64 = 0;
    let mut state: BTreeMap<i64, DfsSegmentState> = BTreeMap::new();
    let mut start_degree: u64 = 0;

    let mut stack = vec![DfsStackFrame {
        segment: start_segment_id,
        pred: i64::MIN,
        maybe_neighbor: None,
    }];
    while let Some(frame) = stack.pop() {
        if frame.maybe_neighbor.is_none() {
            // first arrival at segment
            timestamp += 1;
            let mut segment_state = DfsSegmentState {
                t_in: timestamp,
                t_low: timestamp,
                cuts_component: false,
            };
            // enumerate segment's neighborhood (except itself & the predecessor)
            let mut neighbors_cursor = neighbors.query(params!(frame.segment))?;
            while let Some(nrow) = neighbors_cursor.next()? {
                let neighbor: i64 = nrow.get(0)?;
                if neighbor == frame.segment
                    || (frame.segment != start_segment_id && neighbor == frame.pred)
                {
                    continue;
                }
                // if we already visited this neighbor,
                if let Some(neighbor_state) = state.get(&neighbor) {
                    // record the earliest timestamp of such visits
                    segment_state.t_low = cmp::min(segment_state.t_low, neighbor_state.t_in)
                } else {
                    // schedule our return to segment after...
                    stack.push(DfsStackFrame {
                        segment: frame.segment,
                        pred: frame.pred,
                        maybe_neighbor: Some(neighbor),
                    });
                    // ...searching the neighbor next
                    stack.push(DfsStackFrame {
                        segment: neighbor,
                        pred: frame.segment,
                        maybe_neighbor: None,
                    });
                }
            }
            state.insert(frame.segment, segment_state);
        } else {
            // just finished searching one of segment's neighbors
            let neighbor = frame.maybe_neighbor.unwrap();
            // record the earliest timestamp seen amongst/beneath them
            let neighbor_low = state.get(&neighbor).unwrap().t_low;
            let ref mut segment_state = state.get_mut(&frame.segment).unwrap(); // to update in-place
            segment_state.t_low = cmp::min(segment_state.t_low, neighbor_low);
            if frame.segment != start_segment_id {
                // If nothing that'd been visited earlier is reachable via this neighbor, then
                // deleting segment would disconnect neighbor -- making it a cut segment.
                if neighbor_low >= segment_state.t_in {
                    segment_state.cuts_component = true;
                }
            } else {
                start_degree += 1;
            }
        }
    }

    // dump results into gfa1_topology
    let mut segments = HashSet::new();
    if timestamp >= 2 {
        for (segment_id, segment_state) in state.iter() {
            let cuts_component = if *segment_id != start_segment_id {
                segment_state.cuts_component
            } else {
                start_degree > 1
            };
            insert.execute(params!(
                segment_id,
                component_id,
                if cuts_component { 1 } else { 0 }
            ))?;
            segments.insert(*segment_id);
        }
    }

    Ok(segments)
}

pub fn has_index(db: &rusqlite::Connection, schema: &str) -> Result<bool> {
    Ok(db
        .query_row(
            &format!(
                "SELECT 1 FROM {}sqlite_master WHERE type='table' AND name='gfa1_topology'",
                schema
            ),
            NO_PARAMS,
            |_| Ok(()),
        )
        .optional()?
        .is_some())
}
