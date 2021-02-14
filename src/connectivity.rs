// Connectivity index: at the end of the load process, traverse the DFS forest to discover
// connected components (treating the segment graph as undirected). Store a relation table
// annotating each segment with its connected component, identified by the smallest connected
// segment_id, and whether it's a "cutpoint" whose individual deletion would increase the number
// of connected components. Disconnected segments are omitted from the table.

use bloomfilter::Bloom;
use rusqlite::{params, OptionalExtension, NO_PARAMS};
use std::cmp;
use std::collections::BTreeMap;

use crate::util::Result;

pub fn index(db: &rusqlite::Connection) -> Result<()> {
    db.execute_batch(include_str!("schema/GFA1.connectivity.sql"))?;

    let mut neighbors = db.prepare(
        // remove directionality from links
        "  SELECT from_segment FROM gfa1_link WHERE to_segment = ?1 AND from_segment != ?1
         UNION
           SELECT to_segment FROM gfa1_link WHERE from_segment = ?1 AND to_segment != ?1",
    )?;
    let mut insert = db.prepare(
        "INSERT INTO gfa1_connectivity(segment_id,component_id,is_cutpoint) VALUES(?,?,?)",
    )?;

    let mut visited_query = db.prepare("SELECT 1 from gfa1_connectivity WHERE segment_id = ?")?;
    // use a bloom filter in front of visited_query
    let approx_segment_count: i64 = db.query_row(
        "SELECT coalesce(max(segment_id),100000) FROM gfa1_segment_meta",
        NO_PARAMS,
        |row| row.get(0),
    )?;
    let mut visited_bloom = Bloom::new_for_fp_rate(approx_segment_count as usize, 0.05);

    // traverse DFS forest to discover connected components
    let mut all_segments = db.prepare("SELECT segment_id FROM gfa1_segment_meta")?;
    let mut all_segments_cursor = all_segments.query(NO_PARAMS)?;
    while let Some(segrow) = all_segments_cursor.next()? {
        let segment_id: i64 = segrow.get(0)?;
        if !(visited_bloom.check(&segment_id)
            && visited_query
                .query_row(params!(segment_id), |_| Ok(()))
                .optional()?
                .is_some())
        {
            component_dfs(segment_id, &mut neighbors, &mut insert, &mut visited_bloom)?
        }
    }

    db.execute_batch(
        "CREATE INDEX gfa1_connectivity_component ON gfa1_connectivity(component_id)",
    )?;
    Ok(())
}

// DFS traversal from given start segment; populate gfa1_connectivity with the discovered connected
// component, also marking is_cutpoint therein. https://cp-algorithms.com/graph/cutpoints.html

// cutpoint algo state for each discovered segment
struct DfsSegmentState {
    // timestamp when we first reach segment
    t_in: u64,
    // lowest timestamp reachable via segment (excluding the edge with the predecessor via which we
    // first reach segment)
    t_low: u64,
    // marked true once this segment is a proven cut point
    is_cutpoint: bool,
}
// stack frames for iterative DFS
enum DfsStackFrame {
    Arrive { segment: i64, parent: i64 },
    Return { segment: i64, child: i64 },
}
fn component_dfs(
    start_segment_id: i64,
    neighbors: &mut rusqlite::Statement,
    insert: &mut rusqlite::Statement,
    visited_bloom: &mut Bloom<i64>,
) -> Result<()> {
    let mut timestamp: u64 = 0;
    let mut state: BTreeMap<i64, DfsSegmentState> = BTreeMap::new();
    let mut start_degree: u64 = 0;

    let mut stack = vec![DfsStackFrame::Arrive {
        segment: start_segment_id,
        parent: i64::MIN, // undefined for start segment
    }];
    while let Some(frame) = stack.pop() {
        match frame {
            DfsStackFrame::Arrive { segment, parent } => {
                if let Some(t_in) = state.get(&segment).map(|segment_state| segment_state.t_in) {
                    // already visited this segment; reduce parent t_low to the first such visit
                    assert!(timestamp > 1);
                    let ref mut pt_state = state.get_mut(&parent).unwrap();
                    pt_state.t_low = cmp::min(pt_state.t_low, t_in)
                } else {
                    // first arrival at segment
                    timestamp += 1;
                    state.insert(
                        segment,
                        DfsSegmentState {
                            t_in: timestamp,
                            t_low: timestamp,
                            is_cutpoint: false,
                        },
                    );
                    // schedule return to parent after...
                    if segment != start_segment_id {
                        stack.push(DfsStackFrame::Return {
                            segment: parent,
                            child: segment,
                        });
                    }
                    // visiting segment's other neighbors
                    let mut neighbors_cursor = neighbors.query(params!(segment))?;
                    while let Some(nrow) = neighbors_cursor.next()? {
                        let neighbor: i64 = nrow.get(0)?;
                        if segment != start_segment_id && neighbor == parent {
                            continue;
                        }
                        stack.push(DfsStackFrame::Arrive {
                            segment: neighbor,
                            parent: segment,
                        })
                    }
                }
            }
            DfsStackFrame::Return { segment, child } => {
                // returning to segment after completing (what turned out to be) the first visit
                // to child; reduce segment's t_low to child's
                let child_low = state.get(&child).unwrap().t_low;
                let ref mut segment_state = state.get_mut(&segment).unwrap();
                segment_state.t_low = cmp::min(segment_state.t_low, child_low);
                if segment != start_segment_id {
                    // If none of segment's ancestors were reachable via child, then deleting
                    // segment would disconnect child.
                    if child_low >= segment_state.t_in {
                        segment_state.is_cutpoint = true;
                    }
                } else {
                    start_degree += 1;
                }
            }
        }
    }

    if timestamp < 2 {
        return Ok(());
    }
    // dump results into gfa1_connectivity
    let mut component_id = None;
    for (segment_id, segment_state) in state.iter() {
        if component_id.is_none() {
            component_id = Some(segment_id) // smallest segment_id
        }
        let is_cutpoint = if *segment_id != start_segment_id {
            segment_state.is_cutpoint
        } else {
            start_degree > 1
        };
        insert.execute(params!(
            segment_id,
            component_id.unwrap(),
            if is_cutpoint { 1 } else { 0 },
        ))?;
        visited_bloom.set(segment_id);
    }
    Ok(())
}

pub fn has_index(db: &rusqlite::Connection, schema: &str) -> Result<bool> {
    Ok(db
        .query_row(
            &format!(
                "SELECT 1 FROM {}sqlite_master WHERE type='table' AND name='gfa1_connectivity'",
                schema
            ),
            NO_PARAMS,
            |_| Ok(()),
        )
        .optional()?
        .is_some())
}
