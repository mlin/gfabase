// Connectivity index: at the end of the load process, traverse the DFS forest to discover
// connected components (treating the segment graph as undirected). Store a relation table
// annotating each segment with its connected component, identified by the smallest connected
// segment_id, and whether it's a "cutpoint" whose individual deletion would increase the number
// of connected components. Disconnected segments are omitted from the table.
//
// In the same pass, also discover biconnected components, sets of >=3 segments which remain
// connected following deletion of any one. Store a relation table annotating which biconnected
// component(s) each segment is part of (possibly multiple for cutpoint segments). The ID of a
// biconnected component is the tuple of its min and max constituent segment IDs.

use bloomfilter::Bloom;
use rusqlite::{params, OptionalExtension};
use std::cmp;
use std::collections::{BTreeMap, BTreeSet};

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
    let mut insert_bicon = db.prepare(
        "INSERT INTO gfa1_biconnectivity(segment_id,bicomponent_min,bicomponent_max) VALUES(?,?,?)",
    )?;

    let mut visited_query = db.prepare("SELECT 1 from gfa1_connectivity WHERE segment_id = ?")?;
    // use a bloom filter in front of visited_query
    let approx_segment_count: i64 = db.query_row(
        "SELECT coalesce(max(segment_id),100000) FROM gfa1_segment_meta",
        [],
        |row| row.get(0),
    )?;
    let mut visited_bloom = Bloom::new_for_fp_rate(approx_segment_count as usize, 0.05);

    // traverse DFS forest to discover connected components
    let mut all_segments = db.prepare("SELECT segment_id FROM gfa1_segment_meta")?;
    let mut all_segments_cursor = all_segments.query([])?;
    while let Some(segrow) = all_segments_cursor.next()? {
        let segment_id: i64 = segrow.get(0)?;
        if !(visited_bloom.check(&segment_id)
            && visited_query
                .query_row(params!(segment_id), |_| Ok(()))
                .optional()?
                .is_some())
        {
            component_dfs(
                segment_id,
                &mut neighbors,
                &mut insert,
                &mut insert_bicon,
                &mut visited_bloom,
            )?
        }
    }

    // index each Walk to the associated connected component. By definition, all segments in a Walk
    // must be in one connected component, so it suffices just to look up one exemplar segment.
    // Also, checking all of them would be costly.
    db.execute_batch(
        "INSERT INTO gfa1_walk_connectivity(walk_id,component_id)
         SELECT walk_id, component_id
         FROM gfa1_walk INNER JOIN gfa1_connectivity ON gfa1_walk.min_segment_id = gfa1_connectivity.segment_id"
    )?;

    db.execute_batch(
        "CREATE INDEX gfa1_connectivity_component ON gfa1_connectivity(component_id);
         CREATE INDEX gfa1_walk_connectivity_component ON gfa1_walk_connectivity(component_id);
         CREATE INDEX gfa1_biconnectivity_component ON gfa1_biconnectivity(bicomponent_min,bicomponent_max,segment_id)",
    )?;
    Ok(())
}

// DFS traversal from given start segment; populate gfa1_connectivity with the discovered connected
// component, also marking its cutpoints and biconnected components. refs:
//     https://cp-algorithms.com/graph/cutpoints.html
//     https://www.cs.cmu.edu/~avrim/451f12/lectures/biconnected.pdf

// cutpoint algo state for each discovered segment
struct DfsSegmentState {
    // timestamp when we first reach segment
    t_in: u64,
    // lowest timestamp reachable via segment (excluding the edge with the predecessor via which we
    // first reach segment)
    t_low: u64,
    // marked true once this segment is a proven cut point
    is_cutpoint: bool,
    // biconnected component(s) in which this segment resides (possibly zero or multiple for
    // cutpoints). the id of a biconnected component is the tuple of its min and max segment IDs.
    bicon_components: BTreeSet<(i64, i64)>,
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
    insert_bicon: &mut rusqlite::Statement,
    visited_bloom: &mut Bloom<i64>,
) -> Result<()> {
    let mut timestamp: u64 = 0;
    let mut state: BTreeMap<i64, DfsSegmentState> = BTreeMap::new();
    let mut start_returns: u64 = 0;

    let mut stack = vec![DfsStackFrame::Arrive {
        segment: start_segment_id,
        parent: i64::MIN, // undefined for start segment
    }];
    let mut bicon_stack: Vec<(i64, i64)> = vec![(i64::MIN, start_segment_id)];
    while let Some(frame) = stack.pop() {
        match frame {
            DfsStackFrame::Arrive { segment, parent } => {
                assert_ne!(segment, i64::MIN);
                if let Some(t_in) = state.get(&segment).map(|segment_state| segment_state.t_in) {
                    // previously visited segment
                    assert!(timestamp > 1 && parent > i64::MIN);
                    let ref mut pt_state = state.get_mut(&parent).unwrap();
                    if t_in < pt_state.t_in {
                        // cycle back to ancestor of parent; update parent t_low
                        pt_state.t_low = cmp::min(pt_state.t_low, t_in);
                        bicon_stack.push((parent, segment))
                    }
                } else {
                    // first visit to segment
                    timestamp += 1;
                    state.insert(
                        segment,
                        DfsSegmentState {
                            t_in: timestamp,
                            t_low: timestamp,
                            is_cutpoint: false,
                            bicon_components: BTreeSet::new(),
                        },
                    );
                    bicon_stack.push((parent, segment));
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
                    // segment would disconnect child, therefore segment is a cutpoint.
                    if child_low >= segment_state.t_in {
                        segment_state.is_cutpoint = true;
                        pop_bicon_component((segment, child), &mut bicon_stack, &mut state);
                    }
                } else {
                    start_returns += 1;
                }
            }
        }
    }

    if timestamp < 2 {
        return Ok(());
    }

    // postprocess the start segment
    if start_returns > 1 {
        state.get_mut(&start_segment_id).unwrap().is_cutpoint = true;
    }
    pop_bicon_component((i64::MIN, start_segment_id), &mut bicon_stack, &mut state);

    // dump results into gfa1_connectivity
    let mut component_id = i64::MIN;
    for (segment_id, segment_state) in state {
        if component_id == i64::MIN {
            component_id = segment_id // smallest segment_id
        }
        insert.execute(params!(segment_id, component_id, segment_state.is_cutpoint))?;
        visited_bloom.set(&segment_id);
        for (bicomponent_min, bicomponent_max) in segment_state.bicon_components {
            insert_bicon.execute(params!(segment_id, bicomponent_min, bicomponent_max))?;
        }
    }
    Ok(())
}

// subroutine of component_dfs: after discovering a cutpoint (or returning to the start segment),
// recover the just-traversed biconnected component, if any, and mark it in state
fn pop_bicon_component(
    entry: (i64, i64),
    bicon_stack: &mut Vec<(i64, i64)>,
    state: &mut BTreeMap<i64, DfsSegmentState>,
) {
    // pop from edge stack to get a segment set
    let mut segments = BTreeSet::new();
    loop {
        let (s0, s1) = bicon_stack.pop().unwrap();
        if s0 > i64::MIN {
            segments.insert(s0);
        }
        segments.insert(s1);
        if (s0, s1) == entry {
            break;
        }
    }
    let bcc_size = segments.len();
    // we found a biconnected component iff bcc_size > 2
    if bcc_size > 2 {
        let mut bcc_min: i64 = i64::MAX;
        let mut bcc_max: i64 = i64::MIN;
        for elt in segments.iter() {
            bcc_min = cmp::min(bcc_min, *elt);
            bcc_max = cmp::max(bcc_max, *elt);
        }
        assert!(i64::MIN < bcc_min && bcc_min < bcc_max && bcc_max < i64::MAX);
        assert_eq!(segments.len(), bcc_size);
        // update state
        for elt in segments {
            state
                .get_mut(&elt)
                .unwrap()
                .bicon_components
                .insert((bcc_min, bcc_max));
        }
    }
}

pub fn has_index(db: &rusqlite::Connection, schema: &str) -> Result<bool> {
    Ok(db
        .query_row(
            &format!(
                "SELECT 1 FROM {}sqlite_master WHERE type='table' AND name='gfa1_connectivity'",
                schema
            ),
            [],
            |_| Ok(()),
        )
        .optional()?
        .is_some())
}
