use std::{env, time::Instant};

use rectangles::parse_input;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Point {
    y: i64,
    x: i64,
    inner_diags: [bool; 4], // [UpLeft, UpRight, DownLeft, DownRight]
}

#[cfg(test)]
impl Point {
    fn new(y: i64, x: i64) -> Self {
        Point {
            y,
            x,
            inner_diags: [false; 4],
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Dir {
    Up,
    Down,
    Left,
    Right,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DiagDir {
    UpLeft = 0,
    UpRight = 1,
    DownLeft = 2,
    DownRight = 3,
}

#[inline]
fn dir(coord: (i64, i64), other: (i64, i64)) -> Dir {
    if other.0 < coord.0 {
        Dir::Up
    } else if other.0 > coord.0 {
        Dir::Down
    } else if other.1 < coord.1 {
        Dir::Left
    } else {
        Dir::Right
    }
}

/// Returns +1 for clockwise (right turn), -1 for counter-clockwise (left turn), 0 for straight
#[inline]
fn turn_direction(from: Dir, to: Dir) -> i32 {
    use Dir::*;

    match (from, to) {
        (Up, Right) | (Right, Down) | (Down, Left) | (Left, Up) => 1, // Clockwise
        (Up, Left) | (Left, Down) | (Down, Right) | (Right, Up) => -1, // Counter-clockwise
        (Up, Down) | (Down, Up) | (Left, Right) | (Right, Left) => {
            panic!("U-turn - assuming this isn't possible.")
        }
        _ => 0, // Straight
    }
}

/// Returns the inner side direction when walking in `dir`
/// For clockwise polygons, inner is on the right; for counter-clockwise, inner is on the left
fn inner_side(d: Dir, clockwise: bool) -> Dir {
    use Dir::*;
    if clockwise {
        // Right side
        match d {
            Up => Right,
            Right => Down,
            Down => Left,
            Left => Up,
        }
    } else {
        // Left side
        match d {
            Up => Left,
            Left => Down,
            Down => Right,
            Right => Up,
        }
    }
}

/// Returns the diagonal direction of `other` relative to `point`
fn diag_dir_of(point: Point, other: Point) -> Option<DiagDir> {
    use DiagDir::*;
    if other.y < point.y && other.x < point.x {
        Some(UpLeft)
    } else if other.y < point.y && other.x > point.x {
        Some(UpRight)
    } else if other.y > point.y && other.x < point.x {
        Some(DownLeft)
    } else if other.y > point.y && other.x > point.x {
        Some(DownRight)
    } else {
        None // On same row or column, no diagonal
    }
}

/// Check if `r1`'s inner directions point towards `r2`.
/// For diagonal directions, the exact diagonal must be inner.
/// For cardinal directions (same row/column), either adjacent diagonal being inner is sufficient.
fn inner_dirs_point_towards(r1: Point, r2: Point) -> bool {
    use DiagDir::*;
    let inner = r1.inner_diags;

    if let Some(diag) = diag_dir_of(r1, r2) {
        // Diagonal: must be in inner_diags
        inner[diag as usize]
    } else if r2.y < r1.y {
        // Directly Up: valid if UpLeft or UpRight is inner
        inner[UpLeft as usize] || inner[UpRight as usize]
    } else if r2.y > r1.y {
        // Directly Down: valid if DownLeft or DownRight is inner
        inner[DownLeft as usize] || inner[DownRight as usize]
    } else if r2.x < r1.x {
        // Directly Left: valid if UpLeft or DownLeft is inner
        inner[UpLeft as usize] || inner[DownLeft as usize]
    } else if r2.x > r1.x {
        // Directly Right: valid if UpRight or DownRight is inner
        inner[UpRight as usize] || inner[DownRight as usize]
    } else {
        // Same point
        false
    }
}

/// Returns the diagonal direction from two cardinal directions
fn to_diag(d1: Dir, d2: Dir) -> DiagDir {
    use DiagDir::*;
    use Dir::*;
    match (d1, d2) {
        (Up, Left) | (Left, Up) => UpLeft,
        (Up, Right) | (Right, Up) => UpRight,
        (Down, Left) | (Left, Down) => DownLeft,
        (Down, Right) | (Right, Down) => DownRight,
        _ => panic!("Invalid direction combo for diagonal"),
    }
}

/// Returns the opposite diagonal
fn opposite_diag(d: DiagDir) -> DiagDir {
    use DiagDir::*;
    match d {
        UpLeft => DownRight,
        UpRight => DownLeft,
        DownLeft => UpRight,
        DownRight => UpLeft,
    }
}

/// Get neighbors at index i with wraparound
fn neighbors(coords: &[(i64, i64)], i: usize) -> ((i64, i64), (i64, i64)) {
    let len = coords.len();
    (coords[(i + len - 1) % len], coords[(i + 1) % len])
}

fn preprocess_to_points(coords: &[(i64, i64)]) -> Vec<Point> {
    // Compute turn info about our shape:
    // * the direction of each incoming and outgoing edge
    // * the "turn" type - counterclockwise, clockwise, or straight
    let turn_info: Vec<_> = coords
        .iter()
        .enumerate()
        .map(|(i, &coord)| {
            let (before, after) = neighbors(coords, i);
            let from_dir = dir(before, coord);
            let to_dir = dir(coord, after);
            (from_dir, to_dir, turn_direction(from_dir, to_dir))
        })
        .collect();

    let total_turn: i32 = turn_info.iter().map(|(_, _, t)| t).sum();
    let is_clockwise = total_turn > 0;

    // Now calculate which diagonal directions from each coord point into
    // space *within* our shape: these are valid directions we can form a
    // rectangle into.
    coords
        .iter()
        .zip(turn_info)
        .map(|(&(y, x), (from_dir, to_dir, turn))| {
            // Get the inner side of each edge
            let incoming_inner = inner_side(from_dir, is_clockwise);
            let outgoing_inner = inner_side(to_dir, is_clockwise);

            // The diagonal in the "inner quadrant" formed by both edges
            let inner_diag = to_diag(incoming_inner, outgoing_inner);

            // Convex corner: only the inner diagonal is inside
            // Concave corner: all diagonals except the outer one are inside
            let is_concave = if is_clockwise { turn < 0 } else { turn > 0 };

            let mut inner_diags = [is_concave; 4];
            let toggle_diag = if is_concave {
                opposite_diag(inner_diag)
            } else {
                inner_diag
            };
            inner_diags[toggle_diag as usize] = !is_concave;

            Point { y, x, inner_diags }
        })
        .collect()
}

/// Calculate the area of a rectangle defined by opposite corners r1 and r2.
#[inline]
fn area(r1: Point, r2: Point) -> u64 {
    (((r2.y - r1.y).abs() + 1) * ((r2.x - r1.x).abs() + 1)) as u64
}

/// Check if any point from the x-sorted list is strictly within the rectangle.
/// Uses binary search to find candidate points efficiently.
fn any_point_within(r1: Point, r2: Point, points_by_x: &[(i64, i64)]) -> bool {
    let min_y = r1.y.min(r2.y);
    let max_y = r1.y.max(r2.y);
    let min_x = r1.x.min(r2.x);
    let max_x = r1.x.max(r2.x);

    // Binary search to find points where x > min_x
    let start = points_by_x.partition_point(|&(_, x)| x <= min_x);
    // Check points until x >= max_x
    for &(y, x) in &points_by_x[start..] {
        if x >= max_x {
            break;
        }
        // x is strictly between min_x and max_x, check y
        if min_y < y && y < max_y {
            return true;
        }
    }

    false
}

/// A vertical edge at x coordinate, spanning from min_y to max_y
#[derive(Debug, Clone, Copy)]
struct VerticalEdge {
    x: i64,
    min_y: i64,
    max_y: i64,
}

/// A horizontal edge at y coordinate, spanning from min_x to max_x
#[derive(Debug, Clone, Copy)]
struct HorizontalEdge {
    y: i64,
    min_x: i64,
    max_x: i64,
}

/// Check if any edge from the sorted lists crosses the interior of the rectangle.
/// Uses binary search to find candidate edges efficiently.
fn any_edge_crosses_interior(
    r1: Point,
    r2: Point,
    vertical_edges: &[VerticalEdge],
    horizontal_edges: &[HorizontalEdge],
) -> bool {
    let min_y = r1.y.min(r2.y);
    let max_y = r1.y.max(r2.y);
    let min_x = r1.x.min(r2.x);
    let max_x = r1.x.max(r2.x);

    // Binary search to find vertical edges where x > min_x
    let v_start = vertical_edges.partition_point(|e| e.x <= min_x);
    // Check vertical edges until x >= max_x
    for edge in &vertical_edges[v_start..] {
        if edge.x >= max_x {
            break;
        }
        // x is strictly between min_x and max_x, check y-interval overlap
        if edge.min_y < max_y && edge.max_y > min_y {
            return true;
        }
    }

    // Binary search to find horizontal edges where y > min_y
    let h_start = horizontal_edges.partition_point(|e| e.y <= min_y);
    // Check horizontal edges until y >= max_y
    for edge in &horizontal_edges[h_start..] {
        if edge.y >= max_y {
            break;
        }
        // y is strictly between min_y and max_y, check x-interval overlap
        if edge.min_x < max_x && edge.max_x > min_x {
            return true;
        }
    }

    false
}

fn part2(coords: &[(i64, i64)]) -> u64 {
    let points = preprocess_to_points(&coords);

    // Precompute edges separated by orientation and sorted for binary search
    let mut vertical_edges: Vec<VerticalEdge> = Vec::new();
    let mut horizontal_edges: Vec<HorizontalEdge> = Vec::new();

    for i in 0..coords.len() {
        let e1 = coords[i];
        let e2 = coords[(i + 1) % coords.len()];

        if e1.1 == e2.1 {
            // Vertical edge
            vertical_edges.push(VerticalEdge {
                x: e1.1,
                min_y: e1.0.min(e2.0),
                max_y: e1.0.max(e2.0),
            });
        } else {
            // Horizontal edge
            horizontal_edges.push(HorizontalEdge {
                y: e1.0,
                min_x: e1.1.min(e2.1),
                max_x: e1.1.max(e2.1),
            });
        }
    }

    // Sort by the fixed coordinate for binary search
    vertical_edges.sort_by_key(|e| e.x);
    horizontal_edges.sort_by_key(|e| e.y);

    // Points sorted by x-coordinate for binary search in within check
    let mut points_by_x: Vec<(i64, i64)> = coords.to_vec();
    points_by_x.sort_by_key(|&(_, x)| x);

    let mut best: u64 = 0;
    for (i, &r1) in points.iter().enumerate() {
        // All point indices except i, sorted by area with r1 - these are possible rectangles we can make.
        let mut other_indices: Vec<usize> = (0..points.len()).filter(|&idx| idx != i).collect();
        other_indices.sort_by_key(|&idx| area(r1, points[idx]));

        // Step through possible rectangles (r1, r2) from biggest to smallest...
        for &r2_idx in other_indices.iter().rev() {
            // ... avoiding double counting rectangles we've already looked at...
            if r2_idx <= i {
                continue;
            }
            let r2 = points[r2_idx];

            // ...and ensuring their surround inner area points toward each other.
            if !inner_dirs_point_towards(r1, r2) || !inner_dirs_point_towards(r2, r1) {
                continue;
            }

            let rect_area = area(r1, r2);
            if rect_area <= best {
                // If this is worse than our current best, we can skip to the next r1;
                // the remaining rectangles with this r1 will have smaller area.
                break;
            }

            // Check if any vertex is strictly within our rectangle.
            if any_point_within(r1, r2, &points_by_x) {
                continue;
            }

            // Check whether _any_ polygon edge crosses into our rectangle.
            if any_edge_crosses_interior(r1, r2, &vertical_edges, &horizontal_edges) {
                continue;
            }

            best = rect_area;
            break;
        }
    }
    best
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let filename = &args[1];
    let coords = parse_input(filename).expect("Failed to parse input");

    let start = Instant::now();
    let best = part2(&coords);
    let elapsed = start.elapsed();
    println!("Part 2: {} in {:?}", best, elapsed);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_any_point_within() {
        // Point (3, 7) should be within rectangle from (1, 2) to (5, 11)
        let r1 = Point::new(1, 2);
        let r2 = Point::new(5, 11);

        // Points sorted by x: (y, x)
        let inside_point = vec![(3, 7)];
        assert!(any_point_within(r1, r2, &inside_point));

        // Points on the edge are not considered within the rectangle
        let on_edge_point = vec![(1, 7)];
        assert!(!any_point_within(r1, r2, &on_edge_point));

        // Empty list should return false
        let empty: Vec<(i64, i64)> = vec![];
        assert!(!any_point_within(r1, r2, &empty));
    }

    #[test]
    fn test_preprocess_clockwise_square() {
        // Clockwise square: (0,0) -> (0,2) -> (2,2) -> (2,0)
        //   (0,0)---(0,2)
        //     |       |
        //   (2,0)---(2,2)
        let coords = vec![(0, 0), (0, 2), (2, 2), (2, 0)];
        let points = preprocess_to_points(&coords);

        assert_eq!(points.len(), 4);

        // All corners should be convex (only 1 inner diagonal each)
        for point in &points {
            let inner_count = point.inner_diags.iter().filter(|&&b| b).count();
            assert_eq!(
                inner_count, 1,
                "Convex corner should have exactly 1 inner diagonal"
            );
        }

        // Check specific corners:
        // Clockwise square: (0,0) -> (0,2) -> (2,2) -> (2,0)
        assert!(points[0].inner_diags[DiagDir::DownRight as usize]);
        assert!(points[1].inner_diags[DiagDir::DownLeft as usize]);
        assert!(points[2].inner_diags[DiagDir::UpLeft as usize]);
        assert!(points[3].inner_diags[DiagDir::UpRight as usize]);
    }

    #[test]
    fn test_preprocess_counterclockwise_square() {
        // Counter-clockwise square: (0,0) -> (2,0) -> (2,2) -> (0,2)
        let coords = vec![(0, 0), (2, 0), (2, 2), (0, 2)];
        let points = preprocess_to_points(&coords);

        assert_eq!(points.len(), 4);

        // All corners should be convex (only 1 inner diagonal each)
        for point in &points {
            let inner_count = point.inner_diags.iter().filter(|&&b| b).count();
            assert_eq!(
                inner_count, 1,
                "Convex corner should have exactly 1 inner diagonal"
            );
        }

        // (0,0) -> (2,0) -> (2,2) -> (0,2)
        assert!(points[0].inner_diags[DiagDir::DownRight as usize]);
        assert!(points[1].inner_diags[DiagDir::UpRight as usize]);
        assert!(points[2].inner_diags[DiagDir::UpLeft as usize]);
        assert!(points[3].inner_diags[DiagDir::DownLeft as usize]);
    }

    #[test]
    fn test_preprocess_l_shape_with_concave() {
        // L-shape (clockwise) with one concave corner at (1,1):
        //   x:  0   1   2
        //   y:
        //   0:  *-------*
        //       |       |
        //   1:  *---*   |
        //           |   |
        //   2:      *---*
        //
        // Clockwise: (0,0) -> (0,2) -> (2,2) -> (2,1) -> (1,1) -> (1,0) -> back
        let coords = vec![(0, 0), (0, 2), (2, 2), (2, 1), (1, 1), (1, 0)];
        let points = preprocess_to_points(&coords);

        assert_eq!(points.len(), 6);

        // (1,1) at index 4 is a concave corner - should have 3 inner diagonals
        let concave_point = &points[4]; // (1, 1)
        let inner_count = concave_point.inner_diags.iter().filter(|&&b| b).count();
        assert_eq!(
            inner_count, 3,
            "Concave corner should have 3 inner diagonals"
        );

        // (0,0) -> (0,2) -> (2,2) -> (2,1) -> (1,1) -> (1,0)
        assert!(points[0].inner_diags[DiagDir::DownRight as usize]);
        assert!(points[1].inner_diags[DiagDir::DownLeft as usize]);
        assert!(points[2].inner_diags[DiagDir::UpLeft as usize]);
        assert!(points[3].inner_diags[DiagDir::UpRight as usize]);

        // (1, 1) - concave corner
        assert!(points[4].inner_diags[DiagDir::UpLeft as usize]);
        assert!(points[4].inner_diags[DiagDir::UpRight as usize]);
        assert!(points[4].inner_diags[DiagDir::DownRight as usize]);

        // (1, 0)
        assert!(points[5].inner_diags[DiagDir::UpRight as usize]);
    }
}
