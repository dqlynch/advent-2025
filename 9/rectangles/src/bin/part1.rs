use std::env;

use rectangles::parse_input;

fn part1(coords: &[(i64, i64)]) -> u64 {
    coords
        .iter()
        .map(|&(y1, x1)| {
            coords
                .iter()
                .map(move |&(y2, x2)| (y2 - y1 + 1).abs() as u64 * (x2 - x1 + 1).abs() as u64)
                .max()
                .unwrap_or(0)
        })
        .max()
        .unwrap_or(0)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let filename = &args[1];
    let coords = parse_input(filename).expect("Failed to parse input");

    let best = part1(&coords);
    println!("Part 1: {}", best);
}
