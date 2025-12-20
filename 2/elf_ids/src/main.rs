use std::env;
use std::error::Error;
use std::fs;

fn parse_input(filename: &String) -> Result<Vec<(u64, u64)>, Box<dyn Error>> {
    let contents = fs::read_to_string(filename)?;

    let data: Vec<(u64, u64)> = contents
        .split(',')
        .map(|range| {
            let (start, end) = range.trim().split_once('-').unwrap();
            (start.parse().unwrap(), end.parse().unwrap())
        })
        .collect();

    Ok(data)
}

fn count_digits(n: u64) -> u32 {
    if n == 0 { 1 } else { n.ilog10() + 1 }
}

/// Split a single range if it crosses a digit boundary.
/// E.g., (5, 114) -> [(5, 9), (10, 114)]
fn split_across_digits((start, end): (u64, u64)) -> Vec<(u64, u64)> {
    if count_digits(start) != count_digits(end) {
        vec![
            (start, 10u64.pow(count_digits(start)) - 1),
            (10u64.pow(count_digits(start)), end),
        ]
    } else {
        vec![(start, end)]
    }
}

fn split_all_across_digits(mut ranges: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    loop {
        let len = ranges.len();
        ranges = ranges.into_iter().flat_map(split_across_digits).collect();
        if ranges.len() == len {
            return ranges;
        }
    }
}

fn keep_even_digits((start, end): (u64, u64)) -> Option<(u64, u64)> {
    match count_digits(start) % 2 {
        0 => Some((start, end)),
        _ => None,
    }
}

fn subdivide_ranges((start, end): (u64, u64)) -> (u64, u64, u64, u64) {
    // Subdivide our numbers in "half" by digit count
    let digits = count_digits(start);

    let factor = 10u64.pow(digits / 2);
    // (start_front, end_front, start_back, end_back)
    let res = (start / factor, end / factor, start % factor, end % factor);
    dbg!(res);
    res
}

fn split_to_simple(
    (start_front, end_front, start_back, end_back): (u64, u64, u64, u64),
) -> Vec<(u64, u64, u64, u64)> {
    let digits = count_digits(start_front);
    let factor = 10u64.pow(digits);
    let simple_splits = match end_front - start_front {
        0 => vec![(start_front, end_front, start_back, end_back)],

        1 => vec![
            (start_front, start_front, start_back, factor - 1),
            (end_front, end_front, 0, end_back),
        ],

        _ => vec![
            (start_front, start_front, start_back, factor - 1),
            (start_front + 1, end_front - 1, 0, factor - 1),
            (end_front, end_front, 00, end_back),
        ],
    };
    dbg!(&simple_splits);
    simple_splits
}

fn get_repeated_numbers(
    (start_front, end_front, start_back, end_back): (u64, u64, u64, u64),
) -> Vec<u64> {
    // Now, the intersection of our front range and back range gives us every repeated half
    let d = (start_front.max(start_back)..=end_front.min(end_back))
        .map(|half| {
            let digits = count_digits(half);
            half * 10u64.pow(digits) + half
        })
        .collect();
    dbg!(&d);
    d
}

fn part1() {
    let args: Vec<String> = env::args().collect();
    let filename = &args[1];
    let ranges = parse_input(filename).expect("Should be able to read the given filename");

    // 1. Split across digit boundaries: [5-14] -> [5-9, 10-14]
    // 2. Drop odd-digit ranges (can't have repeated numbers)
    // 3. Subdivide to digit halves: [1050-2150] -> [10-21/50-50]
    // 3. Split to "simple ranges": [1050-2150] -> [1050-1099, 1100-2099, 2100-2150]
    // 4. Calculate repeated number for each simple range
    let ranges = split_all_across_digits(ranges);
    dbg!(&ranges);

    let sum: u64 = ranges
        .into_iter()
        .filter_map(keep_even_digits)
        .map(subdivide_ranges)
        .flat_map(split_to_simple)
        .flat_map(get_repeated_numbers)
        .sum();

    println!("Part 1 sum: {sum}");
}

fn main() {
    part1();
}
