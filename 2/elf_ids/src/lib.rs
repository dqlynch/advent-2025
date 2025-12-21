use std::error::Error;
use std::fs;

pub fn parse_input(filename: &str) -> Result<Vec<(u64, u64)>, Box<dyn Error>> {
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

pub fn count_digits(n: u64) -> u32 {
    if n == 0 { 1 } else { n.ilog10() + 1 }
}

/// Split a single range if it crosses a digit boundary.
/// E.g., (5, 114) -> [(5, 9), (10, 114)]
pub fn split_across_digits((start, end): (u64, u64)) -> Vec<(u64, u64)> {
    if count_digits(start) != count_digits(end) {
        vec![
            (start, 10u64.pow(count_digits(start)) - 1),
            (10u64.pow(count_digits(start)), end),
        ]
    } else {
        vec![(start, end)]
    }
}

/// Recursively split ranges until all ranges have uniform digit counts.
pub fn split_all_across_digits(mut ranges: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    loop {
        let len = ranges.len();
        ranges = ranges.into_iter().flat_map(split_across_digits).collect();
        if ranges.len() == len {
            return ranges;
        }
    }
}
