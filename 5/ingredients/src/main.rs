use std::env;
use std::error::Error;
use std::fs;

fn parse_input(filename: &str) -> Result<(Vec<(u64, u64)>, Vec<u64>), Box<dyn Error>> {
    let contents = fs::read_to_string(filename)?;
    let mut sections = contents.split("\n\n");

    // Parse ranges (first section)
    let ranges: Vec<(u64, u64)> = sections
        .next()
        .unwrap()
        .lines()
        .map(|line| {
            let (start, end) = line.split_once('-').unwrap();
            (start.parse().unwrap(), end.parse().unwrap())
        })
        .collect();

    // Parse individual numbers (second section)
    let numbers: Vec<u64> = sections
        .next()
        .unwrap()
        .lines()
        .map(|line| line.parse().unwrap())
        .collect();

    Ok((ranges, numbers))
}

fn part1() {
    let args: Vec<String> = env::args().collect();
    let filename = &args[1];
    let (mut ranges, numbers) = parse_input(filename).expect("Failed to parse input");
    ranges.sort_by_key(|&(a, _)| a);

    let count = numbers
        .into_iter()
        .map(|num| {
            for &(start, end) in &ranges {
                if num >= start && num <= end {
                    return 1;
                }
                if num < start {
                    return 0;
                }
            }
            0
        })
        .sum::<u64>();

    println!("Part 1: {count}");
}

fn part2() {
    let args: Vec<String> = env::args().collect();
    let filename = &args[1];
    let (mut ranges, _) = parse_input(filename).expect("Failed to parse input");
    ranges.sort_by_key(|&(a, _)| a);

    let accum_span = ranges
        .into_iter()
        .fold((0u64, 0u64), |(span, last_end), (start, end)| {
            if end <= last_end {
                (span, last_end)
            } else if start > last_end {
                (span + end - start + 1, end)
            } else {
                (span + end - last_end, end)
            }
        })
        .0;

    println!("Part 2: {accum_span}");
}

fn main() {
    part1();
    part2();
}
