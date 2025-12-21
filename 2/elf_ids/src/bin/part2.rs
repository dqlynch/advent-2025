use std::{collections::HashSet, env};

use elf_ids::{count_digits, parse_input, split_all_across_digits};

/// Get repeats of size `k` that reside within [start, end] (inclusive).
/// A "repeat" of size k is a number whose digits are a repeated sequence of length k:
///     ex: 45454545 is a repeat of size 2 (45), but also of size 4 (4545).
///
/// The number of digits of start and end must match.
/// The number of digits must be divisible by k.
/// O(10^k), in the degenerate case (range of 0..n).
fn get_repeats_of_size((start, end): (u64, u64), k: u32) -> Vec<u64> {
    // We'll call the repeated element of our repeat a "part" - the "part that's repeated"
    // We can express a repeat as 1 * part + 10**k * part + 10**2k * part + ...
    //      ex: 123_123_123 is 123 + 123000 + 123000000
    // In other words, it's 1111... * part, in base 10**k.
    //
    // We'll explore the possible repeats by generating from the enumeration of
    // our most significant part, because it always has the most restrictive range.
    let base = 10u64.pow(k);
    let num_parts = count_digits(start) / k;

    // Most significant parts, i.e. the leading k digits
    let (start_msp, end_msp) = (
        start / base.pow(num_parts - 1),
        end / base.pow(num_parts - 1),
    );

    (start_msp..=end_msp)
        .map(|msp| {
            // Generate a repeat from the start part
            (0..num_parts).map(|i| base.pow(i) * msp).sum::<u64>()
        })
        .filter(|val| *val >= start && *val <= end)
        .collect()
}

fn get_repeats((start, end): (u64, u64)) -> HashSet<u64> {
    let digits = count_digits(start); // We know start and end have the same digits
    (1..=digits / 2)
        .filter(|k| digits % k == 0)
        .flat_map(|k| get_repeats_of_size((start, end), k))
        .collect()
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let filename = &args[1];

    let ranges = parse_input(filename).expect("Should be able to read the given filename");

    // Split across digit boundaries: [5-14] -> [5-9, 10-14]
    let ranges = split_all_across_digits(ranges);

    let sum = ranges.into_iter().flat_map(get_repeats).sum::<u64>();

    println!("Part 2 sum: {sum}");
}
