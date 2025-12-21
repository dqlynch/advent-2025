use std::env;

use elf_ids::{count_digits, parse_input, split_all_across_digits};

fn keep_even_digits((start, end): (u64, u64)) -> Option<(u64, u64)> {
    match count_digits(start) % 2 {
        0 => Some((start, end)),
        _ => None,
    }
}

fn subdivide_ranges((start, end): (u64, u64)) -> (u64, u64, u64, u64) {
    // Subdivide our numbers in "half" by digit count
    let digits = count_digits(start);

    let base = 10u64.pow(digits / 2); // as in "base 10" - we're extracting these in "base 10^(k/2)"

    // ((start_front, end_front), (start_back, end_back))
    (start / base, end / base, start % base, end % base)
}

fn split_to_simple(
    (start_front, end_front, start_back, end_back): (u64, u64, u64, u64),
) -> Vec<(u64, u64, u64, u64)> {
    let digits = count_digits(start_front);
    let base = 10u64.pow(digits);
    match end_front - start_front {
        0 => vec![(start_front, end_front, start_back, end_back)],

        1 => vec![
            (start_front, start_front, start_back, base - 1),
            (end_front, end_front, 0, end_back),
        ],

        _ => vec![
            (start_front, start_front, start_back, base - 1),
            (start_front + 1, end_front - 1, 0, base - 1),
            (end_front, end_front, 0, end_back),
        ],
    }
}

fn get_repeated_numbers(
    (start_front, end_front, start_back, end_back): (u64, u64, u64, u64),
) -> Vec<u64> {
    // The intersection of our front range and back range gives us every repeated half
    (start_front.max(start_back)..=end_front.min(end_back))
        .map(|half| {
            let digits = count_digits(half);
            half * 10u64.pow(digits) + half
        })
        .collect()
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let filename = &args[1];
    let ranges = parse_input(filename).expect("Should be able to read the given filename");

    // 1. Split across digit boundaries: [5-14] -> [5-9, 10-14]
    // 2. Drop odd-digit ranges (can't have repeated numbers)
    // 3. Subdivide to digit halves: [1050-2150] -> [10-21/50-50]
    // 4. Split to "simple ranges": [1050-2150] -> [1050-1099, 1100-2099, 2100-2150]
    // 5. Intersect each simple range's front and back ranges:
    //      ex: for range [199000-253999], we have:
    //          front range: [199-254], back range [000-999]
    //          intersection = [199-253]
    //              => 199, 200, 201, ..., 253 -> have repeats
    //              => 199, 200200, ..., 253253 are repeats
    //
    //      ex: tail range [254000-254709] -> [254-254], [000-709]
    //          intersection = [254] => 254254 is a repeat
    //
    let ranges = split_all_across_digits(ranges);

    // Worst case O(sqrt(n)) for each range with upper bound of n
    let sum: u64 = ranges
        .into_iter()
        .filter_map(keep_even_digits)
        .map(subdivide_ranges)
        .flat_map(split_to_simple)
        .flat_map(get_repeated_numbers)
        .sum();

    println!("Part 1 sum: {sum}");

    // If we generalized this algorithm for all k-repeats where 2 <= k <= log_10(n) / 2, then:
    // 1. Our "split to simple" becomes much more complicated - for a k-repeat, we split across
    //      k-1 boundaries. => create O(3^k) = 3(log10(n)) ~= O(n^.477) ~= O(sqrt(n)) splits (less)
    // 2. For each simple range:
    //  * Midde range must enumerate the full range [0, base) => O(n^(1/k)) values
    //  * Boundary ranges have at least one constant => (O(1)) work
    // 3. For each _starting_ range, we have only ONE middle range
    //  * 1 middle range * (O(n^1/k))
    //  * < O(sqrt(n)) boundary ranges * O(1)
    //      => O(n^1/2) > O(n^.477), so this is dominated by the expansion of the k=2 middle range.
    // => Overall complexity is O(sqrt(n)), asymptotically the same as our n=2 case!
    //
    // However, the algorithm gets complex to implement, and I thought of a MUCH easier way!
}
