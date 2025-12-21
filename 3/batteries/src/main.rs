use std::env;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};

fn parse_input(filename: &str) -> Result<Vec<Vec<u32>>, Box<dyn Error>> {
    let f = File::open(filename)?;
    let f = BufReader::new(f);
    let mut data = vec![];

    for line in f.lines() {
        let line = line?;
        data.push(
            line.chars()
                .map(|c| c.to_digit(10).unwrap() as u32)
                .collect::<Vec<u32>>(),
        );
    }

    Ok(data)
}

fn eval_bank_n(bank: Vec<u32>, n: u32) -> u64 {
    let mut digits = vec![];

    let mut lbound = 0;
    for rbound in (0..n).rev() {
        let (i, val) = bank[lbound..bank.len() - rbound as usize]
            .iter()
            .enumerate()
            .rev()
            .max_by_key(|(_, val)| *val)
            .unwrap();
        lbound += i + 1;
        digits.push(*val);
    }
    digits
        .iter()
        .enumerate()
        .map(|(i, val)| 10u64.pow(n - i as u32 - 1) * *val as u64)
        .sum()
}

fn part2() {
    let args: Vec<String> = env::args().collect();
    let filename = &args[1];
    let battery_banks = parse_input(filename).expect("Should have been able to read input");
    let result = battery_banks
        .into_iter()
        .map(|bank| eval_bank_n(bank, 12))
        .sum::<u64>();
    println!("Part 2: {result}");
}

fn part1() {
    let args: Vec<String> = env::args().collect();
    let filename = &args[1];
    let battery_banks = parse_input(filename).expect("Should have been able to read input");
    let result = battery_banks
        .into_iter()
        .map(|bank| eval_bank_n(bank, 2))
        .sum::<u64>();

    println!("Part 1: {result}");
}

fn main() {
    part1();
    part2();
}
