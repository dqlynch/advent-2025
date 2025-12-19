use std::env;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};

fn parse_input(filename: &String) -> Result<Vec<i32>, Box<dyn Error>> {
    let f = File::open(filename)?;
    let f = BufReader::new(f);
    let mut data = vec![];

    for line in f.lines() {
        let line = line?;
        let mult = if line[..1].parse::<char>()? == 'L' {
            -1
        } else {
            1
        };

        let val = line[1..].parse::<i32>()?;

        data.push(mult * val)
    }

    Ok(data)
}

fn count_0_vals(instructions: &Vec<i32>, start: i32) -> u32 {
    let mut val = start;
    let mut count = 0;

    for inst in instructions {
        val = (val + inst).rem_euclid(100); // Rust's `%` is remainder, not modulo
        if val == 0 {
            count += 1;
        }
    }

    count
}

fn count_0_passes(instructions: &Vec<i32>, start: i32) -> i32 {
    let mut val = start;
    let mut count = 0;

    for inst in instructions {
        count += match *inst {
            n if n >= 0 => (n + val) / 100,
            // If we're turning left, just pretend we're turning right
            n => {
                let val = (100 - val) % 100;
                let n = -n;
                (n + val) / 100
            }
        };

        val = (val + inst).rem_euclid(100);
    }

    count
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let filename = &args[1];
    let instructions = parse_input(filename).expect("Should be able to read the given filename");
    let count = count_0_vals(&instructions, 50);
    println!("Part 1 on {filename}: {count}");

    let count = count_0_passes(&instructions, 50);
    println!("Part 2 on {filename}: {count}");
}
