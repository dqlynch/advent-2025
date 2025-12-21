use std::env;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};

fn parse_input(filename: &str) -> Result<Vec<Vec<u8>>, Box<dyn Error>> {
    let f = File::open(filename)?;
    let f = BufReader::new(f);
    let mut data = vec![];

    for line in f.lines() {
        let line = line?;
        data.push(
            line.chars()
                .map(|c| match c {
                    '.' => 0,
                    '@' => 1,
                    _ => panic!("Invalid char"),
                })
                .collect::<Vec<u8>>(),
        );
    }

    Ok(data)
}

fn part1() {
    let args: Vec<String> = env::args().collect();
    let filename = &args[1];
    let mut grid = parse_input(filename).expect("Should have been able to read input");

    let h = grid.len();

    let mut sum: u64 = 0;
    loop {
        let mut pass_sum = sum;
        for i in 0..h {
            let w = grid[i].len();
            for j in 0..w {
                if grid[i][j] == 0 {
                    continue;
                }

                let mut adj: u64 = 0;
                for y in 0.max(i as isize - 1)..h.min(i + 2) as isize {
                    for x in 0.max(j as isize - 1)..w.min(j + 2) as isize {
                        if (y as usize, x as usize) != (i, j) {
                            adj += grid[y as usize][x as usize] as u64;
                        }
                        //println!("[{y},{x}] - {}", grid[y as usize][x as usize]);
                    }
                }
                pass_sum += match adj < 4 {
                    true => {
                        grid[i][j] = 0;
                        1
                    }
                    false => 0,
                };
            }
        }
        if pass_sum == sum {
            break;
        }
        sum = pass_sum;
    }

    println!("Sum: {sum}");
}

fn main() {
    part1();
}
