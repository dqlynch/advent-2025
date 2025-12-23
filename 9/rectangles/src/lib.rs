use std::error::Error;
use std::fs;

pub fn parse_input(filename: &str) -> Result<Vec<(i64, i64)>, Box<dyn Error>> {
    let contents = fs::read_to_string(filename)?;

    let coords = contents
        .lines()
        .map(|line| {
            let (x, y) = line.split_once(',').unwrap();
            (x.parse().unwrap(), y.parse().unwrap())
        })
        .collect();

    Ok(coords)
}
