use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;

use anyhow::Result;

use crate::core::Rheology;
use crate::utils::gamma::gamma;

pub fn read_model(fname: &str) -> Result<()> {
    let path = Path::new(fname);
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Skip comment lines
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('!') {
            continue;
        }
        // Split at whitespace
        let parts: Vec<&str> = line.split_whitespace().collect();
        let radius = parts[0].parse::<f64>()?;
        let density = parts[1].parse::<f64>()?;
        let rigidity = parts[2].parse::<f64>()?;
        let viscosity = parts[3].parse::<f64>()?;
        let rheology = parts[4].parse::<Rheology>()?;

        match rheology {
            Rheology::Burgers => {
                let par = (parts[5].parse::<f64>()?, parts[6].parse::<f64>()?);
                println!("{:?}", par);
            }
            Rheology::Andrade => {
                let p1 = parts[5].parse::<f64>()?;
                let p2 = gamma(p1 + 1.0);
                let par = (p1, p2);
                println!("{:?}", par);
            }
            _ => {
                // Do nothing
            }
        }
    }

    Ok(())
}
