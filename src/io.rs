use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;

use anyhow::Result;

use crate::core::Rheology;
use crate::numeric::Precision;
use crate::numeric::Real;
use crate::numeric::special::GammaFn;
use crate::numeric::traits::Owned;
use crate::numeric::traits::ParseNumeric;

pub fn read_model(fname: &str, precision: Precision) -> Result<()> {
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
        let radius = <Real as ParseNumeric>::parse(parts[0], Some(precision))?;
        let density = <Real as ParseNumeric>::parse(parts[1], Some(precision))?;
        let rigidity = <Real as ParseNumeric>::parse(parts[2], Some(precision))?;
        let viscosity = <Real as ParseNumeric>::parse(parts[3], Some(precision))?;
        let rheology = parts[4].parse::<Rheology>()?;

        match rheology {
            Rheology::Burgers => {
                let p1 = <Real as ParseNumeric>::parse(parts[5], Some(precision))?;
                let p2 = <Real as ParseNumeric>::parse(parts[6], Some(precision))?;
                let par = (p1, p2);
                println!("{:?}", par);
            }
            Rheology::Andrade => {
                let p1: Real = ParseNumeric::parse(parts[5], Some(precision))?;
                let p2 = GammaFn::gamma(&(p1.owned() + 1.0));
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
