use std::fmt;
use std::str::FromStr;

use anyhow::anyhow;
use anyhow::{Error, Result};

pub struct Model {}

pub enum Rheology {
    Fluid,
    Elastic,
    Maxwell,
    Newtonian,
    Kelvin,
    Burgers,
    Andrade,
}

// Define a custom error type for parsing failures

impl FromStr for Rheology {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        match s.to_lowercase().as_str() {
            // Convert to lowercase for case-insensitive matching
            "fluid" => Ok(Rheology::Fluid),
            "elastic" => Ok(Rheology::Elastic),
            "maxwell" => Ok(Rheology::Maxwell),
            "newtonian" => Ok(Rheology::Newtonian),
            "kelvin" => Ok(Rheology::Kelvin),
            "burgers" => Ok(Rheology::Burgers),
            "andrade" => Ok(Rheology::Andrade),
            _ => Err(anyhow!("Invalid rheology: '{}'", s)), // If no match, return an error
        }
    }
}

// Optional: Implement Display for Rheology to easily convert it back to a string
impl fmt::Display for Rheology {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Rheology::Fluid => write!(f, "Fluid"),
            Rheology::Elastic => write!(f, "Elastic"),
            Rheology::Maxwell => write!(f, "Maxwell"),
            Rheology::Newtonian => write!(f, "Newtonian"),
            Rheology::Kelvin => write!(f, "Kelvin"),
            Rheology::Burgers => write!(f, "Burgers"),
            Rheology::Andrade => write!(f, "Andrade"),
        }
    }
}

pub fn get_salzer_weights() {}
