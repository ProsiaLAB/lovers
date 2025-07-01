use std::f64::consts::PI;
use std::fmt;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;
use std::str::FromStr;

use anyhow::anyhow;
use anyhow::{Error, Result};

use crate::constants::GRAVITATIONAL_CONSTANT;
use crate::numeric::Precision;
use crate::numeric::Real;
use crate::numeric::special::GammaFn;
use crate::numeric::traits::{FromFloat, Owned, ParseNumeric, Pow};

#[derive(Debug)]
pub struct Model {
    pub nlayers: usize,
    pub radii: Vec<Real>,
    pub densities: Vec<Real>,
    pub rigidities: Vec<Real>,
    pub viscosities: Vec<Real>,
    pub pars: Vec<Option<(Real, Real)>>,
    pub gravity: Vec<Real>,
}

impl Model {
    pub fn load(fname: &str, precision: Precision) -> Result<Self> {
        let path = Path::new(fname);
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        let mut radii = Vec::new();
        let mut densities = Vec::new();
        let mut rigidities = Vec::new();
        let mut viscosities = Vec::new();
        let mut pars = Vec::new();

        radii.push(Real::from_f64(0.0, precision.to_bits()));

        // Skip comment lines
        for line in reader.lines() {
            let line = line?;
            if line.starts_with('!') {
                continue;
            }
            // Split at whitespace
            let parts: Vec<&str> = line.split_whitespace().collect();
            let radius: Real = ParseNumeric::parse(parts[0], Some(precision))?;
            let density: Real = ParseNumeric::parse(parts[1], Some(precision))?;
            let rigidity: Real = ParseNumeric::parse(parts[2], Some(precision))?;
            let viscosity: Real = ParseNumeric::parse(parts[3], Some(precision))?;
            let rheology = parts[4].parse::<Rheology>()?;

            let par = match rheology {
                Rheology::Burgers => {
                    let p1: Real = ParseNumeric::parse(parts[5], Some(precision))?;
                    let p2: Real = ParseNumeric::parse(parts[6], Some(precision))?;
                    Some((p1, p2))
                }
                Rheology::Andrade => {
                    let p1: Real = ParseNumeric::parse(parts[5], Some(precision))?;
                    let p2 = GammaFn::gamma(&(p1.owned() + 1.0));
                    Some((p1, p2))
                }
                _ => None,
            };

            radii.push(radius);
            densities.push(density);
            rigidities.push(rigidity);
            viscosities.push(viscosity);
            pars.push(par);
        }

        let nlayers = densities.len();

        // Compute the mass of the planet
        let zero = Real::from_f64(0.0, precision.to_bits());
        let coeff = Real::from_f64(4.0, precision.to_bits())
            / Real::from_f64(3.0, precision.to_bits())
            * PI;
        let mass: Real = densities
            .iter()
            .zip(radii.windows(2).skip(1))
            .map(|(rho, window)| {
                let [r1, r2] = window else {
                    panic!("Expected window of exactly 2 elements");
                };
                rho.owned()
                    * (r1.pow_u32(3, precision.to_bits()) - r2.pow_u32(3, precision.to_bits()))
            })
            .map(|m| m * coeff.owned())
            .fold(zero, |acc, x| acc + x);

        // Compute gravity at the interface boundaries
        let gravity: Vec<Real> = radii
            .iter()
            .enumerate()
            .map(|(i, r)| {
                if i == 0 {
                    Real::from_f64(0.0, precision.to_bits())
                } else {
                    Real::from_f64(GRAVITATIONAL_CONSTANT, precision.to_bits()) * mass.owned()
                        / r.pow_u32(2, precision.to_bits())
                }
            })
            .collect();

        Ok(Self {
            nlayers,
            radii,
            densities,
            rigidities,
            viscosities,
            pars,
            gravity,
        })
    }

    pub fn normalize(&mut self, precision: Precision) {
        let r0 = self.radii[self.nlayers];
        let rho0 = self
            .densities
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .cloned()
            .expect("Empty vector");
        let mu0 = self
            .rigidities
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .cloned()
            .expect("Empty vector");
        let t0 = Real::from_f64(1000.0, precision.to_bits())
            * Real::from_f64(365.25, precision.to_bits())
            * Real::from_f64(24.0, precision.to_bits())
            * Real::from_f64(3600.0, precision.to_bits());
        let eta0 = mu0 * t0;
        let mass0 = rho0 * r0.pow_u32(3, precision.to_bits());

        self.radii = self.radii.iter().map(|r| r / r0).collect();
        self.densities = self.densities.iter().map(|d| d / rho0).collect();
        self.rigidities = self.rigidities.iter().map(|r| r / mu0).collect();
        self.viscosities = self.viscosities.iter().map(|v| v / eta0).collect();

        let grav_normalized = Real::from_f64(GRAVITATIONAL_CONSTANT, precision.to_bits())
            * rho0.pow_u32(2, precision.to_bits())
            * r0.pow_u32(2, precision.to_bits())
            / mu0;
        self.gravity = self
            .gravity
            .iter()
            .map(|g| {
                g * r0.pow_u32(2, precision.to_bits())
                    * (grav_normalized / Real::from_f64(mass0, precision.to_bits()))
            })
            .collect();
    }
}

/// The rheological constitutive law
pub enum Rheology {
    /// Inviscid fluid
    Fluid,
    /// Hooke's rheology
    Elastic,
    /// Maxwell's rheology
    Maxwell,
    /// Newtonian rheology
    Newtonian,
    /// Kelvin rheology
    Kelvin,
    /// Burgers rheology
    Burgers,
    /// Andrade rheology
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
