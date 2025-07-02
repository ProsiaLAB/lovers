use std::f64::consts::PI;
use std::fmt;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;
use std::str::FromStr;

use anyhow::anyhow;
use anyhow::{Error, Result};

use crate::config::Config;
use crate::config::TimeScale;
use crate::constants::GRAVITATIONAL_CONSTANT;
use crate::numeric::Precision;
use crate::numeric::Real;
use crate::numeric::combinatorics::{Combinations, Factorial};
use crate::numeric::special::GammaFn;
use crate::numeric::traits::{FromFloat, Owned, ParseNumeric, PowI32, PowReal};

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
                    * (r1.pow_i32(3, precision.to_bits()) - r2.pow_i32(3, precision.to_bits()))
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
                        / r.pow_i32(2, precision.to_bits())
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
        let r0 = self.radii[self.nlayers].owned();
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
        let eta0 = mu0.owned() * t0.owned();
        let mass0 = rho0.owned() * r0.pow_i32(3, precision.to_bits());

        self.radii = self.radii.iter().map(|r| r / r0.owned()).collect();
        self.densities = self.densities.iter().map(|d| d / rho0.owned()).collect();
        self.rigidities = self.rigidities.iter().map(|r| r / mu0.owned()).collect();
        self.viscosities = self.viscosities.iter().map(|v| v / eta0.owned()).collect();

        let grav_normalized = Real::from_f64(GRAVITATIONAL_CONSTANT, precision.to_bits())
            * rho0.pow_i32(2, precision.to_bits())
            * r0.pow_i32(2, precision.to_bits())
            / mu0;
        self.gravity = self
            .gravity
            .iter()
            .map(|g| {
                g * r0.pow_i32(2, precision.to_bits()) * (grav_normalized.owned() / mass0.owned())
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

pub fn get_time_steps(config: &Config) -> Result<Vec<Real>> {
    match config.time_scale {
        TimeScale::Linear => {
            let m1 = config.time_range.0;
            let m2 = config.time_range.1;
            let t1 = Real::from_f64(10.0, config.precision.to_bits())
                .pow_i32(m1, config.precision.to_bits());
            let t2 = Real::from_f64(10.0, config.precision.to_bits())
                .pow_i32(m2, config.precision.to_bits());
            if config.time_points == 0 {
                Ok(vec![t1])
            } else {
                let dt = (t2 - t1.owned())
                    / Real::from_f64(config.time_points as f64, config.precision.to_bits());
                Ok((0..config.time_points + 1)
                    .map(|i| {
                        t1.owned()
                            + dt.owned()
                                * (Real::from_f64(i as f64, config.precision.to_bits())
                                    - Real::from_f64(1.0, config.precision.to_bits()))
                    })
                    .collect::<Vec<Real>>())
            }
        }
        TimeScale::Log => {
            let a1 = Real::from_f64(config.time_range.0.into(), config.precision.to_bits());
            let a2 = Real::from_f64(config.time_range.1.into(), config.precision.to_bits());
            if config.time_points == 0 {
                let t1 = Real::from_f64(10.0, config.precision.to_bits()).pow_real(a1);
                Ok(vec![t1])
            } else {
                let da = (a2 - a1.owned())
                    / Real::from_f64(config.time_points as f64, config.precision.to_bits());
                Ok((0..config.time_points + 1)
                    .map(|i| {
                        a1.owned()
                            + da.owned()
                                * (Real::from_f64(i as f64, config.precision.to_bits())
                                    - Real::from_f64(1.0, config.precision.to_bits()))
                    })
                    .map(|a| Real::from_f64(10.0, config.precision.to_bits()).pow_real(a))
                    .collect::<Vec<Real>>())
            }
        }
        TimeScale::External => {
            let path = Path::new("time_steps.dat");
            let file = File::open(path)?;
            let reader = BufReader::new(file);

            let mut time_steps = Vec::new();
            for line in reader.lines() {
                let line = line?;
                let time_step: Real = ParseNumeric::parse(line.trim(), Some(config.precision))?;
                time_steps.push(time_step);
            }

            Ok(time_steps)
        }
    }
}

pub fn get_salzer_weights(order: u32, precision: Precision) -> Vec<Real> {
    let m = order;
    let mut zeta = Vec::new();
    for k in 0..(2 * m) {
        let j1 = ((k as f64 + 1.0) / 2.0).floor() as u32;
        let j2 = k.min(m);
        zeta.push(Real::from_f64(0.0, precision.to_bits()));
        for j in j1..j2 {
            let fattm = m.factorial();
            let q1 = m.combinations(j);
            let q2 = (2 * j).combinations(j);
            let q3 = j.combinations(k - j);
            zeta[k as usize] +=
                (j as f64).powf(m as f64 + 1.0) / fattm as f64 * q1 as f64 * q2 as f64 * q3 as f64;
            if k % 2 != 0 {
                zeta[k as usize] = -zeta[k as usize].owned();
            }
        }
    }
    zeta
}
