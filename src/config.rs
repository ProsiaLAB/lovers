use std::fs;

use anyhow::Result;
use serde::Deserialize;

use crate::numeric::Precision;

#[derive(Debug, Deserialize)]
pub struct Config {
    /// Order of the Gaver sequence
    pub order: usize,
    /// Number of digits of precision
    pub precision: Precision,
    /// Love number kind
    pub kind: Kind,
    /// Minimum degree
    pub min_degree: usize,
    /// Maximum degree
    pub max_degree: usize,
    /// Degree step
    pub degree_step: usize,
    /// Time scale
    pub time_scale: TimeScale,
    /// Time points (minus one)
    pub time_points: usize,
    /// Time range
    pub time_range: (f64, f64),
    /// Load function
    pub load_function: LoadFunction,
    /// Ramp length
    pub ramp_length: f64,
    /// Number of layers
    pub num_layers: usize,
    /// Models
    pub models: Vec<String>,
    /// Log file
    pub log_file: String,
    /// Output type
    pub output_type: OutputType,
    /// Output format
    pub output_format: OutputFormat,
    /// Output files
    pub output_files: Vec<String>,
}

impl Config {
    pub fn load(fname: &str) -> Result<Config> {
        let contents = fs::read_to_string(fname)?;
        let config: Config = toml::from_str(&contents)?;
        Ok(config)
    }
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Kind {
    Loading,
    Tidal,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum TimeScale {
    Log,
    Linear,
    External,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum LoadFunction {
    Step,
    Ramp,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum OutputType {
    Real,
    Complex,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum OutputFormat {
    Time,
    Degree,
}
