use anyhow::{Result, anyhow};
#[cfg(feature = "high_precision")]
use rug::ops::CompleteRound;

use crate::types::Real;

/// Trait for parsing strings into the Real type
pub trait ParseReal {
    /// Parse a string into the Real type
    fn parse_real(s: &str, precision: Option<u32>) -> Result<Real>;
}

pub struct Parser;

impl ParseReal for Parser {
    fn parse_real(s: &str, precision: Option<u32>) -> Result<Real> {
        // If precision is not specified, use the default

        let precision = precision.unwrap_or(53); // By default use 53 bits of precision
        #[cfg(feature = "high_precision")]
        {
            // For high precision, use rug::Float
            // You can customize precision here if needed
            Ok(rug::Float::parse(s)
                .map_err(|_| anyhow!("Failed to parse '{}' as high precision float", s))?
                .complete(precision))
        }

        #[cfg(not(feature = "high_precision"))]
        {
            // For standard precision, use f64
            s.parse::<f64>()
                .map_err(|e| anyhow!("Failed to parse '{}' as f64: {}", s, e))
        }
    }
}

pub trait Owned {
    fn owned(&self) -> Self;
}

#[cfg(not(feature = "high_precision"))]
impl Owned for f64 {
    fn owned(&self) -> Self {
        *self
    }
}

#[cfg(feature = "high_precision")]
impl Owned for rug::Float {
    fn owned(&self) -> Self {
        self.clone()
    }
}
