#[cfg(feature = "high_precision")]
pub use rug::Float as Real;

#[cfg(not(feature = "high_precision"))]
pub type Real = f64;

#[cfg(not(feature = "high_precision"))]
#[allow(unused_imports)]
use std::str::FromStr;

use crate::numeric::precision::Precision;
use crate::numeric::traits::{Owned, ParseNumeric};
use anyhow::{Result, anyhow};

#[cfg(feature = "high_precision")]
use rug::{Float, ops::CompleteRound};

#[cfg(feature = "high_precision")]
impl ParseNumeric for Float {
    fn parse(s: &str, precision: Option<Precision>) -> Result<Self> {
        let bits = precision.unwrap_or(Precision::Bits(53)).to_bits();
        Ok(Float::parse(s)
            .map_err(|_| anyhow!("Failed to parse '{}' as Float", s))?
            .complete(bits))
    }
}

#[cfg(not(feature = "high_precision"))]
impl ParseNumeric for f64 {
    fn parse(s: &str, _: Option<Precision>) -> Result<Self> {
        s.parse::<f64>()
            .map_err(|e| anyhow!("Failed to parse '{}' as f64: {}", s, e))
    }
}

#[cfg(feature = "high_precision")]
impl Owned for Float {
    fn owned(&self) -> Self {
        self.clone()
    }
}

#[cfg(not(feature = "high_precision"))]
impl Owned for f64 {
    fn owned(&self) -> Self {
        *self
    }
}
