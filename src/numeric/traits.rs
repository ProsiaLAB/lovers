use crate::numeric::precision::Precision;
use anyhow::Result;
#[allow(unused_imports)]
#[cfg(feature = "high_precision")]
use rug::ops::Pow as RugPow;

/// Generic trait for parsing a numeric type from a string
pub trait ParseNumeric: Sized {
    fn parse(s: &str, precision: Option<Precision>) -> Result<Self>;
}

/// Trait to abstract over cheap vs. expensive clones
pub trait Owned {
    fn owned(&self) -> Self;
}

pub trait FromFloat {
    fn from_f64(val: f64, precision: u32) -> Self;
}

#[cfg(feature = "high_precision")]
impl FromFloat for rug::Float {
    fn from_f64(val: f64, precision: u32) -> Self {
        rug::Float::with_val(precision, val)
    }
}

#[cfg(not(feature = "high_precision"))]
impl FromFloat for f64 {
    fn from_f64(val: f64, _: u32) -> Self {
        val
    }
}

pub trait Pow {
    fn pow_u32(&self, exp: u32, prec: u32) -> Self;
}

#[cfg(not(feature = "high_precision"))]
impl Pow for f64 {
    fn pow_u32(&self, exp: u32, _: u32) -> Self {
        self.powi(exp as i32)
    }
}

#[cfg(feature = "high_precision")]
impl Pow for rug::Float {
    fn pow_u32(&self, exp: u32, prec: u32) -> Self {
        use rug::ops::CompleteRound;

        self.pow(exp).complete(prec)
    }
}
