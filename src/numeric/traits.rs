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

pub trait PowI32 {
    fn pow_i32(&self, exp: i32, prec: u32) -> Self;
}

#[cfg(not(feature = "high_precision"))]
impl PowI32 for f64 {
    fn pow_i32(&self, exp: i32, _: u32) -> Self {
        self.powi(exp)
    }
}

#[cfg(feature = "high_precision")]
impl PowI32 for rug::Float {
    fn pow_i32(&self, exp: i32, prec: u32) -> Self {
        use rug::ops::CompleteRound;

        self.pow(exp).complete(prec)
    }
}

pub trait PowReal<Rhs = Self> {
    fn pow_real(&self, exp: Rhs) -> Self;
}

#[cfg(not(feature = "high_precision"))]
impl PowReal for f64 {
    fn pow_real(&self, exp: f64) -> Self {
        self.powf(exp)
    }
}

#[cfg(feature = "high_precision")]
impl PowReal for rug::Float {
    fn pow_real(&self, exp: rug::Float) -> Self {
        use rug::ops::Pow;

        self.clone().pow(exp)
    }
}

pub trait Floor {
    fn floor(&self) -> Self;
}

#[cfg(not(feature = "high_precision"))]
impl Floor for f64 {
    fn floor(&self) -> Self {
        f64::floor(*self)
    }
}

#[cfg(feature = "high_precision")]
impl Floor for rug::Float {
    fn floor(&self) -> Self {
        self.clone().round()
    }
}
