use crate::numeric::precision::Precision;
use anyhow::Result;

/// Generic trait for parsing a numeric type from a string
pub trait ParseNumeric: Sized {
    fn parse(s: &str, precision: Option<Precision>) -> Result<Self>;
}

/// Trait to abstract over cheap vs. expensive clones
pub trait Owned {
    fn owned(&self) -> Self;
}
