#[cfg(feature = "high_precision")]
pub use rug::Integer;

#[cfg(not(feature = "high_precision"))]
pub type Integer = i64;
