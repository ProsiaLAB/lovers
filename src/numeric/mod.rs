pub mod combinatorics;
pub mod complex;
pub mod integer;
pub mod precision;
pub mod real;
pub mod special;
pub mod traits;

// Publicly expose types conditionally
#[cfg(feature = "high_precision")]
pub use complex::Complex;
#[cfg(not(feature = "high_precision"))]
pub use complex::Complex;

#[cfg(feature = "high_precision")]
pub use integer::Integer;
#[cfg(not(feature = "high_precision"))]
pub use integer::Integer;

pub use precision::Precision;

#[cfg(feature = "high_precision")]
pub use real::Real;
#[cfg(not(feature = "high_precision"))]
pub use real::Real;
