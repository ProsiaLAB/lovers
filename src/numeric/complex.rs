#[cfg(feature = "high_precision")]
pub use rug::Complex;

#[cfg(not(feature = "high_precision"))]
pub use num_complex::Complex64 as Complex;
