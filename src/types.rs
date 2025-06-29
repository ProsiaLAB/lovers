use num_complex::Complex64;

#[cfg(feature = "high_precision")]
pub use rug::Float as Real;

#[cfg(feature = "high_precision")]
pub use rug::Integer;

#[cfg(feature = "high_precision")]
pub use rug::Complex;

#[cfg(not(feature = "high_precision"))]
pub type Real = f64;

#[cfg(not(feature = "high_precision"))]
pub type Integer = i64;

#[cfg(not(feature = "high_precision"))]
pub type Complex = Complex64;
