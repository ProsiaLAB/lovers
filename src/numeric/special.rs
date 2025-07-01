#[allow(unused_imports)]
use crate::numeric::Real;

#[cfg(not(feature = "high_precision"))]
use crate::utils::gamma::gamma;

pub trait GammaFn {
    fn gamma(&self) -> Self;
}

#[cfg(not(feature = "high_precision"))]
impl GammaFn for f64 {
    fn gamma(&self) -> Self {
        gamma(*self)
    }
}

#[cfg(feature = "high_precision")]
impl GammaFn for Real {
    fn gamma(&self) -> Self {
        self.clone().gamma()
    }
}
