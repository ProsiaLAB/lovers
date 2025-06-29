mod combinatorics {
    use crate::types::Integer;
    use rug::Complete;
    pub trait Factorial {
        fn factorial(self) -> Self;
    }

    impl Factorial for i64 {
        fn factorial(self) -> Self {
            (1..=self).product()
        }
    }

    impl Factorial for rug::Integer {
        fn factorial(self) -> Self {
            let mut result = rug::Integer::from(1);
            let mut i = rug::Integer::from(2);
            while i <= self {
                result *= &i;
                i += 1;
            }
            result
        }
    }

    pub trait Combinations {
        fn combinations(self, k: Self) -> Self;
    }

    impl Combinations for i64 {
        fn combinations(self, k: Self) -> Self {
            if k > self {
                0
            } else {
                (1..=k).fold(1, |acc, val| acc * (self - val + 1) / val)
            }
        }
    }

    impl Combinations for rug::Integer {
        fn combinations(self, k: Self) -> Self {
            if k > self {
                return Integer::from(0);
            }

            let mut acc = Integer::from(1);
            let mut val = Integer::from(1);

            while val <= k {
                // Complete subtraction before adding
                let num = (&self - &val).complete() + 1;
                acc *= num;
                acc /= &val;
                val += 1;
            }

            acc
        }
    }

    pub trait Binomial {
        fn binomial(self, k: Self) -> Self;
    }

    impl Binomial for i64 {
        fn binomial(self, k: Self) -> Self {
            self.combinations(k)
        }
    }

    impl Binomial for rug::Integer {
        fn binomial(self, k: Self) -> Self {
            self.combinations(k)
        }
    }

    pub trait Permutations {
        fn permutations(self, k: Self) -> Self;
    }
    impl Permutations for i64 {
        fn permutations(self, k: Self) -> Self {
            (self - k + 1..=self).product()
        }
    }
    impl Permutations for rug::Integer {
        fn permutations(self, k: Self) -> Self {
            if k > self {
                return Integer::from(0);
            }

            let mut result = Integer::from(1);
            let mut i = (&self - &k).complete() + 1;

            while i <= self {
                result *= &i;
                i += 1;
            }

            result
        }
    }
}
pub mod constants {
    //! Defines mathematical expressions commonly used when computing distribution
    //! values as constants
    //!
    //! This module is directly copied from the `statrs` crate.
    //! The original code is licensed under the MIT license.
    //! It was purely done to avoid a dependency on the `statrs` crate.
    //!
    //! The original code can be found here:
    //! <https://github.com/statrs-dev/statrs/blob/master/src/constants.rs>
    #![allow(clippy::excessive_precision)]

    /// Constant value for `ln(pi)`
    pub const LN_PI: f64 = 1.1447298858494001741434273513530587116472948129153;

    /// Constant value for `ln(2 * sqrt(e / pi))`
    pub const LN_2_SQRT_E_OVER_PI: f64 = 0.6207822376352452223455184457816472122518527279025978;

    /// Constant value for `2 * sqrt(e / pi)`
    pub const TWO_SQRT_E_OVER_PI: f64 = 1.8603827342052657173362492472666631120594218414085755;

    /// Default accuracy for `f64`, equivalent to `0.0 * F64_PREC`
    pub const DEFAULT_F64_ACC: f64 = 0.0000000000000011102230246251565;
}

pub mod special {
    use approx::AbsDiffEq;

    /// Compares if two floats are close via `approx::abs_diff_eq`
    /// using a maximum absolute difference (epsilon) of `acc`.
    ///
    /// This is directly copied from the `statrs` crate.
    /// The original code is licensed under the MIT license.
    /// It was purely done to avoid a dependency on the `statrs` crate.
    ///
    /// The original code can be found here:
    /// <https://github.com/statrs-dev/statrs/blob/master/src/prec.rs>
    pub fn almost_eq(a: f64, b: f64, acc: f64) -> bool {
        if a.is_infinite() && b.is_infinite() {
            return a == b;
        }
        a.abs_diff_eq(&b, acc)
    }
}

pub mod gamma {
    //! Gamma function related utilities
    //!
    //! This module is directly copied from the `statrs` crate.
    //! The original code is licensed under the MIT license.
    //! It was purely done to avoid a dependency on the `statrs` crate.
    //!
    //! The original code can be found here:
    //! <https://github.com/statrs-dev/statrs/blob/master/src/function/gamma.rs>

    #![allow(clippy::excessive_precision)]
    use std::f64::consts::{E, PI};

    use anyhow::Result;
    use anyhow::anyhow;
    use approx::ulps_eq;

    use crate::types::Real;

    use super::constants::{DEFAULT_F64_ACC, LN_2_SQRT_E_OVER_PI, LN_PI, TWO_SQRT_E_OVER_PI};
    use super::special::almost_eq;

    /// Auxiliary variable when evaluating the `gamma_ln` function
    const GAMMA_R: f64 = 10.900511;

    /// Polynomial coefficients for approximating the `gamma_ln` function
    const GAMMA_DK: &[f64] = &[
        2.48574089138753565546e-5,
        1.05142378581721974210,
        -3.45687097222016235469,
        4.51227709466894823700,
        -2.98285225323576655721,
        1.05639711577126713077,
        -1.95428773191645869583e-1,
        1.70970543404441224307e-2,
        -5.71926117404305781283e-4,
        4.63399473359905636708e-6,
        -2.71994908488607703910e-9,
    ];

    /// Computes the logarithm of the gamma function
    /// with an accuracy of 16 floating point digits.
    /// The implementation is derived from
    /// "An Analysis of the Lanczos Gamma Approximation",
    /// Glendon Ralph Pugh, 2004 p. 116
    pub fn ln_gamma(x: f64) -> f64 {
        if x < 0.5 {
            let s = GAMMA_DK
                .iter()
                .enumerate()
                .skip(1)
                .fold(GAMMA_DK[0], |s, t| s + t.1 / (t.0 as f64 - x));

            LN_PI
                - (PI * x).sin().ln()
                - s.ln()
                - LN_2_SQRT_E_OVER_PI
                - (0.5 - x) * ((0.5 - x + GAMMA_R) / E).ln()
        } else {
            let s = GAMMA_DK
                .iter()
                .enumerate()
                .skip(1)
                .fold(GAMMA_DK[0], |s, t| s + t.1 / (x + t.0 as f64 - 1.0));

            s.ln() + LN_2_SQRT_E_OVER_PI + (x - 0.5) * ((x - 0.5 + GAMMA_R) / E).ln()
        }
    }

    /// Computes the gamma function with an accuracy
    /// of 16 floating point digits. The implementation
    /// is derived from "An Analysis of the Lanczos Gamma Approximation",
    /// Glendon Ralph Pugh, 2004 p. 116
    pub fn gamma(x: f64) -> f64 {
        if x < 0.5 {
            let s = GAMMA_DK
                .iter()
                .enumerate()
                .skip(1)
                .fold(GAMMA_DK[0], |s, t| s + t.1 / (t.0 as f64 - x));

            PI / ((PI * x).sin() * s * TWO_SQRT_E_OVER_PI * ((0.5 - x + GAMMA_R) / E).powf(0.5 - x))
        } else {
            let s = GAMMA_DK
                .iter()
                .enumerate()
                .skip(1)
                .fold(GAMMA_DK[0], |s, t| s + t.1 / (x + t.0 as f64 - 1.0));

            s * TWO_SQRT_E_OVER_PI * ((x - 0.5 + GAMMA_R) / E).powf(x - 0.5)
        }
    }

    /// Computes the upper incomplete gamma function
    /// `Gamma(a,x) = int(exp(-t)t^(a-1), t=0..x) for a > 0, x > 0`
    /// where `a` is the argument for the gamma function and
    /// `x` is the lower intergral limit.
    ///
    /// # Panics
    ///
    /// if `a` or `x` are not in `(0, +inf)`
    pub fn gamma_ui(a: f64, x: f64) -> f64 {
        checked_gamma_ui(a, x).unwrap()
    }

    /// Computes the upper incomplete gamma function
    /// `Gamma(a,x) = int(exp(-t)t^(a-1), t=0..x) for a > 0, x > 0`
    /// where `a` is the argument for the gamma function and
    /// `x` is the lower intergral limit.
    ///
    /// # Errors
    ///
    /// if `a` or `x` are not in `(0, +inf)`
    pub fn checked_gamma_ui(a: f64, x: f64) -> Result<f64> {
        checked_gamma_ur(a, x).map(|x| x * gamma(a))
    }

    /// Computes the lower incomplete gamma function
    /// `gamma(a,x) = int(exp(-t)t^(a-1), t=0..x) for a > 0, x > 0`
    /// where `a` is the argument for the gamma function and `x`
    /// is the upper integral limit.
    ///
    ///
    /// # Panics
    ///
    /// if `a` or `x` are not in `(0, +inf)`
    pub fn gamma_li(a: f64, x: f64) -> f64 {
        checked_gamma_li(a, x).unwrap()
    }

    /// Computes the lower incomplete gamma function
    /// `gamma(a,x) = int(exp(-t)t^(a-1), t=0..x) for a > 0, x > 0`
    /// where `a` is the argument for the gamma function and `x`
    /// is the upper integral limit.
    ///
    ///
    /// # Errors
    ///
    /// if `a` or `x` are not in `(0, +inf)`
    pub fn checked_gamma_li(a: f64, x: f64) -> Result<f64> {
        checked_gamma_lr(a, x).map(|x| x * gamma(a))
    }

    /// Computes the upper incomplete regularized gamma function
    /// `Q(a,x) = 1 / Gamma(a) * int(exp(-t)t^(a-1), t=0..x) for a > 0, x > 0`
    /// where `a` is the argument for the gamma function and
    /// `x` is the lower integral limit.
    ///
    /// # Remarks
    ///
    /// Returns `f64::NAN` if either argument is `f64::NAN`
    ///
    /// # Panics
    ///
    /// if `a` or `x` are not in `(0, +inf)`
    pub fn gamma_ur(a: f64, x: f64) -> f64 {
        checked_gamma_ur(a, x).unwrap()
    }

    /// Computes the upper incomplete regularized gamma function
    /// `Q(a,x) = 1 / Gamma(a) * int(exp(-t)t^(a-1), t=0..x) for a > 0, x > 0`
    /// where `a` is the argument for the gamma function and
    /// `x` is the lower integral limit.
    ///
    /// # Remarks
    ///
    /// Returns `f64::NAN` if either argument is `f64::NAN`
    ///
    /// # Errors
    ///
    /// if `a` or `x` are not in `(0, +inf)`
    pub fn checked_gamma_ur(a: f64, x: f64) -> Result<f64> {
        if a.is_nan() || x.is_nan() {
            return Ok(f64::NAN);
        }
        if a <= 0.0 || a == f64::INFINITY {
            return Err(anyhow!("AInvalid"));
        }
        if x <= 0.0 || x == f64::INFINITY {
            return Err(anyhow!("XInvalid"));
        }

        let eps = 0.000000000000001;
        let big = 4503599627370496.0;
        let big_inv = 2.22044604925031308085e-16;

        if x < 1.0 || x <= a {
            return Ok(1.0 - gamma_lr(a, x));
        }

        let mut ax = a * x.ln() - x - ln_gamma(a);
        if ax < -709.78271289338399 {
            return if a < x { Ok(0.0) } else { Ok(1.0) };
        }

        ax = ax.exp();
        let mut y = 1.0 - a;
        let mut z = x + y + 1.0;
        let mut c = 0.0;
        let mut pkm2 = 1.0;
        let mut qkm2 = x;
        let mut pkm1 = x + 1.0;
        let mut qkm1 = z * x;
        let mut ans = pkm1 / qkm1;
        loop {
            y += 1.0;
            z += 2.0;
            c += 1.0;
            let yc = y * c;
            let pk = pkm1 * z - pkm2 * yc;
            let qk = qkm1 * z - qkm2 * yc;

            pkm2 = pkm1;
            pkm1 = pk;
            qkm2 = qkm1;
            qkm1 = qk;

            if pk.abs() > big {
                pkm2 *= big_inv;
                pkm1 *= big_inv;
                qkm2 *= big_inv;
                qkm1 *= big_inv;
            }

            if qk != 0.0 {
                let r = pk / qk;
                let t = ((ans - r) / r).abs();
                ans = r;

                if t <= eps {
                    break;
                }
            }
        }
        Ok(ans * ax)
    }

    /// Computes the lower incomplete regularized gamma function
    /// `P(a,x) = 1 / Gamma(a) * int(exp(-t)t^(a-1), t=0..x) for real a > 0, x > 0`
    /// where `a` is the argument for the gamma function and `x` is the upper
    /// integral limit.
    ///
    /// # Remarks
    ///
    /// Returns `f64::NAN` if either argument is `f64::NAN`
    ///
    /// # Panics
    ///
    /// if `a` or `x` are not in `(0, +inf)`
    pub fn gamma_lr(a: f64, x: f64) -> f64 {
        checked_gamma_lr(a, x).unwrap()
    }

    /// Computes the lower incomplete regularized gamma function
    /// `P(a,x) = 1 / Gamma(a) * int(exp(-t)t^(a-1), t=0..x) for real a > 0, x > 0`
    /// where `a` is the argument for the gamma function and `x` is the upper
    /// integral limit.
    ///
    /// # Remarks
    ///
    /// Returns `f64::NAN` if either argument is `f64::NAN`
    ///
    /// # Errors
    ///
    /// if `a` or `x` are not in `(0, +inf)`
    pub fn checked_gamma_lr(a: f64, x: f64) -> Result<f64> {
        if a.is_nan() || x.is_nan() {
            return Ok(f64::NAN);
        }
        if a <= 0.0 || a == f64::INFINITY {
            return Err(anyhow!("AInvalid"));
        }
        if x <= 0.0 || x == f64::INFINITY {
            return Err(anyhow!("XInvalid"));
        }

        let eps = 0.000000000000001;
        let big = 4503599627370496.0;
        let big_inv = 2.22044604925031308085e-16;

        if almost_eq(a, 0.0, DEFAULT_F64_ACC) {
            return Ok(1.0);
        }
        if almost_eq(x, 0.0, DEFAULT_F64_ACC) {
            return Ok(0.0);
        }

        let ax = a * x.ln() - x - ln_gamma(a);
        if ax < -709.78271289338399 {
            if a < x {
                return Ok(1.0);
            }
            return Ok(0.0);
        }
        if x <= 1.0 || x <= a {
            let mut r2 = a;
            let mut c2 = 1.0;
            let mut ans2 = 1.0;
            loop {
                r2 += 1.0;
                c2 *= x / r2;
                ans2 += c2;

                if c2 / ans2 <= eps {
                    break;
                }
            }
            return Ok(ax.exp() * ans2 / a);
        }

        let mut y = 1.0 - a;
        let mut z = x + y + 1.0;
        let mut c = 0;

        let mut p3 = 1.0;
        let mut q3 = x;
        let mut p2 = x + 1.0;
        let mut q2 = z * x;
        let mut ans = p2 / q2;

        loop {
            y += 1.0;
            z += 2.0;
            c += 1;
            let yc = y * f64::from(c);

            let p = p2 * z - p3 * yc;
            let q = q2 * z - q3 * yc;

            p3 = p2;
            p2 = p;
            q3 = q2;
            q2 = q;

            if p.abs() > big {
                p3 *= big_inv;
                p2 *= big_inv;
                q3 *= big_inv;
                q2 *= big_inv;
            }

            if q != 0.0 {
                let nextans = p / q;
                let error = ((ans - nextans) / nextans).abs();
                ans = nextans;

                if error <= eps {
                    break;
                }
            }
        }
        Ok(1.0 - ax.exp() * ans)
    }

    /// Computes the Digamma function which is defined as the derivative of
    /// the log of the gamma function. The implementation is based on
    /// "Algorithm AS 103", Jose Bernardo, Applied Statistics, Volume 25, Number 3
    /// 1976, pages 315 - 317
    pub fn digamma(x: f64) -> f64 {
        let c = 12.0;
        let d1 = -0.57721566490153286;
        let d2 = 1.6449340668482264365;
        let s = 1e-6;
        let s3 = 1.0 / 12.0;
        let s4 = 1.0 / 120.0;
        let s5 = 1.0 / 252.0;
        let s6 = 1.0 / 240.0;
        let s7 = 1.0 / 132.0;

        if x == f64::NEG_INFINITY || x.is_nan() {
            return f64::NAN;
        }
        if x <= 0.0 && ulps_eq!(x.floor(), x) {
            return f64::NEG_INFINITY;
        }
        if x < 0.0 {
            return digamma(1.0 - x) + PI / (-PI * x).tan();
        }
        if x <= s {
            return d1 - 1.0 / x + d2 * x;
        }

        let mut result = 0.0;
        let mut z = x;
        while z < c {
            result -= 1.0 / z;
            z += 1.0;
        }

        if z >= c {
            let mut r = 1.0 / z;
            result += z.ln() - 0.5 * r;
            r *= r;

            result -= r * (s3 - r * (s4 - r * (s5 - r * (s6 - r * s7))));
        }
        result
    }

    pub fn inv_digamma(x: f64) -> f64 {
        if x.is_nan() {
            return f64::NAN;
        }
        if x == f64::NEG_INFINITY {
            return 0.0;
        }
        if x == f64::INFINITY {
            return f64::INFINITY;
        }
        let mut y = x.exp();
        let mut i = 1.0;
        while i > 1e-15 {
            y += i * signum(x - digamma(y));
            i /= 2.0;
        }
        y
    }

    // modified signum that returns 0.0 if x == 0.0. Used
    // by inv_digamma, may consider extracting into a public
    // method
    fn signum(x: f64) -> f64 {
        if x == 0.0 { 0.0 } else { x.signum() }
    }

    trait GammaFn {
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
}
