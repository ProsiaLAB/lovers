#[allow(unused_imports)]
#[cfg(feature = "high_precision")]
use rug::Complete;
#[allow(unused_imports)]
#[cfg(feature = "high_precision")]
use rug::Integer;

pub trait Factorial {
    fn factorial(self) -> Self;
}

impl Factorial for u32 {
    fn factorial(self) -> Self {
        (1..=self).product()
    }
}

#[cfg(feature = "high_precision")]
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

impl Combinations for u32 {
    fn combinations(self, k: Self) -> Self {
        if k > self {
            0
        } else {
            (1..=k).fold(1, |acc, val| acc * (self - val + 1) / val)
        }
    }
}

#[cfg(feature = "high_precision")]
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

impl Binomial for u32 {
    fn binomial(self, k: Self) -> Self {
        self.combinations(k)
    }
}

#[cfg(feature = "high_precision")]
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

#[cfg(feature = "high_precision")]
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
