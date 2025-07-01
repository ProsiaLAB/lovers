use std::f64::consts::LOG2_10;

use serde::Deserialize;

#[derive(Debug, Clone, Copy, Deserialize)]
pub enum Precision {
    Bits(u32),
    Digits(u32),
}

impl Precision {
    pub fn to_bits(self) -> u32 {
        match self {
            Precision::Bits(bits) => bits,
            Precision::Digits(digits) => (digits as f64 * LOG2_10).ceil() as u32,
        }
    }

    pub fn to_digits(self) -> u32 {
        match self {
            Precision::Digits(digits) => digits,
            Precision::Bits(bits) => (bits as f64 / LOG2_10).floor() as u32,
        }
    }
}
