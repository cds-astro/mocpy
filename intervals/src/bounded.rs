use num::{Integer, PrimInt, One};

#[inline(always)]
const fn num_bits<T>() -> usize { 
    std::mem::size_of::<T>() * 8
}

pub trait Bounded<T>
where T: Integer + PrimInt {
    const MAXDEPTH: i8;

    #[inline(always)]
    fn log_2(x: T) -> u32 {
        num_bits::<T>() as u32 - x.leading_zeros() - 1
    }

    #[inline(always)]
    fn pix_depth(u: T) -> (u32, T) {
        let msb = Self::log_2(u) & !0x1;
        
        let depth = (msb >> 1) - 1;
        let t: T = One::one();
        let pix = u - t.unsigned_shl(msb);

        (depth, pix)
    }
}

impl Bounded<u128> for u128 {
    const MAXDEPTH: i8 = 62;
}
impl Bounded<u64> for u64 {
    const MAXDEPTH: i8 = 29;
}
impl Bounded<u32> for u32 {
    const MAXDEPTH: i8 = 13;
}
impl Bounded<u8> for u8 {
    const MAXDEPTH: i8 = 2;
}