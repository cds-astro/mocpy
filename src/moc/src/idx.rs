//! This module contains the trait defining the type `Idx` that can be used to represent the
//! index value of a MOC cell, associated with utility constants and methods.

use std::mem;
use std::io::{Read, Write};
use std::str::FromStr;
use std::fmt::{Debug, Display};
use std::ops::{AddAssign, BitAndAssign};

use num::{Integer, PrimInt, ToPrimitive};
use byteorder::{ByteOrder, ReadBytesExt, WriteBytesExt};

use crate::deser::fits::keywords::TForm1;

// 'static mean that Idx does not contains any reference
pub trait Idx: 'static + Integer + PrimInt + ToPrimitive
+ AddAssign + BitAndAssign
+ FromStr + From<u8> + Send + Sync + Debug + Display + Copy {
  const N_BYTES: u8 = mem::size_of::<Self>() as u8;
  const N_BITS: u8 = Self::N_BYTES << 3;
  /// Associated TFORM for the FITS serializion
  const TFORM: TForm1;
  const MSB_MASK: Self; // mask use to switch on/select the most significant bit
  fn read<R: Read, B: ByteOrder>(reader: &mut R) -> Result<Self, std::io::Error>;
  fn write<W: Write, B: ByteOrder>(self, writer: &mut W) -> Result<(), std::io::Error>;
  fn convert<T: Idx + From<Self>>(self) -> T {
    let to: T = self.into();
    to.unsigned_shl((T::N_BITS - Self::N_BITS) as u32)
  }
  /// Like all cast operation, to be use with caution!
  fn cast_to_f64(self) -> f64;
  fn from_u64(val: u64) -> Self;
  fn to_u64(self) -> u64;
  fn from_u64_idx(idx: u64) -> Self;
  fn to_u64_idx(self) -> u64;
}

impl Idx for u8 {
  const TFORM: TForm1 = TForm1::OneB;
  const MSB_MASK: u8 = 1_u8 << (Self::N_BITS - 1) as u32;
  fn read<R: Read, B: ByteOrder>(reader: &mut R) -> Result<Self, std::io::Error> {
    reader.read_u8()
  }
  fn write<W: Write, B: ByteOrder>(self, writer: &mut W) -> Result<(), std::io::Error> {
    writer.write_u8(self)
  }
  fn cast_to_f64(self) -> f64 {
    self as f64
  }
  fn from_u64(val: u64) -> Self {
    val as u8
  }
  fn to_u64(self) -> u64 {
    self as u64
  }
  fn from_u64_idx(idx: u64) -> Self {
    const SHIFT: usize = (u64::N_BITS - u8::N_BITS) as usize;
    (idx >> SHIFT) as Self
  }
  fn to_u64_idx(self) -> u64 {
    const SHIFT: usize = (u64::N_BITS - u8::N_BITS) as usize;
    (self as u64) << SHIFT
  }
}
impl Idx for u16 {
  const TFORM: TForm1 = TForm1::OneI;
  const MSB_MASK: u16 = 1_u16 << (Self::N_BITS - 1) as u32;
  fn read<R: Read, B: ByteOrder>(reader: &mut R) -> Result<Self, std::io::Error> {
    reader.read_u16::<B>()
  }
  fn write<W: Write, B: ByteOrder>(self, writer: &mut W) -> Result<(), std::io::Error> {
    writer.write_u16::<B>(self)
  }
  fn cast_to_f64(self) -> f64 {
    self as f64
  }
  fn from_u64(val: u64) -> Self {
    val as u16
  }
  fn to_u64(self) -> u64 {
    self as u64
  }
  fn from_u64_idx(idx: u64) -> Self {
    const SHIFT: usize = (u64::N_BITS - u16::N_BITS) as usize;
    (idx >> SHIFT) as Self
  }
  fn to_u64_idx(self) -> u64 {
    const SHIFT: usize = (u64::N_BITS - u16::N_BITS) as usize;
    (self as u64) << SHIFT
  }
}
impl Idx for u32 {
  const TFORM: TForm1 = TForm1::OneJ;
  const MSB_MASK: u32 = 1_u32 << (Self::N_BITS - 1) as u32;
  fn read<R: Read, B: ByteOrder>(reader: &mut R) -> Result<Self, std::io::Error> {
    reader.read_u32::<B>()
  }
  fn write<W: Write, B: ByteOrder>(self, writer: &mut W) -> Result<(), std::io::Error> {
    writer.write_u32::<B>(self)
  }
  fn cast_to_f64(self) -> f64 {
    self as f64
  }
  fn from_u64(val: u64) -> Self {
    val as u32
  }
  fn to_u64(self) -> u64 {
    self as u64
  }
  fn from_u64_idx(idx: u64) -> Self {
    const SHIFT: usize = (u64::N_BITS - u32::N_BITS) as usize;
    (idx >> SHIFT) as Self
  }
  fn to_u64_idx(self) -> u64 {
    const SHIFT: usize = (u64::N_BITS - u32::N_BITS) as usize;
    (self as u64) << SHIFT
  }
}
impl Idx for u64 {
  const TFORM: TForm1 = TForm1::OneK;
  const MSB_MASK: u64 = 1_u64 << (Self::N_BITS - 1) as u32;
  fn read<R: Read, B: ByteOrder>(reader: &mut R) -> Result<Self, std::io::Error> {
    reader.read_u64::<B>()
  }
  fn write<W: Write, B: ByteOrder>(self, writer: &mut W) -> Result<(), std::io::Error> {
    writer.write_u64::<B>(self)
  }
  fn cast_to_f64(self) -> f64 {
    self as f64
  }
  fn from_u64(val: u64) -> Self {
    val
  }
  fn to_u64(self) -> u64 {
    self
  }
  fn from_u64_idx(idx: u64) -> Self {
    idx
  }
  fn to_u64_idx(self) -> u64 {
    self
  }
}
impl Idx for u128 {
  const TFORM: TForm1 = TForm1::TwoK;
  const MSB_MASK: u128 = 1_u128 << (Self::N_BITS - 1) as u32;
  fn read<R: Read, B: ByteOrder>(reader: &mut R) -> Result<Self, std::io::Error> {
    reader.read_u128::<B>()
  }
  fn write<W: Write, B: ByteOrder>(self, writer: &mut W) -> Result<(), std::io::Error> {
    writer.write_u128::<B>(self)
  }
  fn cast_to_f64(self) -> f64 {
    self as f64
  }
  fn from_u64(val: u64) -> Self {
    val as u128
  }
  fn to_u64(self) -> u64 {
    self as u64
  }
  fn from_u64_idx(idx: u64) -> Self {
    const SHIFT: usize = (u128::N_BITS - u64::N_BITS) as usize;
    (idx as Self) << SHIFT
  }
  fn to_u64_idx(self) -> u64 {
    const SHIFT: usize = (u128::N_BITS - u64::N_BITS) as usize;
    (self >> SHIFT) as u64
  }
}

/*
impl Idx for i16 {
  const TFORM: TForm1 = TForm1::OneI;
  const MSB_MASK: i16 = 1_i16 << (Self::N_BITS - 1) as u32;
  fn read<R: Read, B: ByteOrder>(reader: &mut R) -> Result<Self, std::io::Error> {
    reader.read_i16::<B>()
  }
  fn write<W: Write, B: ByteOrder>(self, writer: &mut W) -> Result<(), std::io::Error> {
    writer.write_i16::<B>(self)
  }
}
impl Idx for i32 {
  const TFORM: TForm1 = TForm1::OneJ;
  const MSB_MASK: i32 = 1_i32 << (Self::N_BITS - 1) as u32;
  fn read<R: Read, B: ByteOrder>(reader: &mut R) -> Result<Self, std::io::Error> {
    reader.read_i32::<B>()
  }
  fn write<W: Write, B: ByteOrder>(self, writer: &mut W) -> Result<(), std::io::Error> {
    writer.write_i32::<B>(self)
  }
}
impl Idx for i64 {
  const TFORM: TForm1 = TForm1::OneK;
  const MSB_MASK: i64 = 1_i64 << (Self::N_BITS - 1) as u32;
  fn read<R: Read, B: ByteOrder>(reader: &mut R) -> Result<Self, std::io::Error> {
    reader.read_i64::<B>()
  }
  fn write<W: Write, B: ByteOrder>(self, writer: &mut W) -> Result<(), std::io::Error> {
    writer.write_i64::<B>(self)
  }
}*/

/*impl<T> Idx for T where T: 'static + Integer + PrimInt + ToPrimitive + AddAssign
                                   + FromStr + From<u8> + TryFrom<u64>
                                   + Send + Sync + Debug + Display + Copy {}*/



/*
pub struct IdxConvert<F: Idx, T: Idx + From<F>> {
  n_bits: u32,
  _from_type: PhantomData<F>,
  _to_type: PhantomData<T>,
}

impl<F: Idx, T: Idx + From<F>> IdxConvert<F, T> {
  pub fn new() -> Self {
    IdxConvert {
      n_bits: (T::N_BITS - F::N_BITS) as u32,
      _from_type: PhantomData,
      _to_type: PhantomData
    }
  }
  pub fn convert(&self, from: F) -> T {
    let to: T = from.into();
    to.unsigned_shl(self.n_bits)
  }
}*/

