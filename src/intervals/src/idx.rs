//! This module contains the trait defining the type `Idx` that can be used to represent the
//! index value of a MOC cell, associated with utility constants and methods.

use std::mem;
use std::io::{Read, Write};
use std::str::FromStr;
use std::fmt::{Debug, Display};
use std::ops::AddAssign;
use std::convert::TryFrom;

use num::{Integer, PrimInt, ToPrimitive};
use byteorder::{ByteOrder, ReadBytesExt, WriteBytesExt};

use crate::deser::fits::keywords::TForm1;

// 'static mean that Idx does not contains any reference
pub trait Idx: 'static + Integer + PrimInt + ToPrimitive + AddAssign
+ FromStr + From<u8> + TryFrom<u64>
+ Send + Sync + Debug + Display + Copy {
  const N_BYTES: u8 = mem::size_of::<Self>() as u8;
  const N_BITS: u8 = Self::N_BYTES << 3;
  /// Associated TFORM for the FITS serializion
  const TFORM: TForm1;
  const MSB_MASK: Self; // mask use to switch on/select the most significant bit
  fn read<R: Read, B: ByteOrder>(reader: &mut R) -> Result<Self, std::io::Error>;
  fn write<W: Write, B: ByteOrder>(self, writer: &mut W) -> Result<(), std::io::Error>;
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
}
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
}

/*impl<T> Idx for T where T: 'static + Integer + PrimInt + ToPrimitive + AddAssign
                                   + FromStr + From<u8> + TryFrom<u64>
                                   + Send + Sync + Debug + Display + Copy {}*/
