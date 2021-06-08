//! Contains common method to parse FITS

use std::io::BufRead;
use std::str::{self, FromStr};
use std::slice::ChunksExact;
use std::num::ParseIntError;

use quick_error::ResultExt;

use crate::deser::fits::error::FitsError;

/// # Params
/// - `header_block`: re-usable header block used to avoid multiple alloxations
pub(super) fn consume_primary_hdu<R: BufRead>(
  context: &str,
  reader: &mut R,
  header_block: &mut [u8; 2880]
) -> Result<(), FitsError> {
  let mut chunks_of_80 = next_80_chunks_of_80_bytes(context, reader, header_block)?;
  check_keyword_and_val(context, chunks_of_80.next().unwrap(), b"SIMPLE ", b"T")?; // standard FITS
  chunks_of_80.next().unwrap(); // Do not check for BITPIX
  check_keyword_and_val(context, chunks_of_80.next().unwrap(), b"NAXIS ", b"0")?; // No data in primary HDU
  // Do not check other possible keywords
  while !contains_end(&mut chunks_of_80) {
    // Few chances to enter here (except if someone has a lot of things to say in the header)
    chunks_of_80 = next_80_chunks_of_80_bytes(context, reader, header_block)?;
  }
  Ok(())
}

pub(super) fn next_80_chunks_of_80_bytes<'a, R: BufRead>(
  context: &'a str,
  reader: &'a mut R,
  header_block: &'a mut [u8; 2880]
) -> Result<ChunksExact<'a, u8>, FitsError> {
  reader.read_exact(header_block).context(context)?;
  Ok(header_block.chunks_exact(80))
}

fn contains_end<'a, I: Iterator<Item=&'a [u8]>>(chunks_of_80: &'a mut I) -> bool {
  for kw_rc in chunks_of_80 {
    debug_assert_eq!(kw_rc.len(), 80);
    if &kw_rc[0..4] == b"END " {
      return true;
    }
  }
  return false;
}

pub(super) fn check_keyword_and_val(
  context: &str,
  keyword_record: &[u8],
  expected_kw: &[u8],
  expected_val: &[u8]
) -> Result<(), FitsError> {
  check_expected_keyword(context, keyword_record, expected_kw)?;
  check_for_value_indicator(context, keyword_record);
  check_expected_value(context, keyword_record, expected_val)
}

pub(super) fn check_keyword_and_parse_uint_val<T>(
  context: &str,
  keyword_record: &[u8],
  expected_kw: &[u8]
) -> Result<T, FitsError>
  where
    T: Into<u64> + FromStr<Err=ParseIntError>
{
  check_expected_keyword(context, keyword_record, expected_kw)?;
  check_for_value_indicator(context, keyword_record);
  parse_uint_val::<T>(context, keyword_record)
}

pub(super) fn check_expected_keyword(context: &str, keyword_record: &[u8], expected: &[u8]) -> Result<(), FitsError> {
  debug_assert!(keyword_record.len() == 80);     // length of a FITS keyword-record
  debug_assert!(expected.len() <= 8); // length of a FITS keyword
  if &keyword_record[..expected.len()] == expected {
    Ok(())
  } else {
    let context = String::from(context);
    // We know what we put in it, so unsafe is ok here
    let expected = String::from(unsafe { str::from_utf8_unchecked(expected) }.trim_end());
    // Here, may contains binary data
    let actual = String::from_utf8_lossy(&keyword_record[..expected.len()]).trim_end().to_string();
    // panic!("Ecpected: {}, Actual: {}", expected, String::from_utf8_lossy(&src[..]));
    Err(FitsError::UnexpectedKeyword(context, expected, actual))
  }
}

pub(super) fn check_for_value_indicator(context: &str, keyword_record: &[u8]) -> Result<(), FitsError> {
  debug_assert!(keyword_record.len() == 80);     // length of a FITS keyword-record
  if get_value_indicator(keyword_record) == b"= " {
    Ok(())
  } else {
    let context = String::from(context);
    let keyword_record = String::from_utf8_lossy(&keyword_record[..]).trim_end().to_string();
    Err(FitsError::ValueIndicatorNotFound(context, keyword_record))
  }
}

pub(super) fn get_keyword(keyword_record: &[u8]) -> &[u8] {
  &keyword_record[..8]
}
pub(super) fn get_value_indicator(keyword_record: &[u8]) -> &[u8] {
  &keyword_record[8..10]
}
pub(super) fn get_value(keyword_record: &[u8]) -> &[u8] {
  &keyword_record[10..]
}
pub(super) fn get_left_trimmed_value(keyword_record: &[u8]) -> &[u8] {
  trim_left(get_value(keyword_record))
}

pub(super) fn check_expected_value(context: &str, keyword_record: &[u8], expected: &[u8]) -> Result<(), FitsError> {
  debug_assert!(keyword_record.len() == 80);     // length of a FITS keyword-record
  let src = get_value(keyword_record);
  let lt_src = trim_left(src);
  if lt_src.len() >= expected.len() && &lt_src[..expected.len()] == expected {
    Ok(())
  } else {
    let context = String::from(context);
    let keyword = String::from_utf8_lossy(&src[0..8]).trim_end().to_string();
    // We know what we put in it, so unsae is ok here
    let expected = String::from(unsafe { str::from_utf8_unchecked(expected) });
    // Here, may contains binary data
    let actual = String::from_utf8_lossy(&lt_src[..expected.len()]).to_string();
    Err(FitsError::UnexpectedValue(context, keyword, expected, actual))
  }
}

pub(super) fn trim_left(src: &[u8]) -> &[u8] {
  for (i, c) in src.iter().enumerate() {
    if !c.is_ascii_whitespace() {
      return &src[i..];
    }
  }
  &[]
}

pub(super) fn trim_right(src: &[u8]) -> &[u8] {
  for (i, c) in src.iter().enumerate().rev() {
    if !c.is_ascii_whitespace() {
      return &src[..=i];
    }
  }
  &[]
}

/// We know that the expected value does not contains a simple quote.
/// A trim_end is apply so that the result do not contains leading or trailing spaces.
pub(super) fn get_str_val_no_quote<'a>(
  context: &str,
  keyword_record: &'a [u8]
) -> Result<&'a [u8], FitsError> {
  let mut it = get_left_trimmed_value(keyword_record).split_inclusive(|c| *c == b'\'');
  if let Some([b'\'']) = it.next() {
    if let Some([subslice @ .., b'\'']) = it.next() {
      return Ok(trim_right(subslice));
    }
  }
  let context = String::from(context);
  let keyword_record = String::from_utf8_lossy(&keyword_record[..]).trim_end().to_string();
  Err(FitsError::StringValueNotFound(context, keyword_record))
}


pub(super) fn parse_uint_val<T>(
  context: &str,
  keyword_record: &[u8]
) -> Result<T, FitsError>
  where
    T: Into<u64> + FromStr<Err=ParseIntError>
{
  let src = get_left_trimmed_value(keyword_record);
  let to = index_of_last_digit(src);
  if to == 0 {
    let context = String::from(context);
    let keyword_record = String::from_utf8_lossy(&keyword_record[..]).trim_end().to_string();
    Err(FitsError::UintValueNotFound(context, keyword_record))
  } else {
    // we go unsafe and unwrap because we already tested for regular digits
    let str_val = unsafe { str::from_utf8_unchecked(&src[..to]) };
    let res = str_val.parse::<T>().context(format!("'{}' for uint value '{}'", context , str_val))?;
    Ok(res)
  }
}

pub(super) fn index_of_last_digit(src: &[u8]) -> usize {
  for (i, c) in src.iter().enumerate() {
    if !c.is_ascii_digit() {
      return i;
    }
  };
  src.len()
}
