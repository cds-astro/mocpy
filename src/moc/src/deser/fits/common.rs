//! Contains common method to parse FITS

use std::io::{BufRead, Write};
use std::str::{self, FromStr};
use std::slice::ChunksExact;
use std::num::ParseIntError;

use quick_error::ResultExt;

use crate::deser::fits::error::FitsError;

const VALUE_INDICATOR: &[u8; 2] = b"= ";

pub(super) fn write_primary_hdu<R: Write>(
  writer: &mut R,
) -> Result<(), FitsError> {
  let mut header_block = [b' '; 2880];
     header_block[0..30].copy_from_slice(b"SIMPLE  =                    T");
   header_block[80..110].copy_from_slice(b"BITPIX  =                    8");
  header_block[160..190].copy_from_slice(b"NAXIS   =                    0");
  header_block[240..270].copy_from_slice(b"EXTEND  =                    T");
  header_block[320..323].copy_from_slice(b"END");
  writer.write_all(&header_block[..])?;
  Ok(())
}


/// # Params
/// - `dest` the
#[allow(dead_code)]
pub(super) fn write_uint_keyword_record(dest: &mut [u8], keyword: &[u8; 8], val: u64) {
  write_keyword_record(dest, keyword, &val.to_string())
}

/// # Params
/// - `dest` the
#[allow(dead_code)]
pub(super) fn write_str_keyword_record(dest: &mut [u8], keyword: &[u8; 8], val: &str) {
  write_keyword_record(dest, keyword, &format!("'{}'", val))
}

/// # Params
/// - `dest` destination, must contains 80 bytes
/// - `value_part` is string, must be already quote: `'str_value'`
pub(super) fn write_keyword_record(dest: &mut [u8], keyword: &[u8; 8], value_part: &str) {
  debug_assert_eq!(dest.len(), 80);
  dest[0..8].copy_from_slice(&keyword[..]);
  dest[8..10].copy_from_slice(VALUE_INDICATOR);
  let val_bytes = value_part.as_bytes();
  dest[10..10 + val_bytes.len()].copy_from_slice(val_bytes);
}


pub(super) fn write_uint_mandatory_keyword_record(dest: &mut [u8], keyword: &[u8; 8], val: u64) {
  debug_assert_eq!(dest.len(), 80);
  dest[0..8].copy_from_slice(&keyword[..]);
  dest[8..10].copy_from_slice(VALUE_INDICATOR);
  let val = val.to_string();
  let val_bytes = val.as_bytes();
  dest[30 - val_bytes.len()..30].copy_from_slice(val_bytes);
}


// READ PART

/// # Params
/// - `header_block`: re-usable header block used to avoid multiple allocations
pub(super) fn consume_primary_hdu<R: BufRead>(
  reader: &mut R,
  header_block: &mut [u8; 2880]
) -> Result<(), FitsError> {
  let mut chunks_of_80 = next_36_chunks_of_80_bytes(reader, header_block)?;
  // SIMPLE = 'T' => file compliant with the FITS standard
  check_keyword_and_val(chunks_of_80.next().unwrap(), b"SIMPLE ", b"T")?;
  chunks_of_80.next().unwrap(); // Do not check for BITPIX (we expect an empty header)
  // NAXIS = 0 => we only support FITS files with no data in the primary HDU
  check_keyword_and_val(chunks_of_80.next().unwrap(), b"NAXIS ", b"0")?;
  // Ignore possible additional keywords
  while !contains_end(&mut chunks_of_80) {
    // Few chances to enter here (except if someone had a lot of things to say in the header)
    chunks_of_80 = next_36_chunks_of_80_bytes(reader, header_block)?;
  }
  Ok(())
}

pub(super) fn next_36_chunks_of_80_bytes<'a, R: BufRead>(
  reader: &'a mut R,
  header_block: &'a mut [u8; 2880]
) -> Result<ChunksExact<'a, u8>, FitsError> {
  reader.read_exact(header_block)?;
  Ok(header_block.chunks_exact(80))
}

fn contains_end<'a, I: Iterator<Item=&'a [u8]>>(chunks_of_80: &'a mut I) -> bool {
  for kw_rc in chunks_of_80 {
    debug_assert_eq!(kw_rc.len(), 80);
    if &kw_rc[0..4] == b"END " {
      return true;
    }
  }
  false
}

pub(super) fn check_keyword_and_val(
  keyword_record: &[u8],
  expected_kw: &[u8],
  expected_val: &[u8]
) -> Result<(), FitsError> {
  check_expected_keyword(keyword_record, expected_kw)?;
  check_for_value_indicator(keyword_record)?;
  check_expected_value(keyword_record, expected_val)
}

pub(super) fn check_keyword_and_parse_uint_val<T>(
  keyword_record: &[u8],
  expected_kw: &[u8]
) -> Result<T, FitsError>
  where
    T: Into<u64> + FromStr<Err=ParseIntError>
{
  check_expected_keyword(keyword_record, expected_kw)?;
  check_for_value_indicator(keyword_record)?;
  parse_uint_val::<T>(keyword_record)
}

pub(super) fn check_expected_keyword(keyword_record: &[u8], expected: &[u8]) -> Result<(), FitsError> {
  debug_assert!(keyword_record.len() == 80);     // length of a FITS keyword-record
  debug_assert!(expected.len() <= 8); // length of a FITS keyword
  if &keyword_record[..expected.len()] == expected {
    Ok(())
  } else {
    // We know what we put in it, so unsafe is ok here
    let expected = String::from(unsafe { str::from_utf8_unchecked(expected) }.trim_end());
    // Here, may contains binary data
    let actual = String::from_utf8_lossy(&keyword_record[..expected.len()]).trim_end().to_string();
    // panic!("Ecpected: {}, Actual: {}", expected, String::from_utf8_lossy(&src[..]));
    Err(FitsError::UnexpectedKeyword(expected, actual))
  }
}

pub(super) fn check_for_value_indicator(keyword_record: &[u8]) -> Result<(), FitsError> {
  debug_assert!(keyword_record.len() == 80);     // length of a FITS keyword-record
  if get_value_indicator(keyword_record) == VALUE_INDICATOR {
    Ok(())
  } else {
    let keyword_record = String::from_utf8_lossy(&keyword_record).trim_end().to_string();
    Err(FitsError::ValueIndicatorNotFound(keyword_record))
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

pub(super) fn check_expected_value(keyword_record: &[u8], expected: &[u8]) -> Result<(), FitsError> {
  debug_assert!(keyword_record.len() == 80);     // length of a FITS keyword-record
  let src = get_value(keyword_record);
  let lt_src = trim_left(src);
  if lt_src.len() >= expected.len() && &lt_src[..expected.len()] == expected {
    Ok(())
  } else {
    let keyword = String::from_utf8_lossy(&src[0..8]).trim_end().to_string();
    // We know what we put in it, so unsae is ok here
    let expected = String::from(unsafe { str::from_utf8_unchecked(expected) });
    // Here, may contains binary data
    let actual = String::from_utf8_lossy(&lt_src[..expected.len()]).to_string();
    Err(FitsError::UnexpectedValue(keyword, expected, actual))
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
pub(super) fn get_str_val_no_quote(
  keyword_record: &[u8]
) -> Result<&[u8], FitsError> {
  let mut it = get_left_trimmed_value(keyword_record).split_inclusive(|c| *c == b'\'');
  if let Some([b'\'']) = it.next() {
    if let Some([subslice @ .., b'\'']) = it.next() {
      return Ok(trim_right(subslice));
    }
  }
  let keyword_record = String::from_utf8_lossy(&keyword_record).trim_end().to_string();
  Err(FitsError::StringValueNotFound(keyword_record))
}


pub(super) fn parse_uint_val<T>(
  keyword_record: &[u8]
) -> Result<T, FitsError>
  where
    T: Into<u64> + FromStr<Err=ParseIntError>
{
  let src = get_left_trimmed_value(keyword_record);
  let to = index_of_last_digit(src);
  if to == 0 {
    let keyword_record = String::from_utf8_lossy(&keyword_record).trim_end().to_string();
    Err(FitsError::UintValueNotFound(keyword_record))
  } else {
    // we go unsafe and unwrap because we already tested for regular digits
    let str_val = unsafe { str::from_utf8_unchecked(&src[..to]) };
    let res = str_val.parse::<T>().context(str_val.to_string())?;
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

#[cfg(test)]
mod tests {
  use crate::deser::fits::common::write_primary_hdu;

  #[test]
  fn test_write_primary_hdu() {
    let mut buf: Vec<u8> = Default::default();
    write_primary_hdu(&mut buf).unwrap();
    assert_eq!(buf.len(), 2880);
  }
}