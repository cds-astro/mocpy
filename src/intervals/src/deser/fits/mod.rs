
use std::str;
use std::ops::Range;
use std::marker::PhantomData;
use std::hint::unreachable_unchecked;
use std::io::{self, Read, BufRead};
use std::str::FromStr;
use std::num::ParseIntError;
use std::collections::HashMap;

use byteorder::BigEndian;
use quick_error::ResultExt;

use crate::ranges::Idx;
use crate::mocqty::{MocQty, Hpx, Time};
use crate::deser::fits::error::FitsError;
use crate::deser::fits::common::{
  consume_primary_hdu, next_80_chunks_of_80_bytes,
  check_keyword_and_val, check_keyword_and_parse_uint_val
};
use crate::deser::fits::keywords::{
  MocKeywordsMap, MocKeywords,
  MocVers, MocDim, Ordering, FitsCard, MocOrder, TForm1, MocOrdS, MocOrdT
};
use crate::mocell::{Cell, MocCell};
use crate::moc::{CellMOC, HasMaxDepth, ZSorted, NonOverlapping, MOCProperties, RangeMOCIterator};
use crate::mocells::{MocCells, Cells};

// FITS standard
// https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf

pub mod common;
pub mod error;
pub mod keywords;

pub enum MocIdxType<R: BufRead> {
  U16(MocQtyType<u16, R>),
  U32(MocQtyType<u32, R>),
  U64(MocQtyType<u64, R>),
}

pub enum MocQtyType<T: Idx, R: BufRead> {
  Hpx(MocType<T, Hpx<T>, R>),
  Time(MocType<T, Time<T>, R>)
}

pub enum MocType<T: Idx, Q: MocQty<T>, R: BufRead> {
  Ranges(RangeMocIterFromFits<T, Q, R>),
  Cells(CellMOC<T, Q>),
  // GUniq(CellMOC<T, Q>), // TODO: replace by an iterator
  // RUniq(),
  // RMixed(),
}


/// # Parameters
/// - `context`: information printed in the error. Typically the file name the reader comes from.
fn from_fits_ivoa<R: BufRead>(context: &str, mut reader: R) -> Result<MocIdxType<R>, FitsError> {
  let mut header_block = [0_u8; 2880];
  consume_primary_hdu(context, &mut reader,  &mut header_block)?;
  // Read the extention HDU
  let mut it80 = next_80_chunks_of_80_bytes(context, &mut reader,  &mut header_block)?;
  // See Table 10 and 17 in https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf
  check_keyword_and_val(context, it80.next().unwrap(), b"XTENSION", b"'BINTABLE'")?;
  check_keyword_and_val(context, it80.next().unwrap(), b"BITPIX  ", b"8")?;
  check_keyword_and_val(context, it80.next().unwrap(), b"NAXIS  ", b"2")?;
  let n_bytes = check_keyword_and_parse_uint_val::<u8>(context, it80.next().unwrap(), b"NAXIS1  ")?;
  let n_elems = check_keyword_and_parse_uint_val::<u64>(context, it80.next().unwrap(), b"NAXIS2 ")?;
  check_keyword_and_val(context, it80.next().unwrap(), b"PCOUNT  ", b"0")?;
  check_keyword_and_val(context, it80.next().unwrap(), b"GCOUNT  ", b"1")?;
  check_keyword_and_val(context, it80.next().unwrap(), b"TFIELDS ", b"1")?;
  // nbits = |BITPIX|xGCOUNTx(PCOUNT+NAXIS1xNAXIS2x...xNAXISn)
  // In our case (bitpix = 8, GCOUNT = 1, PCOUNT = 0) => nbytes = n_cells * size_of(T)
  // let data_size n_bytes as usize * n_cells as usize; // N_BYTES ok since BITPIX = 8
  // Read MOC keywords
  let mut moc_kws = MocKeywordsMap::new();
  'hr: loop {
    while let Some(kw_record) = it80.next() {
      // Parse only MOC related keywords and ignore others
      if let Some(mkw) = MocKeywords::is_moc_kw(context, kw_record) {
        if let Some(previous_mkw) = moc_kws.insert(mkw?) {
          // A FITS keyword MUST BE uniq (I may be more relax here, taking the last one and not complaining)
          // return Err(FitsError::MultipleKeyword(context.to_string(), previous_mkw.keyword_str().to_string()))
          eprintln!("WARNING: Keyword '{}' found more than once in {}", previous_mkw.keyword_str(), context);
        }
        // else keyword added without error
      } else if &kw_record[0..4] == b"END " {
        break 'hr;
      }
    }
    // Read next 2880 bytes
    it80 = next_80_chunks_of_80_bytes(context, &mut reader,  &mut header_block)?;
  }
  // CREATE A GUNIQ => General UNIQ in which the order is the order at the maximum depth
  // (and does not depends on the depth).
  // CREATE A RUNIQ => same order as GUNIQ, but follow the multi-order map (or bmoc) numbering
  // BOTH allow for streaming (like range)!
  // => USE GUNIQ on u32 for Vizier Catalogues ;)
  // CREATE RMIXED
  // 0: runiq 1: borne inf (max depht) 2: borne sup (max depth) => very fast binary search :)
  // println!("{:?}", &moc_kws);
  match moc_kws.get(PhantomData::<MocVers>) {
    Some(MocKeywords::MOCVers(MocVers::V2_0)) => {
      match moc_kws.get(PhantomData::<MocDim>) {
        Some(MocKeywords::MOCDim(MocDim::Space)) => {
          let depth_max = match moc_kws.get(PhantomData::<MocOrdS>) {
            Some(MocKeywords::MOCOrdS(MocOrdS { depth })) => *depth,
            _ => return Err(FitsError::MissingKeyword(context.to_string(), MocOrder::keyword_string())),
          };
          match moc_kws.get(PhantomData::<Ordering>) {
            Some(MocKeywords::Ordering(Ordering::Nuniq)) => load_s_moc_nuniq(context, reader, n_bytes, n_elems, depth_max, &moc_kws),
            Some(MocKeywords::Ordering(Ordering::Range)) => load_s_moc_range(context, reader, n_bytes, n_elems, depth_max, &moc_kws),
            // ADD GUNIQ? RUNIQ? RMIXED?
            _ => Err(FitsError::MissingKeyword(context.to_string(), Ordering::keyword_string())),
          }
        },
        Some(MocKeywords::MOCDim(MocDim::Time)) => {
          let depth_max = match moc_kws.get(PhantomData::<MocOrdT>) {
            Some(MocKeywords::MOCOrdT(MocOrdT { depth })) => *depth,
            _ => return Err(FitsError::MissingKeyword(context.to_string(), MocOrder::keyword_string())),
          };
          match moc_kws.get(PhantomData::<Ordering>) {
            Some(MocKeywords::Ordering(Ordering::Nuniq)) => Err(FitsError::UncompatibleKeywordContent(
              context.to_string(), String::from("MOCDIM  = 'TIME'"), String::from("ORDERING= 'NUNIQ'"))),
            Some(MocKeywords::Ordering(Ordering::Range)) => load_t_moc_range(context, reader, n_bytes, n_elems, depth_max, &moc_kws),
            // ADD GUNIQ? RUNIQ? RMIXED?
            _ => Err(FitsError::MissingKeyword(context.to_string(), Ordering::keyword_string())),
          }
        },
        Some(MocKeywords::MOCDim(MocDim::TimeSpace)) =>  todo!(),
        _ => Err(FitsError::MissingKeyword(context.to_string(), MocDim::keyword_string())),
      }
    },
    _ => {
      // MOC v1.0 => SMOC only
      let depth_max = match moc_kws.get(PhantomData::<MocOrder>) {
        Some(MocKeywords::MOCOrder(MocOrder { depth })) => *depth,
        _ => return Err(FitsError::MissingKeyword(context.to_string(), MocOrder::keyword_string())),
      };
      match moc_kws.get(PhantomData::<Ordering>) {
        Some(MocKeywords::Ordering(Ordering::Nuniq)) => load_s_moc_nuniq(context, reader, n_bytes, n_elems, depth_max, &moc_kws),
        Some(MocKeywords::Ordering(Ordering::Range)) => load_s_moc_range(context, reader, n_bytes, n_elems, depth_max, &moc_kws),
        _ => Err(FitsError::MissingKeyword(context.to_string(), Ordering::keyword_string())),
      }
    },
  }
}

fn load_s_moc_nuniq<R: BufRead>(
  context: &str,
  mut reader: R,
  n_bytes: u8,
  n_elems: u64,
  depth_max: u8,
  moc_kws: &MocKeywordsMap
) -> Result<MocIdxType<R>, FitsError>
{
  match (moc_kws.get(PhantomData::<TForm1>), n_bytes) {
    (Some(MocKeywords::TForm1(TForm1::OneI)), u16::N_BYTES) =>
      Ok(MocIdxType::U16(MocQtyType::Hpx(MocType::Cells(
        from_fits_nuniq::<u16, R>(context, reader, depth_max, n_elems as usize)?
      )))),
    (Some(MocKeywords::TForm1(TForm1::OneJ)), u32::N_BYTES) =>
      Ok(MocIdxType::U32(MocQtyType::Hpx(MocType::Cells(
        from_fits_nuniq::<u32, R>(context, reader, depth_max, n_elems as usize)?
      )))),
    (Some(MocKeywords::TForm1(TForm1::OneK)), u64::N_BYTES) =>
         Ok(MocIdxType::U64(MocQtyType::Hpx(MocType::Cells(
        from_fits_nuniq::<u64, R>(context, reader, depth_max, n_elems as usize)?
      )))),
    (Some(MocKeywords::TForm1(tform)), nb) =>
      Err(FitsError::UncompatibleKeywordContent(context.to_string(), format!("NAXIS1  = {}", nb), tform.to_string())),
    (None, _) =>
      Err(FitsError::MissingKeyword(context.to_string(), TForm1::keyword_string())),
    (_, _) => unreachable!(), // Except if a bug in the code, we are sure to get a TForm1
  }
}

fn load_s_moc_range<R: BufRead>(
  context: &str,
  mut reader: R,
  n_bytes: u8,
  n_elems: u64,
  depth_max: u8,
  moc_kws: &MocKeywordsMap
) -> Result<MocIdxType<R>, FitsError>
{
  match (moc_kws.get(PhantomData::<TForm1>), n_bytes) {
    (Some(MocKeywords::TForm1(TForm1::OneI)), u16::N_BYTES) =>
      Ok(MocIdxType::U16(MocQtyType::Hpx(MocType::Ranges(
        from_fits_range::<u16, Hpx::<u16>, R>(context, reader, depth_max, n_elems >> 1)?
      )))),
    (Some(MocKeywords::TForm1(TForm1::OneJ)), u32::N_BYTES) =>
      Ok(MocIdxType::U32(MocQtyType::Hpx(MocType::Ranges(
        from_fits_range::<u32, Hpx::<u32>, R>(context, reader, depth_max, n_elems >> 1)?
      )))),
    (Some(MocKeywords::TForm1(TForm1::OneK)), u64::N_BYTES) =>
      Ok(MocIdxType::U64(MocQtyType::Hpx(MocType::Ranges(
        from_fits_range::<u64, Hpx::<u64>, R>(context, reader, depth_max, n_elems >> 1)?
      )))),
    (Some(MocKeywords::TForm1(tform)), nb) =>
      return Err(FitsError::UncompatibleKeywordContent(context.to_string(), format!("NAXIS1  = {}", nb), tform.to_string())),
    (None, _) =>
      return Err(FitsError::MissingKeyword(context.to_string(), TForm1::keyword_string())),
    (_, _) => unreachable!(), // Except if a bug in the code, we are sure to get a TForm1
  }
}

fn load_t_moc_range<R: BufRead>(
  context: &str,
  mut reader: R,
  n_bytes: u8,
  n_elems: u64,
  depth_max: u8,
  moc_kws: &MocKeywordsMap
) -> Result<MocIdxType<R>, FitsError>
{
  match (moc_kws.get(PhantomData::<TForm1>), n_bytes) {
    (Some(MocKeywords::TForm1(TForm1::OneI)), u16::N_BYTES) =>
      Ok(MocIdxType::U16(MocQtyType::Time(MocType::Ranges(
        from_fits_range::<u16, Time::<u16>, R>(context, reader, depth_max, n_elems >> 1)?
      )))),
    (Some(MocKeywords::TForm1(TForm1::OneJ)), u32::N_BYTES) =>
      Ok(MocIdxType::U32(MocQtyType::Time(MocType::Ranges(
        from_fits_range::<u32, Time::<u32>, R>(context, reader, depth_max, n_elems >> 1)?
      )))),
    (Some(MocKeywords::TForm1(TForm1::OneK)), u64::N_BYTES) =>
      Ok(MocIdxType::U64(MocQtyType::Time(MocType::Ranges(
        from_fits_range::<u64, Time::<u64>, R>(context, reader, depth_max, n_elems >> 1)?
      )))),
    (Some(MocKeywords::TForm1(tform)), nb) => {
      //println!("u64: {}", u64::N_BYTES);
      return Err(FitsError::UncompatibleKeywordContent(context.to_string(), format!("NAXIS1  = {}", nb), tform.to_string()));
    },
    (None, _) => return Err(FitsError::MissingKeyword(context.to_string(), TForm1::keyword_string())),
    (_, _) => unreachable!(), // Except if a bug in the code, we are sure to get a TForm1
  }
}

/// Official HEALPix Uniq numbering.
/// The file is sorted first by depth and then by cell number.
fn from_fits_nuniq<T, R>(context: &str, mut reader: R, depth_max: u8, n_elems: usize)
                           -> Result<CellMOC<T, Hpx<T>>, FitsError>
  where
    T: Idx,
    R: BufRead
{
  let mut v: Vec<Cell<T>> = Vec::with_capacity(n_elems);
  for _ in 0..n_elems {
    v.push(Cell::from_uniq_hpx(T::read::<_, BigEndian>(&mut reader).context(context)?));
  }
  v.sort_by(|a, b| a.cmp::<Hpx::<T>>(b));
  Ok(CellMOC::new(depth_max, MocCells::<T, Hpx::<T>>::new(Cells(v))))
}

/// Generic numbering using a sentinel bit.
/// The file is sorted in a way which is independent of the depth of each cell so we can return
/// an iterator.
/// TODO: replace return type by an iterator
fn from_fits_guniq<T, Q, R>(context: &str, mut reader: R, depth_max: u8, n_elems: usize)
                         -> Result<CellMOC<T, Q>, FitsError>
  where
    T: Idx,
    Q: MocQty<T>,
    R: BufRead
{
  let mut v: Vec<Cell<T>> = Vec::with_capacity(n_elems);
  for _ in 0..n_elems {
    v.push(Cell::from_uniq::<Q>(T::read::<_, BigEndian>(&mut reader).context(context)?));
  }
  // Check is_osrted (a function already exists in nighlty rust)
  // v.sort_by(|a, b| a.cmp::<Q>(b));
  Ok(CellMOC::new(depth_max, MocCells::<T, Q>::new(Cells(v))))
}

fn from_fits_range<T, Q, R>(_context: &str, mut reader: R, depth_max: u8, n_elems: u64)
                            -> Result<RangeMocIterFromFits<T, Q, R>, FitsError>
  where
    T: Idx,
    Q: MocQty<T>,
    R: BufRead
{
  Ok(RangeMocIterFromFits::new(depth_max, reader, n_elems))
}

pub struct RangeMocIterFromFits<T: Idx, Q: MocQty<T>, R: BufRead> {
  depth_max: u8,
  reader: R,
  n_elems: u64,
  _t_type: PhantomData<T>,
  _t_qty: PhantomData<Q>,
}
impl<T: Idx, Q: MocQty<T>, R: BufRead> RangeMocIterFromFits<T, Q, R> {
  fn new(depth_max: u8, reader: R, n_elems: u64) -> Self {
    Self {
      depth_max, reader, n_elems, _t_type: PhantomData, _t_qty: PhantomData
    }
  }
}
impl<T: Idx, Q: MocQty<T>, R: BufRead> HasMaxDepth for RangeMocIterFromFits<T, Q, R> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<T: Idx, Q: MocQty<T>, R: BufRead> ZSorted for RangeMocIterFromFits<T, Q, R> { }
impl<T: Idx, Q: MocQty<T>, R: BufRead> NonOverlapping for RangeMocIterFromFits<T, Q, R> { }
impl<T: Idx, Q: MocQty<T>, R: BufRead> MOCProperties for RangeMocIterFromFits<T, Q, R> { }
impl<T: Idx, Q: MocQty<T>, R: BufRead> Iterator for RangeMocIterFromFits<T, Q, R> {
  type Item = Range<T>; // Would de better to return a Result<Range, FitError>...
  fn next(&mut self) -> Option<Self::Item> {
    if self.n_elems > 0_u64 {
      let from = T::read::<_, BigEndian>(&mut self.reader);
      let to = T::read::<_, BigEndian>(&mut self.reader);
      if let (Ok(start), Ok(end)) = (from, to) {
        self.n_elems -= 1;
        Some(Range { start, end })
      } else {
        // Early stop due to read error. Better to return a Result!
        None
      }
    } else {
        None
    }
  }
  // Declaring size_hint, a 'collect' can directly allocate the right number of elements
  fn size_hint(&self) -> (usize, Option<usize>) {
    if self.n_elems > usize::MAX as u64 {
      (usize::MAX, None)
    } else {
      (self.n_elems as usize, Some(self.n_elems as usize))
    }
  }
}
impl<T: Idx, Q: MocQty<T>, R: BufRead> RangeMOCIterator<T> for RangeMocIterFromFits<T, Q, R> {
  type Qty = Q;
}

// from_fits_hpxmap

#[cfg(test)]
mod tests {
  use std::io::BufReader;

  use crate::mocqty::Hpx;
  use crate::deser::fits::{from_fits_ivoa, FitsError, MocIdxType, MocQtyType, MocType};
  use std::path::PathBuf;
  use std::fs::File;
  use crate::moc::{HasMaxDepth, CellMOC, CellMOCIntoIterator};
  use crate::mocells::{MocCells, Cells};
  use crate::mocell::Cell;
  use std::ops::Range;

  #[test]
  fn test_err() {
    let mut buff = [0_u8; 10];
    let reader = BufReader::new(&buff[..]);
    match from_fits_ivoa("toto", reader) {
      Err(e) => assert!(matches!(e, FitsError::ReadIo(..))),
      _ => assert!(false),
    }
  }

  #[test]
  fn test_read_v1_smoc_fits() {
    let mut path_buf = PathBuf::from("resources");
    path_buf.push("MOC2.0");
    path_buf.push("SMOC.fits");
    let file = File::open(&path_buf).unwrap();
    let reader = BufReader::new(file);
    match from_fits_ivoa(path_buf.to_str().unwrap(), reader) {
      Ok(MocIdxType::U64(MocQtyType::Hpx(MocType::Cells(moc1)))) => {
        assert_eq!(moc1.depth_max(), 29);
        assert_eq!(moc1.len(), 10);
        let moc2 = CellMOC::<u64, Hpx<u64>>::new(29, MocCells::new(Cells(
          vec![
            Cell::from_uniq_hpx(259_u64),
            Cell::from_uniq_hpx(266_u64),
            Cell::from_uniq_hpx(1040_u64),
            Cell::from_uniq_hpx(1041_u64),
            Cell::from_uniq_hpx(1042_u64),
            Cell::from_uniq_hpx(1046_u64),
            Cell::from_uniq_hpx(4115_u64),
            Cell::from_uniq_hpx(4116_u64),
            Cell::from_uniq_hpx(68719476958_u64),
            Cell::from_uniq_hpx(288230376275168533_u64)
          ])));
        for (c1, c2) in moc1.into_cell_moc_iter().zip(moc2.into_cell_moc_iter()) {
          assert_eq!(c1, c2);
        }
      },
      // Err(e) => println!("{}", e),
      _ => assert!(false),
    }
  }

  #[test]
  fn test_read_v2_tmoc_fits() {
    let mut path_buf = PathBuf::from("resources");
    path_buf.push("MOC2.0");
    path_buf.push("TMOC.fits");
    let file = File::open(&path_buf).unwrap();
    let reader = BufReader::new(file);
    match from_fits_ivoa(path_buf.to_str().unwrap(), reader) {
      Ok(MocIdxType::U64(MocQtyType::Time(MocType::Ranges(mut moc)))) => {
        assert_eq!(moc.depth_max(), 35);
        assert_eq!(moc.size_hint(), (1, Some(1)));
        assert_eq!(moc.next(), Some(Range { start: 1073741824, end: 2684354560}));
        assert_eq!(moc.next(), None);
      },
      // Err(e) => println!("{}", e),
      _ => assert!(false),
    }
  }
}