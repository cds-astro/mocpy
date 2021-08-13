
use std::str::FromStr;
use std::io::{self, Write, BufRead, Lines};
use std::marker::PhantomData;

use quick_error::quick_error;
use byteorder::WriteBytesExt;

use nom::{IResult};
use nom::multi::many1;
use nom::branch::alt;
use nom::sequence::{preceded, tuple, terminated};
use nom::combinator::{cut, map, map_res, opt};
use nom::error::{ParseError, VerboseError, FromExternalError, convert_error};
use nom::character::complete::{digit1, multispace0};
use nom::character::complete::char;

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::moc::{
  HasMaxDepth, ZSorted, NonOverlapping, MOCProperties,
  CellOrCellRangeMOCIterator,
  cellcellrange::CellOrCellRangeMOC
};
use crate::elem::{
  cell::Cell,
  cellrange::CellRange,
  cellcellrange::CellOrCellRange
};
use crate::elemset::cellcellrange::{CellOrCellRanges, MocCellOrCellRanges};
use crate::moc2d::{
  CellOrCellRangeMOC2ElemIt, CellOrCellRangeMOC2Iterator,
  cellcellrange::{CellOrCellRangeMOC2, CellOrCellRangeMOC2Elem}
};

quick_error! {
  #[derive(Debug)]
  pub enum AsciiError {
    /// IO error
    Io(err: io::Error) {
      from()
      display("I/O error: {}", err)
    }
    ParseError(error: String) {
      display("Parse error: {}", error)
    }
    EmptyReader {
      display("Empty reader!")
    }
    NoData {
      display("No data to be read!")
    }
    QtyExpectedAtFirstLine(qty: String, line: String) {
      display("Error at first line. Expected: 'qty={}'. Actual: {}", qty, line)
    }
    DepthExpectedAtSecondLine(line: String) {
      display("Error at second line. Expected: 'depth=DEPTH'. Actual: {}", line)
    }
    WrongFirstTokenDepthExpected(token: String) {
      display("Wrong first token. Expected: depth. Actual: {}", token)
    }
    RemainingData {
      display("No all data have been parsed!")
    }
    WrongDepthType(depth: String) {
      display("Wrong depth type. Expected type: u8. Actual value: {}", depth)
    }
    ElemNotFound(elem: String, line: String) {
      display("Element '{}' not found in '{}'.", elem, line)
    }
  }
}


/// # Parameters
/// - `fold` inspired from the unix command `fold`: wrap (more or less) each input line to fit in
///   specified width (e.g. fold = Some(80)
/// - `use_range_len`: express ranges in `from+len` instead of `from(inclusive)-end(inclusive)`
///   (use `use_range_len=false` to stay compatible with the IVOA standard)
/// # Remark
/// - The full ASCII serialization is put in memory (one buffer per possible depth) due to
///   the hierarchical ASCII representation.
/// - To stay compatible with the IVOA standard, ranges are express with an inclusive upper bound
pub fn to_ascii_ivoa<T, Q, I, W>(it: I, fold: &Option<usize>, use_range_len: bool, mut writer: W)
  -> std::io::Result<()>
  where
    T: Idx,
    Q: MocQty<T>,
    I: CellOrCellRangeMOCIterator<T, Qty=Q>,
    W: Write
{
  // Create n_depth strings (1 per depth)
  let d_max = it.depth_max() as usize;
  let mut s_by_depth: Vec<String> = (0..=d_max).into_iter()
    .map(|i| format!("{}/", i))
    .collect();
  // Fill them
  for e in it {
    let (d, s) = match e {
      CellOrCellRange::Cell(c) => {
        (c.depth as usize, format!("{} ", c.idx))
      },
      CellOrCellRange::CellRange(r ) => {
        let s = if use_range_len {
          // I don't like the inclusive upper bound (=> -1) but I have to follow the standard :o/
          format!("{}+{} ", r.range.start, r.range.end - r.range.start - T::one())
        } else {
          // I don't like the inclusive upper bound (=> -1) but I have to follow the standard :o/
          format!("{}-{} ", r.range.start, r.range.end - T::one())
        };
        (r.depth as usize, s)
      },
    };
    match fold {
      Some(n_chars) => {
        // Using rfind is not very efficient. We should store the length till last '\n'
        let l = s_by_depth[d].rfind('\n').unwrap_or(0);
        if s_by_depth[d].len() - l + s.len() > *n_chars {
          s_by_depth[d].push_str("\n ");
        }
        s_by_depth[d].push_str(&s)
      },
      None => s_by_depth[d].push_str(&s),
    }
  }
  // Finally write the result
  for (d, s) in s_by_depth.into_iter().enumerate() {
    if fold.is_none() {
      if s.ends_with('/') {
        if d == d_max {
          write!(writer, "{} ", &s)?;
        }
      } else {
        writer.write_all(s.as_bytes())?;
      }
    } else if !s.ends_with('/') || d == d_max {
      writeln!(writer, "{}", &s)?;
    }
  }
  Ok(())
}


/*fn parse_val<T: Idx>(buf: &str) -> IResult<&str,  T> {
  map_res(digit1, |s: &str| s.parse::<T>())(buf)
}*/

/// # Info
/// Internally, we use the `nom` parser with the full file content un memory.
/// Because we use intermediary tokens, this function is neither made to be especially fast
/// nor to have a small memory footprint.
pub fn from_ascii_ivoa<T, Q>(input: &str) -> Result<CellOrCellRangeMOC<T, Q>, AsciiError>
  where
    T: Idx,
    Q: MocQty<T>,
{
  // The type of the value we just parsed, plus the associated value (if any)
  enum ValType<T> {
    Depth,
    RangeWithEnd(T), // store the end
    RangeWithLen(T), // store the len
  }
  fn parse_val<'a, T: Idx, E>(buf: &'a str) -> IResult<&'a str, T, E>
    where
      E: ParseError<&'a str> + FromExternalError<&'a str, <T as FromStr>::Err>
  {
    map_res(digit1, |s: &str| s.parse::<T>())(buf)
  }
  fn parse_range_end<'a, T: Idx, E>(buf: &'a str) -> IResult<&'a str, ValType<T>, E>
    where
      E: ParseError<&'a str> + FromExternalError<&'a str, <T as FromStr>::Err>
  {
    preceded(
      char('-'),
      cut(map(parse_val, ValType::RangeWithEnd)), // cut => fail if '-' but no val
    )(buf)
  }
  fn parse_range_len<'a, T: Idx, E>(buf: &'a str) -> IResult<&'a str, ValType<T>, E>
    where
      E: ParseError<&'a str> + FromExternalError<&'a str, <T as FromStr>::Err>
  {
    preceded(
      char('+'),
      cut(map(parse_val, ValType::RangeWithLen)), // cut => fail if '+' but no val
    )(buf)
  }
  fn parse_depth_delim<'a, T: Idx, E>(buf: &'a str) -> IResult<&'a str, ValType<T>, E>
    where
      E: ParseError<&'a str>
  {
    map(char('/'), |_| ValType::Depth)(buf)
  }
  fn parse_depth_delim_or_range_end_or_range_len<'a, T: Idx, E>(buf: &'a str) -> IResult<&'a str, ValType<T>, E>
    where
      E: ParseError<&'a str> + FromExternalError<&'a str, <T as FromStr>::Err>
  {
    alt((parse_depth_delim, parse_range_end, parse_range_len))(buf)
  }
  #[derive(Debug)]
  enum Token<T> {
    Depth(T), // u8
    Cell(T),
    Range{ start: T, end: T}
  }
  // We made this parser in a way that try to minimizes rollbacks:
  // we just test for '/', '-' or '+' after digits, without having to parse again
  // the same digits.
  fn parse_token<'a, T: Idx, E>(buf: &'a str) -> IResult<&'a str, Token<T>, E>
    where
      E: ParseError<&'a str> + FromExternalError<&'a str, <T as FromStr>::Err>
  {
    map(
      tuple((parse_val::<T, E>, opt(parse_depth_delim_or_range_end_or_range_len::<T, E>))),
      |(val, void_or_range)| {
        match void_or_range {
          None => Token::Cell(val),
          Some(ValType::Depth) => Token::Depth(val),
          Some(ValType::RangeWithEnd(end)) => Token::Range{ start: val, end: end + T::one() },
          Some(ValType::RangeWithLen(len)) => Token::Range{ start: val, end: val + len + T::one()},
        }
      }
    )(buf)
  }
  // Tokenize.
  // With 'nom', I don't know if it is possible to retrieve an iterator (I don't think so) to
  // avoid loading all token in memory.
  // But the IVOA ASCII format has not been done to generate large files, so it is ok.
  fn tokenizer<'a, T: Idx, E>(buf: &'a str) -> IResult<&'a str, Vec<Token<T>>, E>
    where
      E: ParseError<&'a str> + FromExternalError<&'a str, <T as FromStr>::Err>
  {
    terminated(
      many1(preceded(multispace0, parse_token::<T, E>)),
      multispace0
    )(buf)
  }

  let (remain, tokens) = tokenizer::<T, VerboseError<&str>>(input).map_err(|err| match err {
    nom::Err::Error(e) | nom::Err::Failure(e) => AsciiError::ParseError(convert_error(input, e)),
    nom::Err::Incomplete(e) => AsciiError::ParseError(format!("Missing data to be parsed: {:?}", e)),
  })?;
  if !remain.is_empty() {
    return Err(AsciiError::RemainingData);
  }
  // Read tokens to build the list (again, we store the full list in memory. We then sort it).
  // We may have create n lists (one per depth) and performed a kind of n-way merge to get an
  // iterator from those n lists.
  // I am not sure the extra complexity is worth for the ASCII serialization
  // (the memory cost is the same), especially because we have no guarantee that the cells/range
  // in the ASCII serialization are ordered.
  let mut tokens = tokens.into_iter();
  // The first token **must** be a depth
  let mut cur_depth = match tokens.next() {
    Some(Token::Depth(depth)) => depth.to_u8().ok_or_else(|| AsciiError::WrongDepthType(depth.to_string())),
    None => Ok(0),
    Some(token) => Err(AsciiError::WrongFirstTokenDepthExpected(format!("{:?}", token))),
  }?;
  let mut depth_max = cur_depth;
  let mut moc: Vec<CellOrCellRange<T>> = Vec::with_capacity(tokens.len());
  for token in tokens {
    match token {
      Token::Depth(depth) => {
        cur_depth = depth.to_u8().ok_or_else(|| AsciiError::WrongDepthType(depth.to_string()))?;
        depth_max = depth_max.max(cur_depth);
      },
      Token::Cell(icell) => moc.push(CellOrCellRange::Cell(Cell::new(cur_depth, icell))),
      Token::Range{ start, end} => moc.push(CellOrCellRange::CellRange(CellRange::new(cur_depth, start, end))),
    }
  }
  // Sort the list
  //println!("MOC before sort {:?}", &moc);
  moc.sort_by(|a, b| a.flat_cmp::<Q>(b));
  /*println!("MOC after sort {:?}", &moc);
  let v: Vec<Range<T>> = moc.iter().map(|e| MocRange::<T, Q>::from(e).0).collect();
  println!("MOC ranges {:?}", &v);*/

  // Check non-overlaping property
  // TODO: add check for non-overlaping ranges?
  // Return the result
  Ok(CellOrCellRangeMOC::new(depth_max, MocCellOrCellRanges::new(CellOrCellRanges::new(moc))))
}

/// We could have returned an iterator, but this ASCII serialization is not really made
/// for streaming, or we should return an iterator over of Results.
/// So far we prefer to parse the full file instead.
/// For a streaming version, see the FITS deserialization.
pub fn moc2d_from_ascii_ivoa<T, Q, U, R>(input: &str)
  -> Result<CellOrCellRangeMOC2<T, Q, U, R>, AsciiError>
  where
    T: Idx,
    Q: MocQty<T>,
    U: Idx,
    R: MocQty<U>,
{
  let mut depth_max_l = 0_u8;
  let mut depth_max_r = 0_u8;
  let mut elems: Vec<CellOrCellRangeMOC2Elem<T, Q, U, R>> = Vec::with_capacity(100);
  for elem in input.trim().split(Q::PREFIX) {
    if elem.is_empty() {
      continue;
    }
    if let Some((l, r)) = elem.split_once(R::PREFIX) {
      let l: CellOrCellRangeMOC<T, Q> = from_ascii_ivoa(l)?;
      let r: CellOrCellRangeMOC<U, R> = from_ascii_ivoa(r)?;
      depth_max_l = depth_max_l.max(l.depth_max());
      depth_max_r = depth_max_r.max(r.depth_max());
      elems.push(CellOrCellRangeMOC2Elem::new(l, r));
    } else {
      return Err(AsciiError::ElemNotFound(R::PREFIX.to_string(), elem.to_string()));
    }
  }
  // TODO: check that the elements are sorted and non overalpping
  Ok(CellOrCellRangeMOC2::new(depth_max_l, depth_max_r, elems))
}

pub fn moc2d_to_ascii_ivoa<T, Q, I, U, R, J, K, L, W>(
  moc2_it: L,
  fold: &Option<usize>,
  use_range_len: bool,
  mut writer: W
) -> Result<(), AsciiError>
  where
    T: Idx,
    Q: MocQty<T>,
    I: CellOrCellRangeMOCIterator<T, Qty=Q>,
    U: Idx,
    R: MocQty<U>,
    J: CellOrCellRangeMOCIterator<U, Qty=R>,
    K: CellOrCellRangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: CellOrCellRangeMOC2Iterator<T, Q, I, U, R, J, K>,
    W: Write
{
  for e in moc2_it {
    writer.write_u8(Q::PREFIX as u8)?;
    let (moc1_it, moc2_it) = e.cellcellrange_mocs_it();
    to_ascii_ivoa(moc1_it, &fold, use_range_len, &mut writer)?;
    writer.write_u8(R::PREFIX as u8)?;
    to_ascii_ivoa(moc2_it, &fold, use_range_len, &mut writer)?;
  }
  Ok(())
}


/// This serialization is less compact than the IVOA ASCII serialization
/// (because of the multiple repetition of a same depth).
/// But:
/// * It is much simple to implement.
/// * It allows for streaming (very small memory footprint)
/// * It uses range exclusive upper bound (no -1 at serialization and +1 at deserialization).
///
/// The elements **MUST** be ordered following the ZOrder curve order.
///
/// The idea is to store the MOC in a file, so I choose the newline
/// (so we can use tools like grep/sed on it, and to simplify the deserialization code).
/// One can still use `tr '\n' ' '` to replace newline by spaces or perform the `\n`->``
/// in post/pre-treatment.
///
/// This representation is compatible with "Multi Order Map".
pub fn to_ascii_stream<T, Q, I, W>(it: I, use_range_len: bool, mut writer: W)
  -> std::io::Result<()>
  where
    T: Idx,
    Q: MocQty<T>,
    I: CellOrCellRangeMOCIterator<T, Qty=Q>,
    W: Write
{
  writeln!(writer, "qty={}", Q::NAME)?;
  writeln!(writer, "depth={}", it.depth_max())?;
  if use_range_len {
    for e in it {
      match e {
        CellOrCellRange::Cell(c) => writeln!(writer, "{}/{}", c.depth, c.idx)?,
        CellOrCellRange::CellRange(r) => writeln!(writer, "{}/{}+{}", r.depth, r.range.start, r.range.end - r.range.start)?,
      }
    }
  } else {
    for e in it {
      match e {
        CellOrCellRange::Cell(c) => writeln!(writer, "{}/{}", c.depth, c.idx)?,
        CellOrCellRange::CellRange(r) => writeln!(writer, "{}/{}-{}", r.depth, r.range.start, r.range.end)?,
      }
    }
  }
  Ok(())
}


// TODO: add a beter error handling with quickerror!
pub fn from_ascii_stream<T, Q, R>(reader: R)
  -> Result<MOCFromAsciiStream<T, Q, R>, AsciiError>
  where
    T: Idx,
    Q: MocQty<T>,
    R: BufRead
{
  let mut lines = reader.lines();
  if let Some(line) = lines.next().transpose()? {
    match line.trim().split_once('=').map(|(l, r)| (l.trim(), r.trim())) {
      Some(("qty", name)) if name == Q::NAME => (),
      _ => return Err(AsciiError::QtyExpectedAtFirstLine(Q::NAME.to_string(), line)),
    }
  } else {
    return Err(AsciiError::EmptyReader);
  }
  if let Some(line) = lines.next().transpose()? {
    let depth = match line.trim().split_once('=').map(|(l, r)| (l.trim(), r.trim().parse::<u8>())) {
      Some(("depth", Ok(depth))) => Ok(depth),
      _ => Err(AsciiError::DepthExpectedAtSecondLine(line))
    }?;
    Ok(
      MOCFromAsciiStream {
        lines,
        depth_max: depth,
        _t_type: PhantomData,
        _q_type: PhantomData
      }
    )
  } else {
    Err(AsciiError::NoData)
  }
}

pub struct MOCFromAsciiStream<T: Idx, Q: MocQty<T>, R: BufRead> {
  lines: Lines<R>,
  depth_max: u8,
  _t_type: PhantomData<T>,
  _q_type: PhantomData<Q>,
}

impl<T: Idx, Q: MocQty<T>, R: BufRead> HasMaxDepth for MOCFromAsciiStream<T, Q, R> {
  fn depth_max(&self) -> u8 {
    self.depth_max
  }
}
impl<T: Idx, Q: MocQty<T>, R: BufRead> ZSorted for MOCFromAsciiStream<T, Q, R> { }
impl<T: Idx, Q: MocQty<T>, R: BufRead> NonOverlapping for MOCFromAsciiStream<T, Q, R> { }
impl<T: Idx, Q: MocQty<T>, R: BufRead> MOCProperties for MOCFromAsciiStream<T, Q, R> { }
impl<T: Idx, Q: MocQty<T>, R: BufRead> Iterator for MOCFromAsciiStream<T, Q, R> {
  // TODO: replace Iterator<CellOrCellRange<T>> by Iterator<Result<CellOrCellRange<T>>>
  //       for better error handling
  type Item = CellOrCellRange<T>;

  fn next(&mut self) -> Option<Self::Item> {
    loop {
      match self.lines.next().transpose() {
        Ok(Some(line)) => {
          let line = line.trim();
          if !line.is_empty() {
            match line.split_once('/') {
              Some((depth, cell_or_range)) => {
                // Would be more efficient to work on a &[u8]
                // to iterate while we get numbers and then check for - or +.
                // To be change if we encounter performance issues
                if let Ok(depth) = depth.parse::<u8>() {
                  if cell_or_range.contains('-') {
                    // Parse range start-end
                    let (start, end) = cell_or_range.split_once('-').unwrap();
                    if let (Ok(start), Ok(end)) = (start.parse::<T>(), end.parse::<T>()) {
                      return Some(CellOrCellRange::CellRange(CellRange::new(depth, start, end)));
                    }
                  } else if cell_or_range.contains('+') {
                    // Parse range start+len
                    let (start, len) = cell_or_range.split_once('+').unwrap();
                    if let (Ok(start), Ok(len)) = (start.parse::<T>(), len.parse::<T>()) {
                      return Some(CellOrCellRange::CellRange(CellRange::new(depth, start, start + len)));
                    }
                  } else if let Ok(idx) = cell_or_range.parse::<T>() {
                    return Some(CellOrCellRange::Cell(Cell::new(depth, idx)));
                  }
                }
                eprintln!("Error parsing ascii stream at line: {:?}", line);
              },
              _ => eprintln!("Error parsing ascii stream at line: {:?}. Separator '/' not found.", line),
            }
          }
        },
        Ok(None) => return None,
        Err(e) => eprintln!("Error reading ascii stream: {:?}", e),
      }
    }
  }
}
impl<T: Idx, Q: MocQty<T>, R: BufRead> CellOrCellRangeMOCIterator<T> for MOCFromAsciiStream<T, Q, R> {
  type Qty = Q;

  fn peek_last(&self) -> Option<&CellOrCellRange<T>> {
    None
  }
}


#[cfg(test)]
mod tests {
  use std::str;

  use crate::qty::{Hpx, Time};
  use crate::elem::cell::Cell;
  use crate::elemset::range::MocRanges;
  use crate::moc::{
    HasMaxDepth,
    RangeMOCIterator, RangeMOCIntoIterator,
    CellMOCIterator, CellOrCellRangeMOCIterator,
    CellOrCellRangeMOCIntoIterator,
    range::RangeMOC
  };
  use crate::moc2d::CellOrCellRangeMOC2IntoIterator;
  use crate::deser::ascii::{from_ascii_stream, from_ascii_ivoa, moc2d_from_ascii_ivoa, moc2d_to_ascii_ivoa};

  #[test]
  fn test_from_ascii_ivoa() {
    let smoc_ascii = "3/3 10 4/16-18 22 5/19-20 17/222 28/123456789 29/";
    let smoc = from_ascii_ivoa::<u64, Hpx::<u64>>(&smoc_ascii).unwrap();
    let mut rit = smoc.into_cellcellrange_moc_iter().ranges();
    assert_eq!(rit.depth_max(), 29);
    assert_eq!(rit.next(), Some(493827156..493827160));
    assert_eq!(rit.next(), Some(3724541952..3741319168));
    assert_eq!(rit.next(), Some(5348024557502464..5910974510923776));
    assert_eq!(rit.next(), Some(13510798882111488..21392098230009856));
    assert_eq!(rit.next(), Some(24769797950537728..25895697857380352));
    assert_eq!(rit.next(), Some(45035996273704960..49539595901075456));
    assert_eq!(rit.next(), None);

    let tmoc_ascii = "31/1 32/4 35/";
    let tmoc = from_ascii_ivoa::<u64, Time::<u64>>(&tmoc_ascii).unwrap();
    let mut cellit = tmoc.into_cellcellrange_moc_iter().ranges();
    assert_eq!(cellit.depth_max(), 35);
    assert_eq!(cellit.next(), Some(1073741824..2684354560));
    assert_eq!(cellit.next(), None);

    let tmoc_ascii = "31/1 32/4 35/";
    let tmoc = from_ascii_ivoa::<u64, Time::<u64>>(&tmoc_ascii).unwrap();
    let mut cellit = tmoc.into_cellcellrange_moc_iter().ranges().cells();
    assert_eq!(cellit.depth_max(), 35);
    assert_eq!(cellit.next(), Some(Cell::new(31, 1)));
    assert_eq!(cellit.next(), Some(Cell::new(32, 4)));
    assert_eq!(cellit.next(), None);
  }

  #[test]
  fn test_fromto_ascii_ivoa() {
    let rm = RangeMOC::new(29,
      MocRanges::<u64, Hpx<u64>>::new_unchecked(vec![
        0..5,
        6..59,
        78..6953,
        12458..55587,
        55787..65587
      ])
    );
    let mut res_ascii_1 = Vec::new();
    (&rm).into_range_moc_iter()
      .cells()
      .cellranges()
      .to_ascii_ivoa(None, false, &mut res_ascii_1)
      .unwrap();
    println!("{}\n", str::from_utf8(&res_ascii_1).unwrap());
    let mut res_ascii_2 = Vec::new();
    from_ascii_ivoa::<u64, Hpx<u64>>(str::from_utf8(&res_ascii_1).unwrap())
      .unwrap()
      .into_cellcellrange_moc_iter()
      .to_ascii_ivoa(None, false, &mut res_ascii_2)
      .unwrap();
    assert_eq!(res_ascii_1, res_ascii_2);

    let mut res_ascii_1 = Vec::new();
    (&rm).into_range_moc_iter()
      .cells()
      .cellranges()
      .to_ascii_ivoa(None,true, &mut res_ascii_1)
      .unwrap();
    // println!("{}\n", str::from_utf8(&res_ascii_1).unwrap());
    let mut res_ascii_2 = Vec::new();
    from_ascii_ivoa::<u64, Hpx<u64>>(str::from_utf8(&res_ascii_1).unwrap())
      .unwrap()
      .into_cellcellrange_moc_iter()
      .to_ascii_ivoa(None, true, &mut res_ascii_2)
      .unwrap();
    assert_eq!(res_ascii_1, res_ascii_2);

    let mut res_ascii_1 = Vec::new();
    (&rm).into_range_moc_iter()
      .cells()
      .cellranges()
      .to_ascii_ivoa(Some(30),false, &mut res_ascii_1)
      .unwrap();
    // println!("{}\n", str::from_utf8(&res_ascii_1).unwrap());
    let mut res_ascii_2 = Vec::new();
    from_ascii_ivoa::<u64, Hpx<u64>>(str::from_utf8(&res_ascii_1).unwrap())
      .unwrap()
      .into_cellcellrange_moc_iter()
      .to_ascii_ivoa(Some(30), false, &mut res_ascii_2)
      .unwrap();
    assert_eq!(res_ascii_1, res_ascii_2);

    let mut res_ascii_1 = Vec::new();
    (&rm).into_range_moc_iter()
      .cells()
      .cellranges()
      .to_ascii_ivoa(Some(30),true, &mut res_ascii_1)
      .unwrap();
    // println!("{}\n", str::from_utf8(&res_ascii_1).unwrap());
    let mut res_ascii_2 = Vec::new();
    from_ascii_ivoa::<u64, Hpx<u64>>(str::from_utf8(&res_ascii_1).unwrap())
      .unwrap()
      .into_cellcellrange_moc_iter()
      .to_ascii_ivoa(Some(30), true, &mut res_ascii_2)
      .unwrap();
    assert_eq!(res_ascii_1, res_ascii_2);
  }

  #[test]
  fn test_fromto_ascii_stream() {
    let rm = RangeMOC::new(29,
      MocRanges::<u64, Hpx<u64>>::new_unchecked(vec![
        0..5,
        6..59,
        78..6953,
        12458..55587,
        55787..65587
      ])
    );
    let mut res_ascii_1 = Vec::new();
    (&rm).into_range_moc_iter()
      .cells()
      .cellranges()
      .to_ascii_stream(false, &mut res_ascii_1)
      .unwrap();
    let mut res_ascii_2 = Vec::new();
    from_ascii_stream::<u64, Hpx<u64>, _>(&res_ascii_1[..])
      .unwrap()
      .to_ascii_stream(false, &mut res_ascii_2)
      .unwrap();
    assert_eq!(res_ascii_1, res_ascii_2);

    let mut res_ascii_1 = Vec::new();
    rm.into_range_moc_iter()
      .cells()
      .cellranges()
      .to_ascii_stream(true, &mut res_ascii_1)
      .unwrap();
    let mut res_ascii_2 = Vec::new();
    from_ascii_stream::<u64, Hpx<u64>, _>(&res_ascii_1[..])
      .unwrap()
      .to_ascii_stream(true, &mut res_ascii_2)
      .unwrap();
    assert_eq!(res_ascii_1, res_ascii_2);
  }

  #[test]
  fn test_moc2d_fromto_ascii_ivoa() {
    let input = "t61/1 3 5 s3/1-3 t61/50 52 s4/25";
    let stmoc = moc2d_from_ascii_ivoa::<u64, Time<u64>, u64, Hpx<u64>>(input).unwrap();
    let mut res_ascii_1 = Vec::new();
    moc2d_to_ascii_ivoa((&stmoc).into_cellcellrange_moc2_iter(), &Some(20), false, &mut res_ascii_1).unwrap();
    println!("{}\n", str::from_utf8(&res_ascii_1).unwrap());

    let mut res_ascii_1 = Vec::new();
    moc2d_to_ascii_ivoa((&stmoc).into_cellcellrange_moc2_iter(), &None, true, &mut res_ascii_1).unwrap();
    println!("{}\n", str::from_utf8(&res_ascii_1).unwrap());
  }
}
