//! This module contains the MOC specific FITS keywords

use std::str;
use std::fmt;
use std::marker::PhantomData;
use std::slice::ChunksMut;

use crate::deser::fits::common::{get_str_val_no_quote, get_keyword, parse_uint_val, write_keyword_record};
use crate::deser::fits::error::FitsError;
use std::str::FromStr;


pub trait FitsCard: Sized {
  const KEYWORD: &'static [u8; 8];

  fn keyword_str() -> &'static str {
    unsafe { str::from_utf8_unchecked(Self::KEYWORD) }
  }

  fn keyword_string() -> String {
    unsafe { String::from_utf8_unchecked(Self::KEYWORD.to_vec()) }
  }

  fn parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    Self::specific_parse_value(keyword_record)
  }

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError>;

  fn write_keyword_record(&self, keyword_record: &mut [u8]) -> Result<(), FitsError> {
    write_keyword_record(keyword_record, Self::KEYWORD, &self.to_fits_value());
    Ok(())
  }

  /// Must be in quotes `'val'` is value type is string
  fn to_fits_value(&self) -> String;

  /// Generate an error in case the parsed value does not match a pre-define list of possible values
  /// To be called in `specific_parse_value`.
  /// Essentially, it converts &str in String (because once the error is raised, the str in the
  /// read buffer are out-of-scope.
  fn predefine_val_err(
    parsed_value: &[u8],
    expected_values: &[&[u8]],
  ) -> FitsError {
    FitsError::UnexpectedValue(
      Self::keyword_string(),
      format!("{:?}", expected_values.iter()
        .map(|v| unsafe { String::from_utf8_unchecked(v.to_vec()) })
        .collect::<Vec<String>>()
      ),
      String::from_utf8_lossy(parsed_value).to_string()
    )
  }
}

#[derive(Debug)]
pub enum MocVers {
  V1_1,
  V2_0
}
impl FitsCard for MocVers {
  const KEYWORD: &'static [u8; 8] = b"MOCVERS ";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    match get_str_val_no_quote(keyword_record)? {
      b"1.1" => Ok(MocVers::V1_1),
      b"2.0" => Ok(MocVers::V2_0),
      parsed_val => Err(Self::predefine_val_err(parsed_val, &[b"1.1", b"2.0"])),
    }
  }

  fn to_fits_value(&self) -> String {
    String::from(match self {
      MocVers::V1_1 => "'1.1'",
      MocVers::V2_0 => "'2.0'",
    })
  }
}

#[derive(Debug)]
pub enum MocDim {
  Space,
  Time,
  TimeSpace,
}
impl FitsCard for MocDim {
  const KEYWORD: &'static [u8; 8] = b"MOCDIM  ";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    match get_str_val_no_quote(keyword_record)? {
      b"TIME" => Ok(MocDim::Time),
      b"SPACE" => Ok(MocDim::Space),
      b"TIME.SPACE" => Ok(MocDim::TimeSpace),
      parsed_val => Err(Self::predefine_val_err(parsed_val, &[b"TIME", b"SPACE", b"TIME.SPACE"])),
    }
  }

  fn to_fits_value(&self) -> String {
    String::from(match self {
      MocDim::Time => "'TIME'",
      MocDim::Space => "'SPACE'",
      MocDim::TimeSpace => "'TIME.SPACE'",
    })
  }
}

#[derive(Debug)]
pub enum Ordering {
  Nuniq,   // v1.2 v2.0
  Range,   //      v2.0
  Range29, //  pre v2.0
}
impl FitsCard for Ordering {
  const KEYWORD: &'static [u8; 8] = b"ORDERING";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    match get_str_val_no_quote(keyword_record)? {
      b"NUNIQ" => Ok(Ordering::Nuniq),
      b"RANGE" => Ok(Ordering::Range),
      b"RANGE29" => Ok(Ordering::Range29),
      parsed_val => Err(Self::predefine_val_err(parsed_val, &[b"NUNIQ", b"RANGE"])),
    }
  }

  fn to_fits_value(&self) -> String {
    String::from(match self {
      Ordering::Nuniq => "'NUNIQ'",
      Ordering::Range => "'RANGE'",
      Ordering::Range29 => "'RANGE29'",
    })
  }
}

#[derive(Debug)]
pub enum CoordSys {
  ICRS, // C
}
impl FitsCard for CoordSys {
  const KEYWORD: &'static [u8; 8] = b"COORDSYS";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    match get_str_val_no_quote(keyword_record)? {
      b"C" => Ok(CoordSys::ICRS),
      parsed_val => Err(Self::predefine_val_err(parsed_val, &[b"C"])),
    }
  }

  fn to_fits_value(&self) -> String {
    String::from("'C'")
  }
}

#[derive(Debug)]
pub enum TimeSys {
  TCB, // TCB
  JD, // pre V2.0
}
impl FitsCard for TimeSys {
  const KEYWORD: &'static [u8; 8] = b"TIMESYS ";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    match get_str_val_no_quote(keyword_record)? {
      b"TCB" => Ok(TimeSys::TCB),
      b"JD" => Ok(TimeSys::JD),
      parsed_val => Err(Self::predefine_val_err(parsed_val, &[b"TCB", b"JD"])),
    }
  }

  fn to_fits_value(&self) -> String {
    match self {
      TimeSys::TCB => String::from("'TCB'"),
      TimeSys::JD => String::from("'JD'"),
    }
  }
}


#[derive(Debug)]
pub enum MocType {
  Image,
  Catalog,
}
impl FitsCard for MocType {
  const KEYWORD: &'static [u8; 8] = b"MOCTYPE ";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    match get_str_val_no_quote(keyword_record)? {
      b"IMAGE" => Ok(MocType::Image),
      b"CATALOG" => Ok(MocType::Catalog),
      parsed_val => Err(Self::predefine_val_err(parsed_val, &[b"IMAGE", b"CATALOG"])),
    }
  }

  fn to_fits_value(&self) -> String {
    String::from(match self {
      MocType::Image => "'IMAGE'",
      MocType::Catalog => "'CATALOG'",
    })
  }
}
impl FromStr for MocType {
  type Err = FitsError;

  fn from_str(s: &str) -> Result<Self, Self::Err> {
    MocType::specific_parse_value(s.as_bytes())
  }
}
#[derive(Debug)]
pub enum PixType {
  Healpix,
}
impl FitsCard for PixType {
  const KEYWORD: &'static [u8; 8] = b"PIXTYPE ";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    match get_str_val_no_quote(keyword_record)? {
      b"HEALPIX" => Ok(PixType::Healpix),
      parsed_val => Err(Self::predefine_val_err(parsed_val, &[b"TCB"])),
    }
  }

  fn to_fits_value(&self) -> String {
    String::from("'HEALPIX'")
  }
}

#[derive(Debug)]
pub struct MocId {
  pub id: String,
}
impl FitsCard for MocId {
  const KEYWORD: &'static [u8; 8] = b"MOCID   ";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    get_str_val_no_quote(keyword_record)
      .map(|s| MocId { id: String::from_utf8_lossy(s).to_string() })
  }

  fn to_fits_value(&self) -> String {
    format!("'{}'", &self.id)
  }
}

#[derive(Debug)]
pub struct MocTool {
  pub tool: String,
}
impl FitsCard for MocTool {
  const KEYWORD: &'static [u8; 8] = b"MOCTOOL ";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    get_str_val_no_quote(keyword_record)
      .map(|s| MocTool { tool: String::from_utf8_lossy(s).to_string() })
  }

  fn to_fits_value(&self) -> String {
    format!("'{}'", &self.tool)
  }
}

#[derive(Debug)]
pub struct MocOrder {
  pub depth: u8,
}
impl FitsCard for MocOrder {
  const KEYWORD: &'static [u8; 8] = b"MOCORDER";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    parse_uint_val::<u8>(keyword_record).map(|depth| MocOrder { depth })
  }

  fn to_fits_value(&self) -> String {
    format!("{}", &self.depth)
  }
}

#[derive(Debug)]
pub struct MocOrdS {
  pub depth: u8,
}
impl FitsCard for MocOrdS {
  const KEYWORD: &'static [u8; 8] = b"MOCORD_S";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    parse_uint_val::<u8>(keyword_record).map(|depth| MocOrdS { depth })
  }

  fn to_fits_value(&self) -> String {
    format!("{}", &self.depth)
  }
}

#[derive(Debug)]
pub struct MocOrdT {
  pub depth: u8,
}
impl FitsCard for MocOrdT {
  const KEYWORD: &'static [u8; 8] = b"MOCORD_T";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    parse_uint_val::<u8>(keyword_record).map(|depth| MocOrdT { depth })
  }

  fn to_fits_value(&self) -> String {
    format!("{}", &self.depth)
  }
}

#[derive(Debug)]
pub enum TForm1 {
  OneB, // for u8
  OneI, // for i/u16
  OneJ, // for i/u32
  OneK, // for i/u64
  TwoK, // for i/u128 (invented!)
}
impl fmt::Display for TForm1 {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    write!(f, "{} = {}", Self::keyword_str(), self.to_fits_value())
  }
}
impl FitsCard for TForm1 {
  const KEYWORD: &'static [u8; 8] = b"TFORM1  ";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    match get_str_val_no_quote(keyword_record)? {
      b"1B" => Ok(TForm1::OneB),
      b"1I" => Ok(TForm1::OneI),
      b"1J" => Ok(TForm1::OneJ),
      b"1K" => Ok(TForm1::OneK),
      b"2K" => Ok(TForm1::TwoK),
      parsed_val => Err(Self::predefine_val_err(parsed_val, &[b"1I", b"1J", b"1K"])),
    }
  }

  fn to_fits_value(&self) -> String {
    String::from(
      match self {
        TForm1::OneB => "'1B'",
        TForm1::OneI => "'1I'",
        TForm1::OneJ => "'1J'",
        TForm1::OneK => "'1K'",
        TForm1::TwoK => "'2K'",
      }
    )
  }
}

#[derive(Debug)]
pub struct TType1 {
  pub ttype: String,
}
impl FitsCard for TType1 {
  const KEYWORD: &'static [u8; 8] = b"TTYPE1  ";

  fn specific_parse_value(keyword_record: &[u8]) -> Result<Self, FitsError> {
    get_str_val_no_quote(keyword_record)
      .map(|s| TType1 { ttype: String::from_utf8_lossy(s).to_string() })
  }

  fn to_fits_value(&self) -> String {
    format!("'{}'", &self.ttype)
  }
}

// Usse the index in an array of Option(MocKeywords) for fast retrieving of the Card :)
pub trait MocCard: FitsCard { const INDEX: u8; }
impl MocCard for MocVers  { const INDEX: u8 = 0; }
impl MocCard for MocDim   { const INDEX: u8 = 1; }
impl MocCard for Ordering { const INDEX: u8 = 2; }
impl MocCard for CoordSys { const INDEX: u8 = 3; }
impl MocCard for TimeSys  { const INDEX: u8 = 4; }
impl MocCard for MocId    { const INDEX: u8 = 5; }
impl MocCard for MocTool  { const INDEX: u8 = 6; }
impl MocCard for MocType  { const INDEX: u8 = 7; }
impl MocCard for MocOrdS  { const INDEX: u8 = 8; }
impl MocCard for MocOrdT  { const INDEX: u8 = 9; }
impl MocCard for MocOrder { const INDEX: u8 = 10; }
impl MocCard for PixType  { const INDEX: u8 = 11; }
impl MocCard for TForm1  { const INDEX: u8 = 12; }
impl MocCard for TType1  { const INDEX: u8 = 13; }

#[derive(Debug)]
pub(super) struct MocKeywordsMap {
  entries: [Option<MocKeywords>; 14],
}
impl MocKeywordsMap {
  pub(super) fn new() -> MocKeywordsMap {
    Self {
      entries: [None, None, None, None, None, None, None, None, None, None, None, None, None, None]
    }
  }

  pub(super) fn insert(&mut self, entry: MocKeywords) -> Option<MocKeywords> {
    self.entries[entry.index()].replace(entry)
  }

  pub(super) fn get<T: MocCard>(&self, _phantom: PhantomData<T>) -> Option<&MocKeywords> {
    self.entries[T::INDEX as usize].as_ref()
  }

  pub(super) fn write_all(&self, keyword_records: &mut ChunksMut<u8>) -> Result<(), FitsError> {
    for kw in self.entries.iter().filter_map(|v| v.as_ref()) {
      kw.write_keyword_record(keyword_records.next().unwrap())?;
    }
    Ok(())
  }
}

#[derive(Debug)]
pub enum MocKeywords {
  MOCVers(MocVers),   //      v2.0
  MOCDim(MocDim),     //      v2.0
  Ordering(Ordering), // v1.1 v2.0
  CoordSys(CoordSys), // v1.1 v2.0 if MOCDIM = SPACE
  TimeSys(TimeSys),   //      v2.0 if MOCDIM = TIME
  MOCId(MocId),       // v1.1 v2.0, opt
  MOCTool(MocTool),   // v1.1 v2.0, opt
  MOCType(MocType),   // v1.1 v2.0, opt
  MOCOrdS(MocOrdS),   //      v2.0
  MOCOrdT(MocOrdT),   //      v2.0
  MOCOrder(MocOrder), // v1.1
  PixType(PixType),   // v1.1
  // BINTABLE specific
  TForm1(TForm1),     // bintable
  TType1(TType1),     // bintable
}
impl MocKeywords  {

  pub(super) fn is_moc_kw(keyword_record: &[u8])
    -> Option<Result<Self, FitsError>>
  {
    // I have not yet found how to match on the FitsCard::KEYWORD associated constant :o/
    match get_keyword(keyword_record) {
      b"MOCVERS " => Some(MocVers::parse_value(keyword_record).map(MocKeywords::MOCVers)),
      b"MOCDIM  " => Some(MocDim::parse_value(keyword_record).map(MocKeywords::MOCDim)),
      b"ORDERING" => Some(Ordering::parse_value(keyword_record).map(MocKeywords::Ordering)),
      b"COORDSYS" => Some(CoordSys::parse_value(keyword_record).map(MocKeywords::CoordSys)),
      b"TIMESYS " => Some(TimeSys::parse_value(keyword_record).map(MocKeywords::TimeSys)),
      b"MOCID   " => Some(MocId::parse_value(keyword_record).map(MocKeywords::MOCId)),
      b"MOCTOOL " => Some(MocTool::parse_value(keyword_record).map(MocKeywords::MOCTool)),
      b"MOCTYPE " => Some(MocType::parse_value(keyword_record).map(MocKeywords::MOCType)),
      b"MOCORD_S" => Some(MocOrdS::parse_value(keyword_record).map(MocKeywords::MOCOrdS)),
      b"MOCORD_1" => Some(MocOrdS::parse_value(keyword_record).map(MocKeywords::MOCOrdS)), // To support pre v2 ST-MOC
      b"TORDER  " => Some(MocOrdT::parse_value(keyword_record).map(MocKeywords::MOCOrdT)), // To support pre v2 ST-MOC
      b"MOCORD_T" => Some(MocOrdT::parse_value(keyword_record).map(MocKeywords::MOCOrdT)),
      b"MOCORDER" => Some(MocOrder::parse_value(keyword_record).map(MocKeywords::MOCOrder)),
      b"PIXTYPE " => Some(PixType::parse_value(keyword_record).map(MocKeywords::PixType)),
      // BINTABLE
      b"TFORM1  " => Some(TForm1::parse_value(keyword_record).map(MocKeywords::TForm1)),
      b"TTYPE1  " => Some(TType1::parse_value(keyword_record).map(MocKeywords::TType1)),
      _ => None,
    }
  }

  fn index(&self) -> usize {
    (match self {
      MocKeywords::MOCVers(_) => MocVers::INDEX,
      MocKeywords::MOCDim(_) => MocDim::INDEX,
      MocKeywords::Ordering(_) => Ordering::INDEX,
      MocKeywords::CoordSys(_) => CoordSys::INDEX,
      MocKeywords::TimeSys(_) => TimeSys::INDEX,
      MocKeywords::MOCId(_) => MocId::INDEX,
      MocKeywords::MOCTool(_) => MocTool::INDEX,
      MocKeywords::MOCType(_) => MocType::INDEX,
      MocKeywords::MOCOrdS(_) => MocOrdS::INDEX,
      MocKeywords::MOCOrdT(_) => MocOrdT::INDEX,
      MocKeywords::MOCOrder(_) => MocOrder::INDEX,
      MocKeywords::PixType(_) => PixType::INDEX,
      // BINTABLE
      MocKeywords::TForm1(_) => TForm1::INDEX,
      MocKeywords::TType1(_) => PixType::INDEX,
    }) as usize
  }

  pub(super) fn keyword(&self) -> &'static [u8; 8] {
    match self {
      MocKeywords::MOCVers(_) => MocVers::KEYWORD,
      MocKeywords::MOCDim(_) => MocDim::KEYWORD,
      MocKeywords::Ordering(_) => Ordering::KEYWORD,
      MocKeywords::CoordSys(_) => CoordSys::KEYWORD,
      MocKeywords::TimeSys(_) => TimeSys::KEYWORD,
      MocKeywords::MOCId(_) => MocId::KEYWORD,
      MocKeywords::MOCTool(_) => MocTool::KEYWORD,
      MocKeywords::MOCType(_) => MocType::KEYWORD,
      MocKeywords::MOCOrdS(_) => MocOrdS::KEYWORD,
      MocKeywords::MOCOrdT(_) => MocOrdT::KEYWORD,
      MocKeywords::MOCOrder(_) => MocOrder::KEYWORD,
      MocKeywords::PixType(_) => PixType::KEYWORD,
      // BINTABLE
      MocKeywords::TForm1(_) => TForm1::KEYWORD,
      MocKeywords::TType1(_) => TType1::KEYWORD,
    }
  }

  pub(super) fn keyword_str(&self) -> &str {
    unsafe{ str::from_utf8_unchecked(self.keyword()) }.trim_end()
  }

  fn write_keyword_record(&self, keyword_record: &mut [u8]) -> Result<(), FitsError> {
    match self {
      MocKeywords::MOCVers(kw) => kw.write_keyword_record(keyword_record),
      MocKeywords::MOCDim(kw) => kw.write_keyword_record(keyword_record),
      MocKeywords::Ordering(kw) => kw.write_keyword_record(keyword_record),
      MocKeywords::CoordSys(kw) => kw.write_keyword_record(keyword_record),
      MocKeywords::TimeSys(kw) => kw.write_keyword_record(keyword_record),
      MocKeywords::MOCId(kw) => kw.write_keyword_record(keyword_record),
      MocKeywords::MOCTool(kw) => kw.write_keyword_record(keyword_record),
      MocKeywords::MOCType(kw) => kw.write_keyword_record(keyword_record),
      MocKeywords::MOCOrdS(kw) => kw.write_keyword_record(keyword_record),
      MocKeywords::MOCOrdT(kw) => kw.write_keyword_record(keyword_record),
      MocKeywords::MOCOrder(kw) => kw.write_keyword_record(keyword_record),
      MocKeywords::PixType(kw) => kw.write_keyword_record(keyword_record),
      // BINTABLE
      MocKeywords::TForm1(kw) => kw.write_keyword_record(keyword_record),
      MocKeywords::TType1(kw) => kw.write_keyword_record(keyword_record),
    }
  }
}
