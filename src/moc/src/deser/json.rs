
use std::io::Write;
use std::error::Error;

use byteorder::WriteBytesExt;
use serde_json::{self, Value, value::Value::Array};

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::elem::cell::Cell;
use crate::elemset::cell::{Cells, MocCells};
use crate::moc::{
  HasMaxDepth, CellMOCIterator,
  cell::CellMOC
};
use crate::moc2d::{
  CellMOC2ElemIt, CellMOC2Iterator,
  cell::{CellMOC2Elem, CellMOC2}
};

/// Write a JSON following the Aladin JSON format.
pub fn to_json_aladin<T, Q, I, W>(it: I, fold: &Option<usize>, line_prefix: &str, mut writer: W)
  -> std::io::Result<()>
  where
    T: Idx,
    Q: MocQty<T>,
    I: CellMOCIterator<T, Qty=Q>,
    W: Write
{
  // Create n_depth strings (1 per depth)
  let mut s_by_depth: Vec<String> = (0..=it.depth_max()).into_iter()
    .map(|i| format!("{}  \"{}\": [", line_prefix, i))
    .collect();
  // Fill them
  for c in it {
    let buff = &mut s_by_depth[c.depth as usize];
    let s = format!("{}, ", c.idx);
    match fold {
      Some(n_chars) => {
        let l = buff.rfind('\n').unwrap_or(0);
        if buff.len() - l + s.len() > *n_chars {
          buff.push_str(&format!("\n{}    ", line_prefix));
        }
        buff.push_str(&s)
      },
      None => buff.push_str(&s),
    }
  }
  // Finally write the result
  writer.write_all(b"{\n")?;
  let mut first = true;
  for s in s_by_depth {
    if !s.ends_with('[') {
      if first {
        first = false;
        // remove the last ',' and ' ' characters
        write!(writer, "{}]", &s[..s.len() - 2])?;
      } else {
        // remove the last ',' and ' ' characters
        // We could use again rfind to ckeck wether ']'
        // must be on a newline or not.
        // Let's wait for someone to complain...
        write!(writer, ",\n{}]", &s[..s.len() - 2])?;
      }
    }
  }
  writer.write_all(format!("\n{}}}", line_prefix).as_bytes())?;
  Ok(())
}


pub fn cellmoc2d_to_json_aladin<T, Q, I, U, R, J, K, L, W>(
  moc2_it: L,
  fold: &Option<usize>,
  mut writer: W
) -> std::io::Result<()>
  where
    T: Idx,
    Q: MocQty<T>,
    I: CellMOCIterator<T, Qty=Q>,
    U: Idx,
    R: MocQty<U>,
    J: CellMOCIterator<U, Qty=R>,
    K: CellMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: CellMOC2Iterator<T, Q, I, U, R, J, K>,
    W: Write
{
  writer.write_all(b"[\n")?;
  let mut is_first = true;
  for e in moc2_it {
    let (moc1_it, moc2_it) = e.cell_mocs_it();
    if is_first {
      is_first = false;
      writer.write_all(b"{\n  \"")?;
    } else {
      writer.write_all(b",\n{\n  \"")?;
    }
    writer.write_u8(Q::PREFIX as u8)?;
    writer.write_all(b"\": ")?;
    to_json_aladin(moc1_it, &fold, "  ", &mut writer)?;

    writer.write_all(b",\n  \"")?;
    writer.write_u8(R::PREFIX as u8)?;
    writer.write_all(b"\": ")?;
    to_json_aladin(moc2_it, &fold, "  ",&mut writer)?;
    writer.write_all(b"\n}")?;
  }
  writer.write_all(b"\n]\n")?;
  Ok(())
}

/*pub fn rangemoc2d_to_json_aladin<T, Q, I, U, R, J, K, L, W>(
  moc2_it: L,
  fold: &Option<usize>,
  mut writer: W
) -> std::io::Result<()>
  where
    T: Idx,
    Q: MocQty<T>,
    I: RangeMOCIterator<T, Qty=Q>,
    U: Idx,
    R: MocQty<U>,
    J: RangeMOCIterator<U, Qty=R>,
    K: RangeMOC2ElemIt<T, Q, U, R, It1=I, It2=J>,
    L: RangeMOC2Iterator<T, Q, I, U, R, J, K>,
    W: Write
{
  writer.write_all(b"[\n")?;
  let mut is_first = true;
  for e in moc2_it {
    let (moc1_it, moc2_it) = e.range_mocs_it();
    let (moc1_it, moc2_it) = (moc1_it.cells(), moc2_it.cells());
    if is_first {
      is_first = false;
      writer.write_all(b"{\n  \"")?;
    } else {
      writer.write_all(b",\n{\n  \"")?;
    }
    writer.write_u8(Q::PREFIX as u8)?;
    writer.write_all(b"\": ")?;
    to_json_aladin(moc1_it, &fold, "  ", &mut writer)?;

    writer.write_all(b",\n  \"")?;
    writer.write_u8(R::PREFIX as u8)?;
    writer.write_all(b"\": ")?;
    to_json_aladin(moc2_it, &fold, "  ", &mut writer)?;
    writer.write_all(b"\n}")?;
  }
  writer.write_all(b"\n]\n")?;
  Ok(())
}*/


/// Read a JSON following the Aladin JSON format.
pub fn from_json_aladin<T, Q>(input: &str) -> Result<CellMOC<T, Q>, Box<dyn Error>>
  where
    T: Idx, // + From<u64>,
    Q: MocQty<T>,
{
  // If streamming is needed, see:
  //   https://docs.serde.rs/serde_json/struct.StreamDeserializer.html
  // See https://docs.serde.rs/serde_json/#operating-on-untyped-json-values
  let root: Value = serde_json::from_str(input)?;
  from_json_aladin_internal::<T, Q>(&root)
}

fn from_json_aladin_internal<T, Q>(value: &Value) -> Result<CellMOC<T, Q>, Box<dyn Error>>
  where
    T: Idx,
    Q: MocQty<T>,
{
  match value {
    Value::Object(map) => {
      // First reserve the exact number of cells, summing the len of all arrays
      let mut cells: Vec<Cell<T>> = Vec::with_capacity(
        (0..=Q::MAX_DEPTH).into_iter()
          .filter_map(|d| match map.get(&d.to_string()) {
            Some(Array(vec)) => Some(vec.len()),
            _ => None
          })
          .sum()
      );
      // Convert each array element in a cell and add it to the cell list
      let mut depth_max = 0;
      for depth in 0..=Q::MAX_DEPTH {
        if let Some(Array(vec)) = map.get(&depth.to_string()) {
          for v in vec.iter()
            .filter_map(|v| match v {
              Value::Number(idx) => idx.as_u64().and_then(|v| T::from_u64(v).into() /*v.into()*/ /*T::try_from(v).ok()*/),
              _ => None,
            }) {
            cells.push(Cell::new(depth, v.into()));
          }
          depth_max = depth_max.max(depth);
        }
      }
      // Sort the cell list
      cells.sort_by(|a, b| a.flat_cmp::<Q>(b));
      // Check for non-overlapping cells
      // TODO: add non-overlapping check
      Ok(CellMOC::new(depth_max, MocCells::new(Cells::new(cells))))
    },
    _ => Err(format!("Wrong JSON root type. Expected: Object. Actual: {:?}", &value).into()),
  }
}

/// Read a JSON following the Aladin JSON format.
pub fn cellmoc2d_from_json_aladin<T, Q, U, R>(input: &str) -> Result<CellMOC2<T, Q, U, R>, Box<dyn Error>>
  where
    T: Idx, // + From<u64>,
    Q: MocQty<T>,
    U: Idx, // + From<u64>,
    R: MocQty<U>,
{
  // If streamming is needed, see:
  //   https://docs.serde.rs/serde_json/struct.StreamDeserializer.html
  // See https://docs.serde.rs/serde_json/#operating-on-untyped-json-values
  let root: Value = serde_json::from_str(input)?;
  let mut depth_max_l = 0_u8;
  let mut depth_max_r = 0_u8;
  let mut elems: Vec<CellMOC2Elem<T, Q, U, R>> = Vec::with_capacity(100);
  match &root {
    Value::Array(entries) => {
      for entry in entries {
        if let Value::Object(map) = entry {
          let moc1 = map.get(&Q::PREFIX.to_string());
          let moc2 = map.get(&R::PREFIX.to_string());
          match (moc1, moc2) {
            (Some(obj1), Some(obj2)) => {
              let l: CellMOC<T, Q> = from_json_aladin_internal::<T, Q>(obj1)?;
              let r: CellMOC<U, R> = from_json_aladin_internal::<U, R>(obj2)?;
              depth_max_l = depth_max_l.max(l.depth_max());
              depth_max_r = depth_max_r.max(r.depth_max());
              elems.push(CellMOC2Elem::new(l, r));
            },
            _ => return Err(format!("Wrong JSON array object type. Expected: (Object, Object). Actual: ({:?}, {:?})", &moc1, &moc2).into()),
          }
        } else {
          return Err(format!("Wrong JSON array elem type. Expected: Object. Actual: {:?}", &entry).into());
        }
      }
    },
    _ => return Err(format!("Wrong JSON root type. Expected: Array. Actual: {:?}", &root).into()),
  }
  Ok(CellMOC2::new(depth_max_l, depth_max_r, elems))
}


// json_stream :
// {
//   "qty": "",
//   "depth": "",
//   "cells": [uniq, uniq, ...]
// }

#[cfg(test)]
mod tests {
  
  use std::fs;
  use std::path::PathBuf;

  use crate::qty::{Hpx, Time};
  use crate::elemset::range::{MocRanges, TimeRanges, HpxRanges};
  use crate::moc::{
    RangeMOCIterator, RangeMOCIntoIterator,
    CellMOCIterator, CellMOCIntoIterator,
    range::{RangeMOC}
  };
  use crate::moc2d::CellMOC2IntoIterator;
  use crate::moc2d::range::{RangeMOC2Elem, RangeMOC2};
  use crate::deser::json::{from_json_aladin, cellmoc2d_to_json_aladin, cellmoc2d_from_json_aladin};

  #[test]
  fn test_fromto_json() {
    let rm = RangeMOC::new(29,
      MocRanges::<u64, Hpx<u64>>::new_unchecked(vec![
        0..5,
        6..59,
        78..6953,
        12458..55587,
        55787..65587
      ])
    );
    let mut sink = Vec::new();
    (&rm).into_range_moc_iter()
      .cells()
      .to_json_aladin(Some(40), &mut sink)
      .unwrap();
     let json = String::from_utf8_lossy(&sink);
    // println!("{}\n", &json);
    let mut sink2 = Vec::new();
    from_json_aladin::<u64, Hpx<u64>>(&json)
      .unwrap()
      .into_cell_moc_iter()
      .to_json_aladin(Some(40), &mut sink2)
      .unwrap();
    // let json = String::from_utf8_lossy(&sink2);
    // println!("{}\n", &json);
    assert_eq!(sink, sink2);

    let mut sink = Vec::new();
    (&rm).into_range_moc_iter()
      .cells()
      .to_json_aladin(None, &mut sink)
      .unwrap();
    let json = String::from_utf8_lossy(&sink);
    // println!("{}\n", &json);
    let mut sink2 = Vec::new();
    from_json_aladin::<u64, Hpx<u64>>(&json)
      .unwrap()
      .into_cell_moc_iter()
      .to_json_aladin(None, &mut sink2)
      .unwrap();
    // let json = String::from_utf8_lossy(&sink2);
    // println!("{}\n", &json);
    assert_eq!(sink, sink2);
  }

  #[test]
  fn test_stmoc_tofrom_json() {
    // Build moc
    let mut elems: Vec<RangeMOC2Elem<u64, Time<u64>, u64, Hpx<u64>>> = Default::default();
    elems.push(
      RangeMOC2Elem::new(
        RangeMOC::new(61, TimeRanges::new_unchecked(vec![1..2, 3..4, 5..6])),
        RangeMOC::new(4, HpxRanges::new_unchecked(vec![4503599627370496..18014398509481984]))
      )
    );
    elems.push(
      RangeMOC2Elem::new(
        RangeMOC::new(61, TimeRanges::new_unchecked(vec![50..51, 52..53])),
        RangeMOC::new(4, HpxRanges::new_unchecked(vec![28147497671065600..29273397577908224]))
      )
    );
    let moc2 = RangeMOC2::new(61, 4, elems);
    //  
    let mut sink = Vec::new();
    // rangemoc2d_to_json_aladin(moc2.into_cell_moc2_iter(), &None, &mut sink);
    cellmoc2d_to_json_aladin(moc2.into_cell_moc2_iter(), &None, &mut sink).unwrap();
    let json = String::from_utf8_lossy(&sink);
    let _cellmoc2 = cellmoc2d_from_json_aladin::<u64, Time<u64>, u64, Hpx<u64>>(&json);
    // cellmoc2.into_cell_moc2_iter().
    println!("{}\n", &json);
  }

  #[test]
  fn test_stmoc_from_json_large() {
    let path_buf1 = PathBuf::from("resources/MOC2.0/STMOC-test4FX.json");
    let path_buf2 = PathBuf::from("../resources/MOC2.0/STMOC-test4FX.json");
    let json = fs::read_to_string(&path_buf1).or_else(|_| fs::read_to_string(&path_buf2)).unwrap();
    let cellmoc2 = cellmoc2d_from_json_aladin::<u64, Time<u64>, u64, Hpx<u64>>(&json).unwrap();
    // cellmoc2.into_cell_moc2_iter().
    // println!("{}", cellmoc2.n_entries());
    assert_eq!(cellmoc2.n_entries(), 3532);
    /*let mut sink = Vec::new();
    cellmoc2d_to_json_aladin(cellmoc2.into_cell_moc2_iter(), &None, &mut sink);
    let json = String::from_utf8_lossy(&sink);
    let mut file = File::create("STMOC_P_FX.json").unwrap();
    file.write_all(json.as_ref().as_bytes()).unwrap();
    // println!("{}\n", &json);
    */
  }
}