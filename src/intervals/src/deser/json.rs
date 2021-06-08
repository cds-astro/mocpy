
use std::io::Write;

use serde_json;

use crate::ranges::Idx;
use crate::mocqty::MocQty;
use crate::moc::{CellMOCIterator, CellMOC};
use serde_json::Value;
use std::error::Error;
use crate::mocell::Cell;
use serde_json::value::Value::Array;
use crate::mocells::{MocCells, Cells};

/// Write a JSON following the Aladin JSON format.
pub fn to_json_aladin<T, Q, I, W>(it: I, fold: Option<usize>, mut writer: W)
  -> std::io::Result<()>
  where
    T: Idx,
    Q: MocQty<T>,
    I: CellMOCIterator<T, Qty=Q>,
    W: Write
{
  // Create n_depth strings (1 per depth)
  let mut s_by_depth: Vec<String> = (0..=it.depth_max()).into_iter()
    .map(|i| format!("  \"{}\": [", i))
    .collect();
  // Fill them
  for c in it {
    let buff = &mut s_by_depth[c.depth as usize];
    let s = format!("{}, ", c.idx);
    match fold {
      Some(n_chars) => {
        let l = buff.rfind('\n').unwrap_or(0);
        if buff.len() - l + s.len() > n_chars {
          buff.push_str("\n    ");
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
        // remove the las ',' and ' ' characters
        write!(writer, "{}]", &s[..s.len() - 2])?;
      } else {
        // remove the las ',' and ' ' characters
        // We could use again rfind to ckeck wether ']'
        // must be on a newline or not.
        // Let's wait for someone to complain...
        write!(writer, ",\n{}]", &s[..s.len() - 2])?;
      }
    }
  }
  writer.write_all(b"\n}\n")?;
  Ok(())
}

/// Read a JSON following the Aladin JSON format
pub fn from_json_aladin<T, Q>(input: &str) -> Result<CellMOC<T, Q>, Box<dyn Error>>
  where
    T: Idx,
    Q: MocQty<T>,
{
  // See https://docs.serde.rs/serde_json/#operating-on-untyped-json-values
  let root: Value = serde_json::from_str(input)?;
  match &root {
    Value::Object(map) => {
      let mut cells: Vec<Cell<T>> = Vec::with_capacity(
        (0..=Q::MAX_DEPTH).into_iter()
          .filter_map(|d| match map.get(&d.to_string()) {
            Some(Array(vec)) => Some(vec.len()),
            _ => None
          })
          .sum()
      );
      let mut depth_max = 0;
      for depth in 0..=Q::MAX_DEPTH {
        if let Some(Array(vec)) = map.get(&depth.to_string()) {
          for v in vec.iter()
            .filter_map(|v| match v {
              Value::Number(idx) => idx.as_u64().and_then(|v| T::try_from(v).ok()),
              _ => None,
            }) {
            cells.push(Cell::new(depth, v));
          }
          depth_max = depth_max.max(depth);
        }
      }
      // Sort the vec
      cells.sort_by(|a, b| a.cmp::<Q>(b));
      // Check for cell unicity
      // TODO: add unicity check
      Ok(CellMOC::new(depth_max, MocCells::new(Cells(cells))))
    },
    _ => Err(format!("Wrong JSON root type. Expected: Object. Actual: {:?}", &root).into()),
  }
}


// json_stream :
// {
//   "qty": "",
//   "depth": "",
//   "cells": [uniq, uniq, ...]
// }

#[cfg(test)]
mod tests {
  use crate::moc::{RangeMOC, RangeMOCIntoIterator, RangeMOCIterator, CellMOCIterator, CellMOCIntoIterator};
  use crate::mocranges::MocRanges;
  use crate::mocqty::Hpx;
  use crate::deser::json::from_json_aladin;

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
}