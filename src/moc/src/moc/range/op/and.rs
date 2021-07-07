
use std::ops::Range;
use std::cmp::Ordering;

use crate::idx::Idx;
use crate::qty::MocQty;
use crate::moc::{HasMaxDepth, ZSorted, NonOverlapping, MOCProperties, RangeMOCIterator};

/// Performs a logical `AND` between the two input iterators of ranges.
pub fn and<T, Q, I1, I2>(
  left_it: I1, 
  right_it: I2
) -> AndRangeIter<T, Q, I1, I2>
  where 
    T: Idx,
    Q: MocQty<T>,
    I1: RangeMOCIterator<T, Qty=Q>,
    I2: RangeMOCIterator<T, Qty=Q> 
{
  AndRangeIter::new(left_it, right_it)
}

/// Performs an `AND` operation between two iterators of ranges on-the-fly, while iterating.
pub struct AndRangeIter<T, Q, I1, I2>
  where
    T: Idx,
    Q: MocQty<T>,
    I1: RangeMOCIterator<T, Qty=Q>,
    I2: RangeMOCIterator<T, Qty=Q> 
{
  left_it: I1,
  right_it: I2,
  left: Option<Range<T>>,
  right: Option<Range<T>>,
}

impl <T, Q, I1, I2> AndRangeIter<T, Q, I1, I2>
  where
    T: Idx,
    Q: MocQty<T>,
    I1: RangeMOCIterator<T, Qty=Q>,
    I2: RangeMOCIterator<T, Qty=Q>
{
  
  fn new(mut left_it: I1, mut right_it: I2) -> AndRangeIter<T, Q, I1, I2> {
    let left = left_it.next();
    let right = right_it.next();
    // Quick rejection tests
    if let (Some(up_left), Some(low_right)) = (left_it.peek_last(), &right) {
      if up_left.end <= low_right.start {
        return AndRangeIter { left_it, right_it, left: None, right: None };
      }
    } 
    if let (Some(low_left), Some(up_right)) = (&left, right_it.peek_last()) {
      if up_right.end <= low_left.start {
        return AndRangeIter { left_it, right_it, left: None, right: None };
      }
    }
    // Normal behaviour
    AndRangeIter {
      left_it,
      right_it,
      left,
      right,
    }
  }

}

impl<T, Q, I1, I2> HasMaxDepth for  AndRangeIter<T, Q, I1, I2>
  where
    T: Idx,
    Q: MocQty<T>,
    I1: RangeMOCIterator<T, Qty=Q>,
    I2: RangeMOCIterator<T, Qty=Q>
{
  fn depth_max(&self) -> u8 {
    u8::max(self.left_it.depth_max(), self.right_it.depth_max())
  }
}
impl<T, Q, I1, I2> ZSorted for  AndRangeIter<T, Q, I1, I2>
  where
    T: Idx,
    Q: MocQty<T>,
    I1: RangeMOCIterator<T, Qty=Q>,
    I2: RangeMOCIterator<T, Qty=Q> 
{ }
impl<T, Q, I1, I2> NonOverlapping for AndRangeIter<T, Q, I1, I2>
  where
    T: Idx,
    Q: MocQty<T>,
    I1: RangeMOCIterator<T, Qty=Q>,
    I2: RangeMOCIterator<T, Qty=Q> 
{ }
impl<T, Q, I1, I2> MOCProperties for AndRangeIter<T, Q, I1, I2>
  where
    T: Idx,
    Q: MocQty<T>,
    I1: RangeMOCIterator<T, Qty=Q>,
    I2: RangeMOCIterator<T, Qty=Q>
{ }
impl<T, Q, I1, I2> Iterator for AndRangeIter<T, Q, I1, I2>
  where 
    T: Idx,
    Q: MocQty<T>,
    I1: RangeMOCIterator<T, Qty=Q>,
    I2: RangeMOCIterator<T, Qty=Q>
{

  type Item = Range<T>;

  fn next(&mut self) -> Option<Self::Item> {
    while let (Some(l), Some(r)) = (&self.left, &self.right) {
      /*let from = l.start.max(r.start);
      let l_to = l.end;
      let r_to = r.end;
      let to = match l_to.cmp(&r_to) {
        Ordering::Less => {
          self.left = self.left_it.next();
          l_to
        },
        Ordering::Greater => {
          self.right = self.right_it.next();
          r_to
        },
        Ordering::Equal => {
          self.left = self.left_it.next();
          self.right = self.right_it.next();
          l_to
        },
      };
      if from < to {
        return Some(from..to);
      }*/
      if l.end <= r.start { // |--l--| |--r--|
        self.left = self.left_it.next();
        while let Some(l) = &self.left {
          if l.end <= r.start {
            self.left = self.left_it.next();
          } else {
            break;
          }
        }
      } else if r.end <= l.start { // |--r--| |--l--|
        self.right = self.right_it.next();
        while let Some(r) = &self.right {
          if r.end <= l.start {
            self.right = self.right_it.next();
          } else {
            break;
          }
        }
      } else {
        let from = l.start.max(r.start);
        return Some(match l.end.cmp(&r.end) {
          Ordering::Less => {
            let range = from..l.end;
            self.left = self.left_it.next();
            range
          },
          Ordering::Greater => {
            let range = from..r.end;
            self.right = self.right_it.next();
            range
          },
          Ordering::Equal => {
            let range = from..l.end;
            self.left = self.left_it.next();
            self.right = self.right_it.next();
            range
          }
        })
      }
    }
    None
  }

  fn size_hint(&self) -> (usize, Option<usize>) {
    let size_hint_l = self.left_it.size_hint();
    let size_hint_r = self.right_it.size_hint();
    if let (Some(n1), Some(n2)) = (size_hint_l.1, size_hint_r.1) {
      // Worse case:
      // * 1 range of T1 contains (n2-1) T2 ranges
      // * and 1 range of T2 contains the remaining (n1-1) cells of T1
      (0, Some(1 + n1 + n2))
    } else {
      (0, None)
    }
  }

}


impl<T, Q, I1, I2> RangeMOCIterator<T> for AndRangeIter<T, Q, I1, I2>
  where
    T: Idx,
    Q: MocQty<T>,
    I1: RangeMOCIterator<T, Qty=Q>,
    I2: RangeMOCIterator<T, Qty=Q>
{
  type Qty = Q;

  fn peek_last(&self) -> Option<&Range<T>> {
    // We could have considered the case in which the upper bound is the same for both inputs
    None
  }
}



#[cfg(test)]
mod tests {
  use std::fs::File;
  use std::io::BufReader;
  use std::path::PathBuf;

  use crate::qty::Hpx;
  use crate::moc::range::RangeMOC;
  use crate::deser::fits::{
    from_fits_ivoa,
    MocIdxType, MocQtyType, MocType
  };
  use crate::moc::{HasMaxDepth, RangeMOCIntoIterator, CellMOCIntoIterator, CellMOCIterator};
  use crate::moc::range::op::and::and;
  use crate::ranges::SNORanges;

  fn load_moc(filename: &str) -> RangeMOC<u32, Hpx<u32>> {
    let path_buf1 = PathBuf::from(format!("resources/{}", filename));
    let path_buf2 = PathBuf::from(format!("../resources/{}", filename));
    let file = File::open(&path_buf1).or_else(|_| File::open(&path_buf2)).unwrap();
    let reader = BufReader::new(file);
    match from_fits_ivoa(reader) {
      Ok(MocIdxType::U32(MocQtyType::Hpx(MocType::Ranges(moc)))) => {
        let moc = RangeMOC::new(moc.depth_max(), moc.collect());
        moc
      },
      Ok(MocIdxType::U32(MocQtyType::Hpx(MocType::Cells(moc)))) => {
        let moc = RangeMOC::new(
          moc.depth_max(),
          moc.into_cell_moc_iter().ranges().collect()
        );
        moc
      },
      _ => unreachable!(false),
    }
  }

  fn load_mocs() -> (RangeMOC<u32, Hpx<u32>>, RangeMOC<u32, Hpx<u32>>) {
    let sdss = load_moc("V_147_sdss12.moc.fits");
    let other = load_moc("CDS-I-125A-catalog_MOC.fits");
    (sdss, other)
  }

  fn and_ranges(moc_l: RangeMOC<u32, Hpx<u32>>, moc_r: RangeMOC<u32, Hpx<u32>>) -> RangeMOC<u32, Hpx<u32>> {
    let depth = u8::max(moc_l.depth_max(), moc_r.depth_max());
    let ranges_l = moc_l.into_moc_ranges();
    let ranges_r = moc_r.into_moc_ranges();
    RangeMOC::new(depth, ranges_l.intersection(&ranges_r))
  }

  // we could also perform the operation without having first collected the iteartor we obtain from
  // the FITS file
  fn and_ranges_it(moc_l: RangeMOC<u32, Hpx<u32>>, moc_r: RangeMOC<u32, Hpx<u32>>) -> RangeMOC<u32, Hpx<u32>> {
    let and = and(moc_l.into_range_moc_iter(), moc_r.into_range_moc_iter());
    RangeMOC::new(and.depth_max(), and.collect())
  }

  #[test]
  fn test_and() {
    let (sdss, other) = load_mocs();
    let and_1 = and_ranges(sdss.clone(), other.clone());
    let and_2 = and_ranges_it(sdss.clone(), other.clone());
    println!("and size: {}", and_1.moc_ranges().0.0.len());
    assert_eq!(and_1.moc_ranges().0.0.len(), and_2.moc_ranges().0.0.len());
  }
}
