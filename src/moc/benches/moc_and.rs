
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use criterion::{Criterion, criterion_group, criterion_main};

use healpix::nested::bmoc::BMOC;

use moc::qty::Hpx;
use moc::moc::range::RangeMOC;
use moc::deser::fits::{
  from_fits_ivoa,
  MocIdxType, MocQtyType, MocType
};
use moc::moc::{
  HasMaxDepth,
  IntoBMOC,
  RangeMOCIterator, RangeMOCIntoIterator,
  CellMOCIterator, CellMOCIntoIterator
};
use moc::moc::range::op::and::and;
use moc::ranges::SNORanges;

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

fn test_and_ranges(moc_l: RangeMOC<u32, Hpx<u32>>, moc_r: RangeMOC<u32, Hpx<u32>>) -> RangeMOC<u32, Hpx<u32>> {
  let depth = u8::max(moc_l.depth_max(), moc_r.depth_max());
  let ranges_l = moc_l.into_moc_ranges();
  let ranges_r = moc_r.into_moc_ranges();
  RangeMOC::new(depth, ranges_l.intersection(&ranges_r))
}

// we could also perform the operation without having first collected the iteartor we obtain from
// the FITS file
fn test_and_ranges_it(moc_l: RangeMOC<u32, Hpx<u32>>, moc_r: RangeMOC<u32, Hpx<u32>>) -> RangeMOC<u32, Hpx<u32>> {
  let and = and(moc_l.into_range_moc_iter(), moc_r.into_range_moc_iter());
  RangeMOC::new(and.depth_max(), and.collect())
}

fn test_and_ranges_it_ref(moc_l: RangeMOC<u32, Hpx<u32>>, moc_r: RangeMOC<u32, Hpx<u32>>) -> RangeMOC<u32, Hpx<u32>> {
  let and = and((&moc_l).into_range_moc_iter(), (&moc_r).into_range_moc_iter());
  RangeMOC::new(and.depth_max(), and.collect())
}

fn test_and_bmoc(moc_l: BMOC, moc_r: BMOC) -> BMOC {
  moc_l.and(&moc_r)
}

fn bench_and(c: &mut Criterion) {
  // https://bheisler.github.io/criterion.rs/book/user_guide/comparing_functions.html
  let mut group = c.benchmark_group("and");
  let (sdss, other) = load_mocs();
  // let sdss_bmoc = (&sdss).into_moc_ranges().cells().into_bmoc();
  // let other_bmoc = (&other).into_moc_ranges().cells().into_bmoc();
  group.bench_function("Ranges INTERSECTION",
                       |b| b.iter(|| test_and_ranges(sdss.clone(), other.clone())));
  group.bench_function("Ranges Iter AND",
                       |b| b.iter(|| test_and_ranges_it(sdss.clone(), other.clone())));
  group.bench_function("Ranges Ref Iter AND",
                       |b| b.iter(|| test_and_ranges_it_ref(sdss.clone(), other.clone())));
  group.bench_function("BMOC AND",
                       |b| b.iter(|| test_and_bmoc(
                         (&sdss).into_range_moc_iter().cells().into_bmoc(),
                         (&other).into_range_moc_iter().cells().into_bmoc()
                       )));
  /*group.bench_function("Ranges 2 INTERSECTION",
                       |b| b.iter(|| test_and_ranges(sdss.clone(), other.clone())));
  group.bench_function("Ranges Iter 2 AND",
                       |b| b.iter(|| test_and_ranges_it(sdss.clone(), other.clone())));
  group.bench_function("Ranges Ref Iter 2 AND",
                       |b| b.iter(|| test_and_ranges_it_ref(sdss.clone(), other.clone())));*/

}

criterion_group!(benches, bench_and);
criterion_main!(benches);
