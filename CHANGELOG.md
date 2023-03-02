# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.12.2]

### Added

* `MOC.MAX_ORDER` and `TimeMOC.MAX_ORDER` to replace former `IntervalSet.HPX_MAX_ORDER` and `IntervalSet.TIME_MAX_ORDER`
* `MOC.to_depth29_ranges` (+test) to replace former `IntervalSet.nested`
* `TimeMOC.to_depth61_ranges`
* tests on `MOC.uniq_hpx`

### Bugfix

* :bug: return statement was missing in `MOC.uniq_hpx`


## [0.12.1]

### Added

* add methods `AbstractMOC.__copy__` and `AbstractMOC.__deepcopy__`
* add op parameter `timeout` in `query_simbad` and `query_vizier_table` that defaults to 1000s
* warning if `max_depth` is not set in `MOC.from_valued_cells`

### Changed

* ⚠️ BREAKING: public function `set` in `plot.axis_viewport` module has been renamed into `_set_wcs` and moved to module `plot.utils`

### Bugfix

* :bug: a bug was introduced in [0.12.0] on the functions `query_simbad` and `query_vizier_table` that are compatible with `pre_v2` format only
* :bug: when `max_depth=None` in `MOC.from_valued_cells`


## [0.12.0]

### Removed

* ⚠️ BREAKING: Deserialisation of ST-MOCs from FITs HDU (reading/writting FITS file fully supported
  on the Rust side)
* ⚠️ BREAKING: Remove support of pre_v2 ST-MOC fits file writing (reading still ok)
* ⚠️ BREAKING: internal class `IntervalSet` removed
* ⚠️ BREAKING: `utils` file removed

### Added

* add `hdu_index` optional parameter in `MOC.from_fits_images`
* `TimeMOC.to_time_ranges` to get time ranges from a T-MOC
* `MOC.to_rgba` and `MOC.dsiplay_preview` for a quick S-MOC allsky view
* `uniq_gen` and `uniq_zorder` added to `AbstractMOC`
* `flatten` added to `AbstractMOC`
* constructor `MOC.from_healpix_depth29_ranges`
* parameter `values_are_densities` in `MOC.from_valued_healpix_cells`
* parameter `complement` in `MOC.from_polygon`
* `+`, `|`, `-`, `&`, `~` operators redefinition for `union`, `union`, `difference`, `intersect` and `complement` respectively.
* `contains_skycoords` and `contains_lonlat` to replace `contains`
* add `fold` parameter into `save` and `to_string`
* add `MOC.barycenter` and `MOC.largest_distance_from_coo_to_vertices` of a moc
* add `MOC.wcs` giving an astropy WCS centered around a the barycenter of moc

### Changed

* ⚠️ BREAKING: `times_to_microseconds` and `microseconds_to_times` moved from `utils` to `tmoc`.
* ⚠️ BREAKING: `uniq` removed from `IntervalSet`, but replacing method `uniq_hpx` added to `MOC`
    +  ⚠️ BREAKING: the output of `uniq_hpx` is not sorted (but follow the order of the cells in the internal range list)
* ⚠️ BREAKING: `STMOC.query_by_time` now takes in input a `TimeMOC`
* ⚠️ BREAKING: `STMOC.query_by_space` now returns a `MOC`
* ⚠️ BREAKING: `TimeMOC.contains` does not take any longer a time resolution as input parameter
* ⚠️ BREAKING: `TimeMOC.contains_with_timeresolution` as been added with the previous behaviour of  `TimeMOC.contains`
* add `save` to `AbstractMOC` and remove from sub-classes
* add `to_string` to `AbstractMOC` and remove from sub-classes
* `from_uniq` removed from `IntervalSet` and added to `MOC`
* change the `contains` implementation to be much memory efficient, faster, and thus working at all HEALPix orders.
* update `cdshealpix` and `moc` dependencies
* ⚠️ BREAKING: `MOC.from_healpix_cells`
  * now requires the `max_depth`, the depth of the MOC we want to create
  * optional parameter `fully_covered` removed since it was not used
* ⚠️ BREAKING: in `MOC.from_valued_healpix_cells`, the user now have to degrade manually the resolution
  if max_depth < deepest cell depth.
* ⚠️ BREAKING: `World2ScreenMPL` has been renamed `WCS`

## [0.11.0]

### Added

* Add `MOC.from_ring`
* Option `include_indirect_neighbours` to `split` and `split_count`

### Changed

* Extend the `moc.fill` and `moc.border` to directly accept an astropy wcs object. This solves the issue <https://github.com/cds-astro/mocpy/issues/69>
* Addresses the plotting artefacts when plotting big HEALPix cells. Cells of depth < 3 are subdivided to the depth 3 before being plotted
* Set the default `time_order` from 29 to 61, i.e. to `max_order` in `stmoc.from_spatial_coverages` (29 was the `max_order` in the previous MOC V2.0 version).

### Improvement

* More robust FITS UNIQ deserialization (ignoring 0 values)

## [0.10.0]

### Changed

* Rename `TimeMOC` logical operations taking a deltaT constraint adding the sufix `_with_timeresolution`
* WARNING logical `TimeMOC` logical operations are now at the deepest depth (no time resolution parameter)

### Added

* Deprecate `from_fits`, ...  methods
* Add `MOC.spli_count`, `MOC.split`, `MOC.from_multiordermap_fits_file`
* Add support for u16 and u32 fits MOC and TMOC in 'load'

### Bug correction

* Replace empty moc shape (1, 0) by (0, 2)
* Fix `tmoc.max_time`
* Fix doc (due to an evolution of sphinx)

## [0.9.0]

### Changed

* Add compatibility with MOC2.0: Time depth in now in `[0, 61]` instead of `[0, 29]`
* Add `from_time_ranges_in_microsec_since_jd_origin` in `temporal_coverage.rs`
* Time operations now precise to the microseconds, see:
  * `utils.times_to_microseconds`
    * `utils.microseconds_to_times`
* Add several options to `moc.from_valued_healpix_cells`
* Add `load` and `save` methods for serialization/deserialization in pure Rust (ensuring
  MOC2.0 compatibility).
* Improve performance and some operations (like intersection and union)
* Improve `to_uniq` performances (x5 according to a bench in pure Rust)
* Improve `add_neighbours` and `remove_neighbours` performances (now in pure Rust)

### Internal Python changes

* `complement()`: remove from AbstractMoc / IntervalSet, and add it in Moc and TimeMoc (because complement now depends on the qty)

### Internal changes

* Remove the `moc` crate from MOCPy and publish it as a standalone [here](https://github.com/cds-astro/cds-moc-rust)
  (with a lot of added stuff)
* Add FITS serialization/deserialization on the Rust code.
* Add ASCII and JSON serialization/deserialization on the Rust code.
* Move `rand` from dependencies to dev-dependenciies
* Generalize the code to support various quantities with different dimensions (HEALPix indices, Time, ...)
  * create `MocableQty` and `MocQty` implemented by `Hpx` and `Time`
* Remove depth/qty dependent operations from `Ranges` (depth/complement/degrade), create a trait for generic operations
* Add `MocRange(s)` since we introduced `MocQty` for depth dependent operations, and introduce both: `HpxRange(s)` and `TimeRange(s)`
* Add `MocRanges2D` for depth dependent operations
* Rename `NestedXX` in `HpxXX` to possibly support Ring indices (the code should be the same as far as the NSIDE is a power of 2)
* ...

## [0.8.5]

### Changed

* change the CI: replace travis and appveyor by github actions
* replace setuptools rust by maturin
* update dependencies

## [0.8.2]

### Changed

* remove ',' separator when deserializing MOC from ascii (this follows the MOC 1.1 standard [http://ivoa.net/documents/MOC/20191007/REC-MOC-1.1-20191007.pdf](http://ivoa.net/documents/MOC/20191007/REC-MOC-1.1-20191007.pdf))

## [0.8.1]

### Added

* from_valued_healpix_cells

### Changed

* API Breaking change! from_image -> from_fits_image(hdulist, max_norder)
* WCS -> World2ScreenMPL. It's a context manager class

## [0.7.4]

### Added

* Change API for ST-MOC: query_by_time, query_by_space

## [0.7.0]

### Added

* Space-Time coverages, classmethod from creating them from (time, ra, dec) tuples
* Query a Space-Time coverages with time frames and spatial coverages

## [0.6.0]

### Added

* Rust backend
* Add tests for windows py27 and py3
* from_polygon relies on cdshealpix. spherical_geometry dependency removed!
* change astropy-healpix dependency from mandatory to optional. astropy-healpix is now used in only a few deprecated methods
  (such as the old plot method from `mocpy.MOC` and `get_boundaries` which will soon make use of cdshealpix too).

### Changed

* API CHANGE!: the ``inside`` parameter of from_polygon and from_polygon_skycoord has been removed !
  The inside of the polygon is deduced from the order of the sky coordinates passed.

## [0.5.7] - 2019-04-12

### Changed

* Change from_cells to from_healpix_cells. Its API does change too. It now takes the three ipix, depth and flags numpy arrays separatly instead as a numpy structure. This has the advantage of direcltly passing the arrays returned by `cdshealpix`. Creating a numpy structured array from simple numpy column arrays needed to copy the data from the arrays to the structural array.
* Add cdshealpix as a dep. Will remove astropy-healpix dep. When cdshealpix will be merged into astropy-healpix then the astropy-healpix dep will be restored.

## [0.5.6] - 2019-04-11

### Added

* Serialize to str. Call moc.serialize(format="str")
* Load a MOC from a str (see section 2.3.2 of [MOC IVOA paper](http://ivoa.net/documents/MOC/20190215/WD-MOC-1.1-20190215.pdf)).
* Fix vizualization bug when plotting all the allsky MOC. Biggest cells to plot are limited to depth 2. Cells of depth 0 and 1 are
subdivided into cells of depth 2 for the visualization purpose.

### Changed

* Add of a `overwrite` optional keyword in the write method. Before 0.5.6 the default behaviour was to
always overwrite already existing files. Now it does not overwrite by default. To do that, you have to
set the `overwrite` keyword.

## [0.5.5] - 2019-02-08

### Added

* Plotting a moc with matplotlib axes is faster (concers **fill** and **border** methods). The level of detail of the plotted MOC is function of the FoV. The MOC is degraded to the minimum depth so that a cell of this depth can be contained in 1px at the center of the projection. For small FoVs, we only plot the part of MOC contained in the view (thanks to the speed of logical operation between MOCs).
* The [docs](https://mocpy.readthedocs.io/en/latest/) features more examples on how to plot a MOC, perform logical operations between MOCs, etc...
* The doc of the API has been reviewed and features some test codes that can be run with the sphinx command `make doctest`.

### Removed

* The use of **multiprocessing** in the `fill` method.

## [0.5.4] - 2019-02-06

### Added

* Novel python file hierarchy. moc/ and tmoc/ containing the sources for MOC (resp. TMOC) classes.
* A new mocpy.WCS object type that must be declared in a context (with WCS(...) as wcs:) for use. This facilitates the creation of an astropy.wcs.WCS object for plotting a MOC in a matplotlib axes. This replaces the wcs.spatial.utils.make_wcs method that returned an astropy.wcs.WCS object.
* Use of multiprocessing.Pool in **mocpy.MOC.fill**

### Removed

* **wcs.spatial.utils.make_wcs** has been removed. See **mocpy.WCS** as replacement.

## [0.5.3] - 2019-01-30

### Added

* A new constructor from_cells. It returns a new MOC instance from a numpy structured array representing the cells of a MOC. The structured array passed must contains 3 attributes: "ipix": np.uint64, "depth": np.uint32, "fully_covered": np.uint8 (a flag bit. For the moment its value has no effect for the newly created MOC).

## [0.5.2] - 2019-01-28

### Added

* A new **from_polygon_skycoord** method added where you can pass an astropy.coordinates.SkyCoord describing the polygon coordinates instead of two lon, lat astropy.units.Quantity. The *max_depth*, and *inside* optional arguments remain.

### Changed

* Remove spherical geom from dependency so that astroquery.cds wheel for windows/py3 can be generated. Spherical Geom is only used in MOC.from_polygon. A message is addressed to the user telling him to install sphrical geom if it is not installed and if he wants to create a MOC from a polygon.

## [0.5.1] - 2019-01-25

### Changed

* `pip install mocpy` now installs all the dependencies of mocpy. See the setup.py file. (requires changed to install_requires).

## [0.5.0] - 2019-01-09

### Added

* Two methods `fill` and `border` taking an astropy.wcs.WCS and a matplotlib axis. `fill` projects the MOC into the WCS and draws it on the MPL axis using pathpatches for each HEALPix cell. `border` draws the border the same way and requires a WCS and an MPL axe too. You can pass to these functions additional keyword arguments that will directly be passed to MPL when plotting (e.g. the color, the linewidth, and alpha component etc...). Check the notebooks to see how to use these new methods.
* You can retrieve the list of skycoords describing the border MOC. Each skycoord from the list refers to one border of the MOC (either an external border or the border of a hole in a connexe MOC). The process takes for the moment a quite amount of time and thus may be optimized in the future. A new GALEX boundaries notebook has been added so that you can check how it works. I recommend to decrease the order of GALEX to 5 before computing the boundaries otherwise it will take some time. This add relies on the issue [#29](https://github.com/cds-astro/mocpy/issues/29) initiated by [@ggreco77](https://github.com/ggreco77).
* A new `from_polygon` method taking the vertices (i.e. (skycoords) or (lon, lat) tuples) responsible for setting up a new MOC from a polygon.
An inside SkyCoord point is requested and says to the algorithm which area (2 possible as on the sphere, an area and its complement are both finite) must be chosen. If no inside sky coord is given, we consider the mean of the vertices of the polygon as belonging to the MOC (This is without ambiguity for convex polygons but it may not work for concave ones). Vertices describing a convex or concave polygon are accepted. Vertices describing a self-intersecting polygon are not accepted.
This method does not rely on astropy_healpix as there is for the moment no method returning the set of HEALPix cells belonging to a polygon and is thus implemented purely in Python.
* A new `serialize` public method allows to serialize the MOC in two possible format, FITS and json. For a FITS serialization the method returns an astropy HDU binary table. For a JSON serialization, the method returns a python dictionary storing order-[list of HEALPix cells] as key-value pairs.

### Changed

* `write` method does not take a `write_to_file` argument. For serialization purpose, there is now a new `serialize` method.
* astropy_healpix.HEALPix.lonlat_to_healpix seems to not accept astropy MaskedColumn anymore. For lon as a MaskedColumn, please pass lon.T * u.deg to mocpy.MOC.from_lonlat. We need to transpose the column and then convert it to an astropy.units.Quantity.
* Notebooks have been updated and all the plots now use the new methods `fill` and `border`.
* A new package `spatial`, invisible from the user, but keeping all the code of spatial   MOCs (plotting methods, border computation, special utils for creating WCS...) has been created. tmocs and core functions are still located in the root of the project.
* Add of a changelog

