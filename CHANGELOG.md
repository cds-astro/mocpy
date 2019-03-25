# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.5.7] - 2019-04-12
### Changed
- Change from_cells to from_healpix_cells. Its API does change too. It now takes the three ipix, depth and flags numpy arrays separatly instead as a numpy structure. This has the advantage of direcltly passing the arrays returned by `cdshealpix`. Creating a numpy structured array from simple numpy column arrays needed to copy the data from the arrays to the structural array.
- Add cdshealpix as a dep. Will remove astropy-healpix dep. When cdshealpix will be merged into astropy-healpix then the astropy-healpix dep will be restored.

## [0.5.6] - 2019-04-11
### Added
- Serialize to str. Call moc.serialize(format="str")
- Load a MOC from a str (see section 2.3.2 of [MOC IVOA paper](http://ivoa.net/documents/MOC/20190215/WD-MOC-1.1-20190215.pdf)).
- Fix vizualization bug when plotting all the allsky MOC. Biggest cells to plot are limited to depth 2. Cells of depth 0 and 1 are
subdivided into cells of depth 2 for the visualization purpose.

### Changed
- Add of a `overwrite` optional keyword in the write method. Before 0.5.6 the default behaviour was to
always overwrite already existing files. Now it does not overwrite by default. To do that, you have to
set the `overwrite` keyword.

## [0.5.5] - 2019-02-08
### Added
- Plotting a moc with matplotlib axes is faster (concers **fill** and **border** methods). The level of detail of the plotted MOC is function of the FoV. The MOC is degraded to the minimum depth so that a cell of this depth can be contained in 1px at the center of the projection. For small FoVs, we only plot the part of MOC contained in the view (thanks to the speed of logical operation between MOCs).
- The [docs](https://mocpy.readthedocs.io/en/latest/) features more examples on how to plot a MOC, perform logical operations between MOCs, etc...
- The doc of the API has been reviewed and features some test codes that can be run with the sphinx command `make doctest`.

### Removed
- The use of **multiprocessing** in the `fill` method.

## [0.5.4] - 2019-02-06
### Added
- Novel python file hierarchy. moc/ and tmoc/ containing the sources for MOC (resp. TMOC) classes.
- A new mocpy.WCS object type that must be declared in a context (with WCS(...) as wcs:) for use. This facilitates the creation of an astropy.wcs.WCS object for plotting a MOC in a matplotlib axes. This replaces the wcs.spatial.utils.make_wcs method that returned an astropy.wcs.WCS object.
- Use of multiprocessing.Pool in **mocpy.MOC.fill**

### Removed
- **wcs.spatial.utils.make_wcs** has been removed. See **mocpy.WCS** as replacement.

## [0.5.3] - 2019-01-30
### Added
- A new constructor from_cells. It returns a new MOC instance from a numpy structured array representing the cells of a MOC. The structured array passed must contains 3 attributes: "ipix": np.uint64, "depth": np.uint32, "fully_covered": np.uint8 (a flag bit. For the moment its value has no effect for the newly created MOC).

## [0.5.2] - 2019-01-28
### Added
- A new **from_polygon_skycoord** method added where you can pass an astropy.coordinates.SkyCoord describing the polygon coordinates instead of two lon, lat astropy.units.Quantity. The *max_depth*, and *inside* optional arguments remain.

### Changed
- Remove spherical geom from dependency so that astroquery.cds wheel for windows/py3 can be generated. Spherical Geom is only used in MOC.from_polygon. A message is addressed to the user telling him to install sphrical geom if it is not installed and if he wants to create a MOC from a polygon.

## [0.5.1] - 2019-01-25
### Changed
- `pip install mocpy` now installs all the dependencies of mocpy. See the setup.py file. (requires changed to install_requires).

## [0.5.0] - 2019-01-09
### Added
- Two methods `fill` and `border` taking an astropy.wcs.WCS and a matplotlib axis. `fill` projects the MOC into the WCS and draws it on the MPL axis using pathpatches for each HEALPix cell. `border` draws the border the same way and requires a WCS and an MPL axe too. You can pass to these functions additional keyword arguments that will directly be passed to MPL when plotting (e.g. the color, the linewidth, and alpha component etc...). Check the notebooks to see how to use these new methods.
- You can retrieve the list of skycoords describing the border MOC. Each skycoord from the list refers to one border of the MOC (either an external border or the border of a hole in a connexe MOC). The process takes for the moment a quite amount of time and thus may be optimized in the future. A new GALEX boundaries notebook has been added so that you can check how it works. I recommend to decrease the order of GALEX to 5 before computing the boundaries otherwise it will take some time. This add relies on the issue [#29](https://github.com/cds-astro/mocpy/issues/29) initiated by [@ggreco77](https://github.com/ggreco77).
- A new `from_polygon` method taking the vertices (i.e. (skycoords) or (lon, lat) tuples) responsible for setting up a new MOC from a polygon.
An inside SkyCoord point is requested and says to the algorithm which area (2 possible as on the sphere, an area and its complement are both finite) must be chosen. If no inside sky coord is given, we consider the mean of the vertices of the polygon as belonging to the MOC (This is without ambiguity for convex polygons but it may not work for concave ones). Vertices describing a convex or concave polygon are accepted. Vertices describing a self-intersecting polygon are not accepted.
This method does not rely on astropy_healpix as there is for the moment no method returning the set of HEALPix cells belonging to a polygon and is thus implemented purely in Python.
- A new `serialize` public method allows to serialize the MOC in two possible format, FITS and json. For a FITS serialization the method returns an astropy HDU binary table. For a JSON serialization, the method returns a python dictionary storing order-[list of HEALPix cells] as key-value pairs.

### Changed
- `write` method does not take a `write_to_file` argument. For serialization purpose, there is now a new `serialize` method.
- astropy_healpix.HEALPix.lonlat_to_healpix seems to not accept astropy MaskedColumn anymore. For lon as a MaskedColumn, please pass lon.T * u.deg to mocpy.MOC.from_lonlat. We need to transpose the column and then convert it to an astropy.units.Quantity.
- Notebooks have been updated and all the plots now use the new methods `fill` and `border`.
- A new package `spatial`, invisible from the user, but keeping all the code of spatial   MOCs (plotting methods, border computation, special utils for creating WCS...) has been created. tmocs and core functions are still located in the root of the project.
- Add of a changelog
