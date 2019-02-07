import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

def build_plotting_moc(moc, wcs):
    # Get the WCS cdelt giving the deg.px^(-1) resolution.
    cdelt = wcs.wcs.cdelt
    # Convert in rad.px^(-1)
    cdelt = np.abs((2*np.pi/360)*cdelt[0])
    # Get the minimum depth such as the resolution of a cell is contained in 1px. 
    depth_res = int(np.floor(np.log2(np.sqrt(np.pi/3)/cdelt)))
    depth_res = max(depth_res, 0)
    # Degrade the moc to that depth for plotting purposes. It is not necessary to plot pixels
    # that we will not see because they are contained in 1px.
    moc_plot = moc
    if moc.max_order > depth_res:
        moc_plot = moc.degrade_to_order(depth_res)

    # Get the MOC delimiting the FOV polygon
    width_px = int(wcs.wcs.crpix[0]*2.) # Supposing the wcs is centered in the axis
    heigth_px = int(wcs.wcs.crpix[1]*2.)

    # Compute the sky coordinate path delimiting the viewport.
    # It consists of a closed polygon of (4 - 1)*4 = 12 vertices
    x_px = np.linspace(0, width_px, 4)
    y_px = np.linspace(0, heigth_px, 4)

    X, Y = np.meshgrid(x_px, y_px)

    X_px = np.append(X[0, :-1], X[:-1, -1])
    X_px = np.append(X_px, np.flip(X[-1, 1:]))
    X_px = np.append(X_px, X[:-1, 0])
    
    Y_px = np.append(Y[0, :-1], Y[:-1, -1])
    Y_px = np.append(Y_px, Y[-1, :-1])
    Y_px = np.append(Y_px, np.flip(Y[1:, 0]))

    assert(X_px.shape == Y_px.shape)
    # Inverse projection from pixel coordinate space to the world coordinate space
    radec = wcs.all_pix2world(np.vstack((X_px, Y_px)).T, 0)
    ra, dec = radec[:, 0], radec[:, 1]
    # If one coordinate is a NaN we exit the function and do not go further
    if np.isnan(ra).any() or np.isnan(dec).any():
        return moc_plot

    center_x_px, center_y_px = wcs.wcs.crpix[0], wcs.wcs.crpix[1]
    ra_center, dec_center = wcs.all_pix2world(center_x_px, center_y_px, 0)

    frame = "icrs" if wcs.is_celestial else "galactic"
    viewport = SkyCoord(ra, dec, unit="deg", frame=frame)
    inside = SkyCoord(ra_center * u.deg, dec_center * u.deg, frame=frame, unit="deg")

    from ..moc import MOC
    # Create a rough MOC (depth=3 is sufficient) from the viewport
    moc_viewport = MOC.from_polygon_skycoord(viewport, max_depth=3, inside=inside)
    
    # The moc to plot is the INPUT_MOC & MOC_VIEWPORT. For small FOVs this can reduce
    # a lot the time to draw the MOC along with its borders.
    moc_plot = moc_plot.intersection(moc_viewport)
    return moc_plot

