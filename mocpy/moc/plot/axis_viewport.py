from matplotlib.pyplot import figure, axis

def set(ax, wcs):
    # get the plot dimension the same way as in World2ScreenMPL constructor
    # This allows to pass an astropy wcs instead of a World2ScreenMPL when plotting
    # From issue: https://github.com/cds-astro/mocpy/issues/69
    width_px = wcs.wcs.crpix[0] * 2.0
    height_px = wcs.wcs.crpix[1] * 2.0

    x_min = 0
    x_max = width_px - 1
    y_min = 0
    y_max = height_px - 1
    
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])