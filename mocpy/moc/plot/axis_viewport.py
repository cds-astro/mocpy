from matplotlib.pyplot import figure, axis

def set(ax, wcs):
    fig = ax.get_figure()
    width_px, height_px = fig.get_size_inches() * fig.dpi

    x_min = 0
    x_max = width_px - 1
    y_min = 0
    y_max = height_px - 1
    
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])