import numpy as np

from astropy.coordinates import ICRS
from astropy.wcs.utils import skycoord_to_pixel

import cdshealpix
from astropy.coordinates import SkyCoord

def backface_culling(xp, yp):
    # Remove cells crossing the MOC after projection
    # The remaining HEALPix cells are used for computing the patch of the MOC
    vx = xp
    vy = yp

    def compute_vector_at_index(X, Y, i):
        cur = i
        next = (i + 1) % 4

        x = X[:, next] - X[:, cur]
        y = Y[:, next] - Y[:, cur]
        z = np.zeros(x.shape)

        v = np.vstack((x, y, z)).T
        return v

    # A ----- B
    #  \      |
    #   D-----C
    # Compute the cross product between AB and BC
    # and the cross product between BC and CD
    AB = compute_vector_at_index(vx, vy, 0)
    BC = compute_vector_at_index(vx, vy, 1)
    CD = compute_vector_at_index(vx, vy, 2)
    DA = compute_vector_at_index(vx, vy, 3)

    ABC = np.cross(AB, BC)
    BCD = np.cross(BC, CD)
    CDA = np.cross(CD, DA)
    DAB = np.cross(DA, AB)

    tol = -0.1
    frontface_cells = (ABC[:, 2] <= tol) & (CDA[:, 2] <= tol) & (BCD[:, 2] <= tol) & (DAB[:, 2] <= tol)

    frontface_cells = np.asarray(frontface_cells)
    vx = vx[frontface_cells]
    vy = vy[frontface_cells]

    return vx, vy, frontface_cells

def from_moc(depth_ipix_d, wcs):
    # Create a new MOC that do not contain the HEALPix
    # cells that are backfacing the projection
    depths = [int(depth_str) for depth_str in depth_ipix_d.keys()]
    min_depth = min(depths)
    max_depth = max(depths)
    ipixels = np.asarray(depth_ipix_d[str(min_depth)])

    # Split the cells located at the border of the projection
    # until at least the depth 7
    max_split_depth = max(7, max_depth)

    ipix_d = {}
    for depth in range(min_depth, max_split_depth + 1):
        ipix_lon, ipix_lat = cdshealpix.vertices(ipixels, depth)
        
        ipix_lon = ipix_lon[:, [2, 3, 0, 1]]
        ipix_lat = ipix_lat[:, [2, 3, 0, 1]]
        ipix_vertices = SkyCoord(ipix_lon, ipix_lat, frame=ICRS())

        # Projection on the given WCS
        xp, yp = skycoord_to_pixel(coords=ipix_vertices, wcs=wcs)
        _, _, frontface_id = backface_culling(xp, yp)

        # Get the pixels which are backfacing the projection
        backfacing_ipix = ipixels[~frontface_id]
        frontface_ipix = ipixels[frontface_id]

        depth_str = str(depth)
        ipix_d.update({depth_str: frontface_ipix})

        next_depth = str(depth + 1)
        ipixels = []
        
        if next_depth in depth_ipix_d:
            ipixels = depth_ipix_d[next_depth]

        for bf_ipix in backfacing_ipix:
            child_bf_ipix = bf_ipix << 2
            ipixels.extend([child_bf_ipix,
                child_bf_ipix + 1,
                child_bf_ipix + 2,
                child_bf_ipix + 3])

        ipixels = np.asarray(ipixels)

    return ipix_d