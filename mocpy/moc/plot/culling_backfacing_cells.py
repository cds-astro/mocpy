import numpy as np

from astropy_healpix import HEALPix
from astropy.coordinates import ICRS
from astropy.wcs.utils import skycoord_to_pixel

def backface_culling(xp, yp):
    # Remove cells crossing the MOC after projection
    # The remaining HEALPix cells are used for computing the patch of the MOC
    vx = xp
    vy = yp

    def cross_product(vx, vy, i):
        cur = i
        prev = (i - 1) % 4
        next = (i + 1) % 4

        # Construct the first vector from A to B
        x1 = vx[:, cur] - vx[:, prev]
        y1 = vy[:, cur] - vy[:, prev]
        z1 = np.zeros(x1.shape)

        v1 = np.vstack((x1, y1, z1)).T
        # Construct the second vector from B to C
        x2 = vx[:, next] - vx[:, cur]
        y2 = vy[:, next] - vy[:, cur]
        z2 = np.zeros(x2.shape)

        v2 = np.vstack((x2, y2, z2)).T
        # Compute the cross product between the two
        return np.cross(v1, v2)

    # A ----- B
    #  \      |
    #   D-----C
    # Compute the cross product between AB and BC
    # and the cross product between BC and CD
    ABC = cross_product(vx, vy, 1)
    CDA = cross_product(vx, vy, 3)

    frontface_cells  = (ABC[:, 2] < 0) & (CDA[:, 2] < 0)
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
        hp = HEALPix(nside=(1 << depth), order='nested', frame=ICRS())

        ipix_boundaries = hp.boundaries_skycoord(ipixels, step=1)
        # Projection on the given WCS
        xp, yp = skycoord_to_pixel(coords=ipix_boundaries, wcs=wcs)
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