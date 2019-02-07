import numpy as np

from astropy import units as u
from astropy.coordinates import ICRS, Galactic
from astropy.coordinates import SkyCoord

from astropy_healpix import HEALPix
from astropy_healpix.core import boundaries_lonlat

try:
    from spherical_geometry.polygon import SphericalPolygon
    from spherical_geometry import great_circle_arc
    from spherical_geometry import vector
except ImportError as e:
    print("Unable to find spherical-geometry package installed. See https://github.com/spacetelescope/spherical_geometry")
    raise

class PolygonComputer:
    def polygon_contains_ipix(self, ipix):
        def poly_contains_ipix_vertices(poly, ipix):
            # Get the ipix polygon vertices as a Nx3 numpy array
            # Remove the last vertex as it counts double (closed polygon)
            ipix_vertices = list(ipix.points)[0][:-1, :]

            for ipix_vertex in ipix_vertices:
                if not poly.contains_point(ipix_vertex):
                    return False

            return True

        def poly_crossing_ipix(poly, ipix):
            # ipix and poly are spherical shapes
            poly_points = list(poly.points)[0]
            x_poly, y_poly, z_poly = (poly_points[:, 0], poly_points[:, 1], poly_points[:, 2])

            ipix_points = list(ipix.points)[0]
            x_ipix, y_ipix, z_ipix = (ipix_points[:, 0], ipix_points[:, 1], ipix_points[:, 2])

            A_poly = np.stack((x_poly, y_poly, z_poly)).T
            B_poly = A_poly
            B_poly = np.append(B_poly, [B_poly[1]], axis=0)
            B_poly = B_poly[1:]

            A_poly = A_poly[:-1]
            B_poly = B_poly[:-1]

            for i in range(len(ipix_points) - 1):
                A_ipix = (x_ipix[i], y_ipix[i], z_ipix[i])
                B_ipix = (x_ipix[i+1], y_ipix[i+1], z_ipix[i+1])

                inter = great_circle_arc.intersects(A_poly, B_poly, A_ipix, B_ipix)
                if inter.any():
                    return True
            return False

        return poly_contains_ipix_vertices(self.polygon, ipix) and \
         (not poly_contains_ipix_vertices(ipix, self.polygon)) and \
         (not poly_crossing_ipix(self.polygon, ipix))

    def _get_starting_depth(self):
        def compute_angular_distance(n1, n2):
            return np.arctan(np.linalg.norm(np.cross(n1, n2))/np.dot(n1, n2))

        def to_xyz(lon, lat):
            x = np.cos(lon) * np.cos(lat)
            y = np.cos(lat) * np.sin(lon)
            z = np.sin(lat)

            return np.array([x, y, z], dtype=np.float64)

        def max_distance_center_to_vertex(depth):
            nside = 1 << depth

            lat1 = np.arcsin(2 / 3.0)
            lat2 = np.arcsin(1 - ((1 - 1.0/nside)**2 / 3.0))
            lon1 = np.pi/(4 * nside)
            lon2 = 0

            # Convert lon, lat to xyz, vector
            n1 = to_xyz(lon=lon1, lat=lat1)
            n2 = to_xyz(lon=lon2, lat=lat2)

            dist = compute_angular_distance(n1, n2)
            return dist # in rad

        # Get the polygon vertices as a Nx3 numpy array
        # Remove the last vertex as it counts double (closed polygon)
        p_vertices = np.asarray(list(self.polygon.points))[0][:-1, :]
        # Get the center formed by the vertices
        center = p_vertices.mean(axis=0)

        # Normalize it so that it lies on the unit sphere
        vector.normalize_vector(center, output=center)
        center = np.asarray(center)
        # Compute the maximum angular distance between the polygon vertices
        # and its center. This refers to the Vector version of the
        # Great-circle distance computation.
        #
        # See https://en.wikipedia.org/wiki/Great-circle_distance
        max_d = -1

        # Check if the polygon covers more than an hemisphere
        covers_more_than_one_hemisphere = (self.polygon.area() > 2 * np.pi)

        for vertex in p_vertices:
            d = compute_angular_distance(center, vertex)
            if covers_more_than_one_hemisphere:
                d = np.pi - d

            if d > max_d:
                max_d = d

        # Return the min depth so that max_d > max_center_to_vertex_ipix(depth)
        depth = 0
        while max_distance_center_to_vertex(depth) >= max_d:
            depth = depth + 1

        # Get the ipixels from astropy_healpix covering the cone of (center, radius) = (center, max_d)
        lon_center, lat_center = vector.vector_to_lonlat(x=center[0], y=center[1], z=center[2], degrees=False)
        hp = HEALPix(nside=(1 << depth), order='nested', frame=ICRS())
        starting_iter_ipix = hp.cone_search_lonlat(lon=lon_center * u.rad, lat=lat_center * u.rad, radius=max_d * u.rad)

        return depth, starting_iter_ipix

    def _closes_numpy_2d_array(self, x):
        tmp = np.zeros((x.shape[0], x.shape[1] + 1))
        tmp[:, :-1] = x
        tmp[:, -1] = x[:, 0]
        return tmp

    @property
    def ipix(self):
        return self.ipix_d

    def __init__(self, ra, dec, inside=None, max_depth=10):
        ra = ra.to(u.rad).value
        dec = dec.to(u.rad).value
        # Check if the vertices form a closed polygon
        if ra[0] != ra[-1] or dec[0] != dec[-1]:
            # If not, append the first vertex to ``vertices``
            ra = np.append(ra, ra[0])
            dec = np.append(dec, dec[0])
            vertices = SkyCoord(ra=ra, dec=dec, unit="rad", frame="icrs")

        if inside:
            # Convert it to (x, y, z) cartesian coordinates on the sphere
            inside = (inside.icrs.ra.rad, inside.icrs.dec.rad)

        self.polygon = SphericalPolygon.from_lonlat(lon=ra, lat=dec, center=inside, degrees=False)

        start_depth, ipixels = self._get_starting_depth()
        end_depth = max_depth

        # When the start depth returned is > to the depth requested
        # For that specific case, we only do one iteration at start_depth
        # Thus the MOC will contain the partially intersecting cells with the
        # contained ones at start_depth

        # And we degrade the MOC to the max_depth
        self.degrade_to_max_depth = False
        if start_depth > end_depth:
            end_depth = start_depth
            self.degrade_to_max_depth = True

        self.ipix_d = {str(order): [] for order in range(start_depth, end_depth + 1)}

        ## Iterative version of the algorithm: seems a bit faster than the recursive one
        for depth in range(start_depth, end_depth + 1):
            # Define a HEALPix at the current depth
            hp = HEALPix(nside=(1 << depth), order='nested', frame=ICRS())

            # Get the lon and lat of the corners of the pixels
            # intersecting the polygon
            lon, lat = hp.boundaries_lonlat(ipixels, step=1)
            lon = lon.to(u.rad).value
            lat = lat.to(u.rad).value

            # closes the lon and lat array so that their first and last value matches
            lon = self._closes_numpy_2d_array(lon)
            lat = self._closes_numpy_2d_array(lat)

            num_ipix_inter_poly = ipixels.shape[0]

            # Define a 3d numpy array containing the corners coordinates of the intersecting pixels
            # The first dim is the num of ipixels
            # The second is the number of coordinates (5 as it defines the closed polygon of a HEALPix cell)
            # The last is of size 2 (lon and lat)
            shapes = np.vstack((lon.ravel(), lat.ravel())).T.reshape(num_ipix_inter_poly, 5, -1)

            ipix_in_polygon_l = []
            ipix_inter_polygon_l = []

            for i in range(num_ipix_inter_poly):
                shape = shapes[i]
                # Definition of a SphericalPolygon from the border coordinates of a HEALPix cell
                ipix_shape = SphericalPolygon.from_radec(lon=shape[:, 0], lat=shape[:, 1], degrees=False)
                ipix = ipixels[i]

                if self.polygon.intersects_poly(ipix_shape):
                    # If we are at the max depth then we direcly add to the MOC the intersecting ipixels
                    if depth == end_depth:
                        ipix_in_polygon_l.append(ipix)
                    else:
                        # Check whether polygon contains ipix or not
                        if self.polygon_contains_ipix(ipix_shape):
                            ipix_in_polygon_l.append(ipix)
                        else:
                            # The ipix is just intersecting without being contained in the polygon
                            # We split it in its 4 children
                            child_ipix = ipix << 2
                            ipix_inter_polygon_l.extend([child_ipix,
                                                        child_ipix + 1,
                                                        child_ipix + 2,
                                                        child_ipix + 3])

            self.ipix_d.update({str(depth): ipix_in_polygon_l})
            ipixels = np.asarray(ipix_inter_polygon_l)
