import numpy as np
# A python module handling graph manipulations
import networkx as nx

from .utils import make_wcs

from astropy_healpix import HEALPix

from astropy.coordinates import ICRS, SkyCoord
from astropy.wcs.utils import skycoord_to_pixel

class Boundaries():
    @staticmethod
    def get(moc, order):
        boundaries_l = []

        # Get the ipixels of the MOC at the deepest order
        hp, ipixels = Boundaries._compute_HEALPix_indices(moc, order)
        # Compute a graph of the MOC where each node is an ipix belonging to the MOC
        G = Boundaries._build_graph_from_ipixels(hp, ipixels)

        # Split the global MOC graph into all its non connected subgraphs.
        G_subgraphs = nx.connected_components(G)
        for g in G_subgraphs:
            graph_boundaries = Boundaries._compute_graph_HEALPix_boundaries(hp, g)
            boundaries_l.extend(Boundaries._retrieve_skycoords(graph_boundaries))

        return boundaries_l

    @staticmethod
    def _compute_HEALPix_indices(m, order):
        moc = m
        if order:
            if m.max_order > order:
                moc = m.degrade_to_order(order)

        max_order = moc.max_order
        hp = HEALPix(nside=(1 << max_order), order='nested', frame=ICRS())
        ipixels = moc._best_res_pixels()

        # Take the complement if the MOC covers more than half of the sky => the perimeter(MOC) = perimeter(complement(MOC))
        # but we process less hpx cells
        num_ipixels = 3 << (2*(max_order + 1))
        sky_fraction = ipixels.shape[0] / float(num_ipixels)

        if sky_fraction > 0.5:
            ipixels_all = np.arange(num_ipixels)
            ipixels = np.setdiff1d(ipixels_all, ipixels, assume_unique=True)

        return hp, ipixels

    @staticmethod
    # Faster version for computing the graph of all the connected ipixels
    def _build_graph_from_ipixels(hp, ipixels, dir_connected=[0, 2, 4, 6]):
        """
        Build a graph from a list of ipixels

        The connexion relation between a node and its neighbours can be specified.
        By default, each ipix is only connected to its direct neighbours (e.g.
        west, south, east and north)
        """
        neighbours = hp.neighbours(ipixels)
        # Select only the WEST, SOUTH, EAST and NORTH neighbours (i.e. the direct ones)
        neighbours = neighbours[dir_connected, :]

        # Select only neighbours lying in the ipixels ensemble
        mask = np.isin(neighbours, ipixels)
        edges = np.array([])
        i = 0
        for k in dir_connected:
            mask_k = mask[i, :]
            new_edges = np.vstack((ipixels, neighbours[i, :]))[:, mask_k]
            if edges.size == 0:
                edges = new_edges
                alone_id = ~mask_k
            else:
                edges = np.hstack((edges, new_edges))
                alone_id &= ~mask_k
            i += 1

        edges = edges.T
        ipix_alone = ipixels[alone_id]

        # Graph instanciation
        G = nx.Graph()
        # Add the edges giving the interaction between all the ipixel nodes
        G.add_edges_from(edges)
        # Add nodes connected to nothing
        G.add_nodes_from(ipix_alone)
        return G

    @staticmethod
    def _compute_graph_HEALPix_boundaries(hp, g):
        def insert_edge(G, l1, l2, p1, p2):
            # Nodes are indexed by str(skycoord). When getting ordered nodes, one can retrieve back the skycoord instance
            # by accessing the python dict `pts_d`.
            try:
                # Avoid the special case where holes are touching to each other
                # 'x' belongs to the MOC
                # ' ' is part of the holes in the MOC
                #    |xxx
                #    |xxx
                # ---A---
                # xxx|
                # xxx|
                #
                # If this case occurs we split the node A into 2. One is attached to the bottom left graph and the other to the
                # top right one. When computing the MST (minimal spanning tree) from a graph, we need our graphs to have
                # only nodes of degrees 1 or 2 (i.e. to be lines).
                if G.degree[l1] >= 2:
                    l1 += '_'
            except:
                pass

            try:
                if G.degree[l2] >= 2:
                    l2 += '_'
            except:
                pass
            # Set the skycoord instance as an attribute of the nodes
            G.add_node(l1, ra=p1.ra.deg, dec=p1.dec.deg)
            G.add_node(l2, ra=p2.ra.deg, dec=p2.dec.deg)
            G.add_edge(l1, l2)

        # Phase 1: Retrieve the ipixels located at the border of
        # this connexe MOC component
        ipixels = np.asarray(list(g))

        neighbours = hp.neighbours(ipixels)[[0, 2, 4, 6], :]
        isin = np.isin(neighbours, ipixels)
        border = isin.sum(axis=0) < 4

        ipixels_border = ipixels[border]
        isin_border = isin[:, border]

        # Phase 2: Build the graph from the positions of the ipixels boundaries
        ipixels_border_bounds = hp.boundaries_skycoord(ipixels_border, step=1)
        V = nx.Graph()
        for i in range(ipixels_border.shape[0]):
            ipix = ipixels_border[i]
            ipix_bound = ipixels_border_bounds[i]

            p0 = ipix_bound[0]
            p1 = ipix_bound[1]
            p2 = ipix_bound[2]
            p3 = ipix_bound[3]
            s0 = str(p0)
            s1 = str(p1)
            s2 = str(p2)
            s3 = str(p3)

            # WEST border
            if not isin_border[0, i]:
                insert_edge(V, s1, s2, p1, p2)

            # NORTH border
            if not isin_border[3, i]:
                insert_edge(V, s2, s3, p2, p3)

            # EAST border
            if not isin_border[2, i]:
                insert_edge(V, s3, s0, p3, p0)

            # SOUTH border
            if not isin_border[1, i]:
                insert_edge(V, s0, s1, p0, p1)

        return V

    @staticmethod
    def _retrieve_skycoords(V):
        coords_l = []
        # Accessing the borders one by one. At this step, V_subgraphs contains a list of cycles
        # (i.e. one describing the external border of the MOC component and several describing the holes
        # found in the MOC component).
        V_subgraphs = nx.connected_component_subgraphs(V)
        for v in V_subgraphs:
            # Compute the MST for each cycle
            v = nx.convert_node_labels_to_integers(v)
            mst = nx.minimum_spanning_tree(v)
            # Get one end of the span tree by looping over its node and checking if the degree is one
            src = None
            for (node, deg) in mst.degree():
                if deg == 1:
                    src = node
                    break

            # Get the unordered lon and lat
            ra = np.asarray(list(nx.get_node_attributes(v, 'ra').values()))
            dec = np.asarray(list(nx.get_node_attributes(v, 'dec').values()))
            coords = np.vstack((ra, dec)).T
            # Get the ordering from the MST
            ordering = np.asarray(list(nx.dfs_preorder_nodes(mst, src)))
            # Order the coords
            coords = coords[ordering]
            # Get a skycoord containing N coordinates computed from the Nx2 `coords` array
            coords = SkyCoord(coords, unit="deg")
            coords_l.append(coords)

        return coords_l