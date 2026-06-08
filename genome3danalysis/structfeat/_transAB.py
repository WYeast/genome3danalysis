import numpy as np
from alabtools.utils import Index
import h5py
from scipy.spatial import cKDTree
from .. import utils

DEFAULT_DIST_CUTOFF = 500  # nm
DEFAULT_RADIUS_FACTOR = None


def _find_prox_beads(tree: cKDTree,
                     coord: np.ndarray,
                     radii: np.ndarray,
                     i: int,
                     radius_factor: float,
                     dist_sts_thresh: float,
                     r_max: float) -> np.ndarray:
    """Find proximal beads for one bead with KDTree candidate search + exact filter."""

    if radius_factor is not None:
        # Exact threshold: d(i,j) < radius_factor * (r_i + r_j)
        # KDTree prequery upper bound (fixed j<=r_max): radius_factor * (r_i + r_max)
        r_query = float(radius_factor * (radii[i] + r_max))
        dcap_vec = radius_factor * (radii[i] + radii)
    else:
        # Exact threshold: d(i,j) < r_i + r_j + dist_sts_thresh
        # KDTree prequery upper bound: r_i + r_max + dist_sts_thresh
        r_query = float(radii[i] + r_max + dist_sts_thresh)
        dcap_vec = radii[i] + radii + dist_sts_thresh

    cand = tree.query_ball_point(coord[i], r=r_query)
    if len(cand) == 0:
        return np.empty((0,), dtype=int)

    cand = np.asarray(cand, dtype=int)
    cand = cand[cand != i]
    if cand.size == 0:
        return cand

    # Exact threshold check (strictly '<' to match previous behavior)
    d = np.linalg.norm(coord[cand] - coord[i], axis=1)
    keep = d < dcap_vec[cand]
    return cand[keep]

def run(struct_id: int, hss_opt: h5py.File, params: dict) -> np.ndarray:
    """ Calculate the A/B ratio for each bead in the structure.

    By default (``include_intra=False``, the original ``transAB`` semantics):
        transAB[i] = N_trans_A[i] / (N_trans_A[i] + N_trans_B[i])
    where N_trans_{A,B}[i] count *inter-chromosomal* proximal neighbors of
    bead i with the corresponding A/B label.

    When ``include_intra=True`` (the ``abratio`` semantics, all neighbors):
        abratio[i] = N_all_A[i] / (N_all_A[i] + N_all_B[i])
    where N_all_{A,B}[i] count *all* proximal neighbors of bead i regardless
    of chromosome / copy. The intra-chromosomal A/B-labelled neighbors are
    pooled in alongside the trans ones; copies of bead i itself are still
    excluded via the KDTree neighbor selection.

    Args:
        struct_id (int): The index of the structure in the HSS file.
        hss_opt (h5py.File): The optimized HSS file, with coordinates of
            different structures in separate datasets.
        params (dict): A dictionary of parameters. Supported keys:
            - ``filename`` (str): A/B compartment BED (required).
            - ``radius_factor`` (float): proximity threshold as multiple of
              bead radii. If provided, takes priority over ``dist_cutoff``.
            - ``dist_cutoff`` (float, nm): surface-to-surface distance
              threshold when ``radius_factor`` is not used. Defaults to 500.
            - ``include_intra`` (bool): when True, count all A/B-labelled
              proximal neighbors (intra + inter). Defaults to False
              (backward-compatible trans-only behavior).

    Returns:
        (np.ndarray): per-bead A/B ratio (n_intra+n_trans selection
            depending on ``include_intra``).
    """
    
    # get the index
    index = Index(hss_opt)
    
    # Read A/B compartment file name
    try:
        filename = params['filename']
    except KeyError:
        raise KeyError('AB-compartment filename must be specified')
    # Load the A/B compartment file with pandas
    try:
        _, _, _, ab = utils.read_bed(filename, val_type=str)
    except KeyError:
        raise KeyError('AB-compartment file not found')
    # Assert that AB track is correct
    assert isinstance(ab, np.ndarray), 'AB-compartment track must be a numpy array'
    assert len(ab) == np.sum(index.copy == 0),\
        'AB-compartment track must have the same length as the number of haploid beads in the structure'
    assert len(np.unique(ab)) == 2 or len(np.unique(ab)) == 3,\
        'AB-compartment track must contain only A and B and optionally NA (unspecified format)'
    assert 'A' in np.unique(ab) and 'B' in np.unique(ab), 'AB-compartment track must contain A and B'    
    # Adapt AB track to multi-ploid index
    ab = utils.adapt_haploid_to_index(ab, index)
    
    # Distance threshold can be defined either as an absolute nm cutoff
    # (dist_cutoff) or as a multiple of bead radii (radius_factor), i.e.
    # d(i,j) <= radius_factor * (r_i + r_j)
    # Priority: radius_factor if provided; otherwise dist_cutoff.
    dist_sts_thresh = params.get('dist_cutoff', DEFAULT_DIST_CUTOFF)
    radius_factor = params.get('radius_factor', DEFAULT_RADIUS_FACTOR)

    # include_intra controls whether intra-chromosomal proximal beads are
    # also counted toward N_A / N_B. False (default) preserves the original
    # trans-only ``transAB`` semantics; True yields the all-neighbor
    # ``abratio`` semantics.
    include_intra = bool(params.get('include_intra', False))

    # get coordinates of struct_id
    coord = hss_opt['coordinates'][str(struct_id)][:]

    # get the radii of the beads
    radii = hss_opt['radii'][:]

    # Initialize the ratio array
    ab_ratio = np.zeros(len(index)).astype(float)

    # Build KDTree once per structure for fast neighbor queries
    tree = cKDTree(coord)
    r_max = float(np.max(radii))

    # Loop over all beads
    for i in range(len(index)):

        # FIND PROXIMAL BEADS
        prox_beads = _find_prox_beads(tree,
                                      coord,
                                      radii,
                                      i,
                                      radius_factor,
                                      dist_sts_thresh,
                                      r_max)

        if include_intra:
            # All neighbors counted; no chromosome filter.
            prox_ab_beads = prox_beads
        else:
            # Filter to inter-chromosomal neighbors only (different chrom OR copy).
            chrom_prox_beads = index.chrom[prox_beads]
            copy_prox_beads = index.copy[prox_beads]
            inter_mask = np.logical_or(chrom_prox_beads != index.chrom[i],
                                       copy_prox_beads != index.copy[i])
            prox_ab_beads = prox_beads[inter_mask]
            del chrom_prox_beads, copy_prox_beads, inter_mask

        # COMPUTE A/B RATIO over the selected neighbors
        ab_prox = ab[prox_ab_beads]
        n_A = np.sum(ab_prox == 'A')
        n_B = np.sum(ab_prox == 'B')
        if (n_A + n_B) == 0:
            # No A/B-labelled neighbors under current threshold/selection.
            ab_ratio[i] = np.nan
        else:
            ab_ratio[i] = n_A / (n_A + n_B)

        del prox_beads, prox_ab_beads, ab_prox

    return ab_ratio
