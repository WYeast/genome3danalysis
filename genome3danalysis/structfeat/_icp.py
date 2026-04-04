import numpy as np
from alabtools.utils import Index
import h5py
from scipy.spatial import cKDTree

DEFAULT_DIST_CUTOFF = 240  # nm
DEFAULT_RADIUS_FACTOR = 4


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
    """ Calculate the inter-chromosomal contact probability (ICP) of a structure.
    
    ICP is defined as the ratio of the number of inter-chromosomal contacts to the total number of contacts:
        ICP[i] = N_inter[i] / N_total[i],
        where N_inter[i] is the number of inter-chromosomal beads within a distance threshold of bead i,
        and N_total[i] is the total number of beads within the distance threshold of bead i.

    Args:
        struct_id (int): The index of the structure in the HSS file.
        hss_opt (h5py.File): The optimized HSS file, with coordinates of different structures in separate datasets.
        params (dict): A dictionary containing the parameters for the analysis.

    Returns:
        (np.ndarray): The inter-chromosomal contact probability of each bead in the structure.
    """
    
    # get coordinates of struct_id
    coord = hss_opt['coordinates'][str(struct_id)][:]
    
    # get the radii of the beads
    radii = hss_opt['radii'][:]
    
    # get the index
    index = Index(hss_opt)
    
    # get the surface-to-surface distance threshold
    # distance threshold can be defined either as an absolute nm cutoff (dist_cutoff)
    # or as a multiple of bead radii (radius_factor), i.e.
    # d(i,j) <= radius_factor * (r_i + r_j)
    dist_sts_thresh = params.get('dist_cutoff', DEFAULT_DIST_CUTOFF)
    radius_factor = params.get('radius_factor', DEFAULT_RADIUS_FACTOR)
    
    # initialize the inter-chromosomal contact ratio
    inter_ratio = np.zeros(len(index)).astype(float)
    
    # Build KDTree once per structure for fast neighbor queries
    tree = cKDTree(coord)
    r_max = float(np.max(radii))

    # loop over the beads to calculate the inter-chromosomal contact ratio
    for i in range(len(index)):
        prox_beads = _find_prox_beads(tree,
                                      coord,
                                      radii,
                                      i,
                                      radius_factor,
                                      dist_sts_thresh,
                                      r_max)
        
        # FILTER INTER-CHROMOSOMAL FROM PROXIMAL BEADS
        # Get the chromosome and copy of the proximal beads
        chrom_prox_beads = index.chrom[prox_beads]
        copy_prox_beads = index.copy[prox_beads]
        # Get a mask that filters only the proximal beads that are inter-chromosomal (different chromosomes or different copies)
        inter_mask = np.logical_or(chrom_prox_beads != index.chrom[i], copy_prox_beads != index.copy[i])
        # Get the proximal beads that are inter-chromosomal
        prox_inter_beads = prox_beads[inter_mask]
        
        # GET INTER-CHROMOSOMAL CONTACT RATIO
        if len(prox_beads) == 0:
            # No proximal neighbors under current threshold;
            # ICP is undefined for this bead.
            inter_ratio[i] = np.nan
        else:
            inter_ratio[i] = len(prox_inter_beads) / len(prox_beads)
        
        del prox_beads, chrom_prox_beads, copy_prox_beads, inter_mask, prox_inter_beads
    
    return inter_ratio
