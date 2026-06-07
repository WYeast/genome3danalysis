"""Cis-ratio feature: per-bead intra-chromosomal contact fraction.

cisratio[i] = N_intra[i] / (N_intra[i] + N_inter[i])

where N_intra[i] and N_inter[i] are the number of intra- and inter-chromosomal
beads within a proximity threshold of bead i, respectively. The intra-
chromosomal definition depends on ``phased_intra``:

  - ``phased_intra=False`` (default): a neighbor is intra-chromosomal if it
    shares the same chromosome regardless of which copy. This matches the
    single-cell Hi-C convention where homologous copies are not phased, and
    is the mode used by ``intra_count_icp_cisratio_five_cells_4r_unphased``
    in the production analysis pipeline.

  - ``phased_intra=True``: a neighbor is intra-chromosomal only if it shares
    the same chromosome AND the same copy. This treats each homologous copy
    as an independent chromosome.

The kernel returns the per-bead intra/inter neighbor counts (uint32), as well
as the diploid-level ratio (float64, NaN when no neighbors are found). The
counts are required by the haploid sum-then-ratio aggregation in
``feature_extractor.run_feature``.
"""

import numpy as np
import h5py
from alabtools.utils import Index
from scipy.spatial import cKDTree

DEFAULT_DIST_CUTOFF = 240  # nm
DEFAULT_RADIUS_FACTOR = 4
DEFAULT_PHASED_INTRA = False


def _find_prox_beads(tree: cKDTree,
                     coord: np.ndarray,
                     radii: np.ndarray,
                     i: int,
                     radius_factor: float,
                     dist_sts_thresh: float,
                     r_max: float) -> np.ndarray:
    """Find proximal beads for bead i (KDTree candidate search + exact filter)."""

    if radius_factor is not None:
        # Exact threshold: d(i,j) < radius_factor * (r_i + r_j)
        # KDTree prequery upper bound (fixed j<=r_max): radius_factor * (r_i + r_max)
        r_query = float(radius_factor * (radii[i] + r_max))
        dcap_vec = radius_factor * (radii[i] + radii)
    else:
        # Exact threshold: d(i,j) < r_i + r_j + dist_sts_thresh
        r_query = float(radii[i] + r_max + dist_sts_thresh)
        dcap_vec = radii[i] + radii + dist_sts_thresh

    cand = tree.query_ball_point(coord[i], r=r_query)
    if len(cand) == 0:
        return np.empty((0,), dtype=int)

    cand = np.asarray(cand, dtype=int)
    cand = cand[cand != i]
    if cand.size == 0:
        return cand

    d = np.linalg.norm(coord[cand] - coord[i], axis=1)
    keep = d < dcap_vec[cand]
    return cand[keep]


def run(struct_id: int, hss_opt: h5py.File, params: dict) -> dict:
    """Compute per-bead intra/inter neighbor counts and the diploid cisratio.

    Args:
        struct_id: The index of the structure in the HSS file.
        hss_opt: The optimized HSS file, with coordinates of different
            structures stored in separate datasets.
        params: Feature parameters. Supported keys:
            - ``dist_cutoff`` (float, nm): surface-to-surface distance threshold
              when ``radius_factor`` is not used. Defaults to 240.
            - ``radius_factor`` (float): when set, neighbors satisfy
              ``d(i,j) < radius_factor * (r_i + r_j)``. Defaults to 4.
            - ``phased_intra`` (bool): see module docstring. Defaults to False.

    Returns:
        dict with keys:
            ``intra_count``: (nbead,) uint32, per-bead intra neighbor count.
            ``inter_count``: (nbead,) uint32, per-bead inter neighbor count.
            ``ratio``: (nbead,) float64, diploid-level cisratio (NaN if no
                       proximal neighbors).
    """

    coord = hss_opt['coordinates'][str(struct_id)][:]
    radii = hss_opt['radii'][:]
    index = Index(hss_opt)

    dist_sts_thresh = params.get('dist_cutoff', DEFAULT_DIST_CUTOFF)
    radius_factor = params.get('radius_factor', DEFAULT_RADIUS_FACTOR)
    phased_intra = params.get('phased_intra', DEFAULT_PHASED_INTRA)

    nbead = len(index)
    intra_count = np.zeros(nbead, dtype=np.uint32)
    inter_count = np.zeros(nbead, dtype=np.uint32)

    tree = cKDTree(coord)
    r_max = float(np.max(radii))

    chrom_arr = np.asarray(index.chrom)
    copy_arr = np.asarray(index.copy)

    for i in range(nbead):
        prox_beads = _find_prox_beads(tree, coord, radii, i,
                                      radius_factor, dist_sts_thresh, r_max)
        if len(prox_beads) == 0:
            continue

        chrom_prox = chrom_arr[prox_beads]
        if phased_intra:
            copy_prox = copy_arr[prox_beads]
            intra_mask = (chrom_prox == chrom_arr[i]) & (copy_prox == copy_arr[i])
        else:
            intra_mask = (chrom_prox == chrom_arr[i])

        n_intra = int(np.count_nonzero(intra_mask))
        intra_count[i] = n_intra
        inter_count[i] = len(prox_beads) - n_intra

    total = intra_count.astype(np.float64) + inter_count.astype(np.float64)
    ratio = np.full(nbead, np.nan, dtype=np.float64)
    valid = total > 0
    ratio[valid] = intra_count[valid] / total[valid]

    return {
        'intra_count': intra_count,
        'inter_count': inter_count,
        'ratio': ratio,
    }
