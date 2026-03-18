import numpy as np
import pandas as pd
from alabtools.utils import Index
import h5py
from scipy.spatial import cKDTree

DEFAULT_DIST_CUTOFF = 240  # nm
DEFAULT_RADIUS_FACTOR = 4


def _load_enhancer_counts(index: Index, enh_counts_file: str) -> np.ndarray:
    """Load enhancer counts and align them to bead order in ``index``.

    The enhancer-count BED is expected to have at least 4 columns:
      chrom, start, end, enhancer_count
    If a 5th column exists (e.g. domain label), the last column is used as count.
    Missing beads are filled with 0.
    """

    enh_df = pd.read_csv(enh_counts_file, sep='\t', header=None)
    if enh_df.shape[1] < 4:
        raise ValueError(
            "enh_counts_file must have at least 4 columns: chr, start, end, enhancer_count"
        )

    if enh_df.shape[1] == 4:
        enh_df = enh_df.iloc[:, [0, 1, 2, 3]].copy()
    else:
        enh_df = enh_df.iloc[:, [0, 1, 2, enh_df.shape[1] - 1]].copy()
    enh_df.columns = ['chr', 'start', 'end', 'enhancer_count']

    idx_df = pd.DataFrame({
        'chr': index.chromstr,
        'start': index.start,
        'end': index.end,
    })

    merged = idx_df.merge(enh_df, on=['chr', 'start', 'end'], how='left')
    counts = pd.to_numeric(merged['enhancer_count'], errors='coerce').fillna(0).to_numpy(dtype=float)
    return counts


def run(struct_id: int, hss_opt: h5py.File, params: dict) -> np.ndarray:
    """Calculate enhancer density (ENHD) for one structure.

    For each bead i, ENHD is the sum of enhancer counts carried by proximal beads
    within a distance threshold. Proximity follows the same logic as ICP:
      - if ``radius_factor`` is provided: d(i, j) <= radius_factor * (r_i + r_j)
      - otherwise: d(i, j) <= r_i + r_j + dist_cutoff

    Required params:
      - enh_counts_file: BED path with enhancer counts per bead/domain.

    Optional params:
      - radius_factor (default 4)
      - dist_cutoff (default 240; used only when radius_factor is None)
    """

    if 'enh_counts_file' not in params:
        raise KeyError("enhd requires 'enh_counts_file' in config features.enhd")

    coord = hss_opt['coordinates'][str(struct_id)][:]
    radii = hss_opt['radii'][:]
    index = Index(hss_opt)

    dist_sts_thresh = params.get('dist_cutoff', DEFAULT_DIST_CUTOFF)
    radius_factor = params.get('radius_factor', DEFAULT_RADIUS_FACTOR)

    enhancer_counts = _load_enhancer_counts(index, params['enh_counts_file'])
    tree = cKDTree(coord)
    r_max = float(np.max(radii))

    enhd = np.zeros(len(index), dtype=float)

    for i in range(len(index)):
        if radius_factor is not None:
            r_query = float(radius_factor * (radii[i] + r_max))
            dcap = radius_factor * (radii[i] + radii)
        else:
            r_query = float(radii[i] + r_max + dist_sts_thresh)
            dcap = radii[i] + radii + dist_sts_thresh

        cand = tree.query_ball_point(coord[i], r=r_query)
        if len(cand) == 0:
            continue
        cand = np.asarray(cand, dtype=int)
        cand = cand[cand != i]
        if cand.size == 0:
            continue

        d = np.linalg.norm(coord[cand] - coord[i], axis=1)
        prox_beads = cand[d < dcap[cand]]

        enhd[i] = np.sum(enhancer_counts[prox_beads])

    return enhd
