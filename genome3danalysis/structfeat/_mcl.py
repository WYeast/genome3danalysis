import numpy as np
import scipy.sparse as sp
from alabtools.utils import Index
from .. import utils

try:
    import markov_clustering as mc
except ImportError:
    mc = None


def _load_haploid_labels_from_bed(filename: str, index: Index, default: str = 'NA') -> np.ndarray:
    chroms, starts, ends, labels = utils.read_bed(filename, val_type=str)
    label_map = {(str(c), int(s), int(e)): str(v) for c, s, e, v in zip(chroms, starts, ends, labels)}

    hap_idx = np.where(index.copy == 0)[0]
    hap_labels = np.full(len(hap_idx), default, dtype=object)
    for i, bi in enumerate(hap_idx):
        key = (str(index.chromstr[bi]), int(index.start[bi]), int(index.end[bi]))
        if key in label_map:
            hap_labels[i] = label_map[key]
    return hap_labels


def _cluster_beads_by_mcl(coords: np.ndarray,
                          radii: np.ndarray,
                          tagged_beads: np.ndarray,
                          dist_threshold: float,
                          inflation: float,
                          min_cluster_size: int) -> list:
    if len(tagged_beads) == 0:
        return []
    if mc is None:
        raise ImportError('markov_clustering is required for MCL clustering features.')

    sub_xyz = coords[tagged_beads, :]
    dists = np.linalg.norm(sub_xyz[:, None, :] - sub_xyz[None, :, :], axis=2)

    if dist_threshold is None:
        dcap = 2.0 * (radii[tagged_beads][:, None] + radii[tagged_beads][None, :])
        adj = dists <= dcap
    else:
        adj = dists <= dist_threshold

    mat = sp.csr_matrix(adj.astype(np.int8))
    result = mc.run_mcl(mat, inflation=inflation)
    clusters = mc.get_clusters(result)

    out = []
    for cluster in clusters:
        if len(cluster) >= min_cluster_size:
            out.append(np.array([tagged_beads[i] for i in cluster], dtype=int))
    return out


def compute_mcl_cluster_coms(struct_id: int, hss_opt, params: dict, body_type: str) -> np.ndarray:
    index = Index(hss_opt)
    coords = hss_opt['coordinates'][str(struct_id)][:]
    radii = hss_opt['radii'][:]

    tag_bed = params['regions_file']
    tag_name = params['regions_label']

    hap_labels = _load_haploid_labels_from_bed(tag_bed, index)
    labels_mtp = utils.adapt_haploid_to_index(hap_labels, index)
    tagged_beads = np.where(labels_mtp == tag_name)[0]

    dist_thr = params.get('mcl_dist_threshold', None)
    mcl_inflation = float(params.get('mcl_inflation', 1.4))
    mcl_min_cluster_size = int(params.get('mcl_min_cluster_size', 2))

    clusters = _cluster_beads_by_mcl(coords,
                                     radii,
                                     tagged_beads,
                                     dist_thr,
                                     mcl_inflation,
                                     mcl_min_cluster_size)

    if body_type == 'speckle':
        min_size = int(params.get('min_cluster_size', 5))
        clusters = [c for c in clusters if len(c) >= min_size]
    elif body_type == 'nucleoli':
        top_percentage = float(params.get('top_percentage', 95.0))
        if len(clusters) > 0:
            sizes = np.array([len(c) for c in clusters], dtype=float)
            thr = int(np.ceil(np.percentile(sizes, top_percentage)))
            thr = max(thr, 1)
            clusters = [c for c in clusters if len(c) >= thr]

    if len(clusters) == 0:
        return np.empty((0, 3))

    coms = []
    for c in clusters:
        xyz = coords[c, :]
        w = radii[c].astype(float) ** 3
        if np.sum(w) == 0:
            coms.append(np.mean(xyz, axis=0))
        else:
            coms.append(np.sum(xyz * w[:, None], axis=0) / np.sum(w))

    return np.array(coms)