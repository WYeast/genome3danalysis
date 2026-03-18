import numpy as np
import scipy.sparse as sp
from alabtools.utils import Index
from .. import utils
import os
import h5py
import fcntl

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
                          inflation: float) -> list:
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

    out = [np.array([tagged_beads[i] for i in cluster], dtype=int) for cluster in clusters]
    return out


def _clusters_to_centroids(coords: np.ndarray, radii: np.ndarray, clusters: list) -> np.ndarray:
    """Compute size-weighted cluster centroids.

    Weights are proportional to bead volume (r^3), consistent with current code.
    """
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


def _resolve_struct_outfile(out_base: str, struct_id: int) -> str:
    """Resolve output filename for per-structure cluster dump.

    If out_base contains '{struct_id}', format it directly.
    Otherwise write to '<root>_struct<id><ext>'.
    """
    out_base = os.path.abspath(out_base)
    if '{struct_id}' in out_base:
        return out_base.format(struct_id=struct_id)
    root, ext = os.path.splitext(out_base)
    if ext == '':
        ext = '.npz'
    return f"{root}_struct{struct_id}{ext}"


def _save_cluster_dump(struct_id: int,
                       out_base: str,
                       body_type: str,
                       coords: np.ndarray,
                       all_clusters: list,
                       kept_mask: np.ndarray,
                       all_centroids: np.ndarray,
                       kept_centroids: np.ndarray) -> str:
    """Save all MCL clusters for one structure to compressed NPZ.

    Stored content includes all clusters before body-specific filtering.
    """
    out_file = _resolve_struct_outfile(out_base, struct_id)
    out_dir = os.path.dirname(out_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Ragged per-cluster arrays saved as object arrays for compactness and easy loading.
    all_idx_obj = np.array([np.asarray(c, dtype=np.int32) for c in all_clusters], dtype=object)
    all_xyz_obj = np.array([coords[c, :].astype(np.float32) for c in all_clusters], dtype=object)
    sizes = np.array([len(c) for c in all_clusters], dtype=np.int32)

    np.savez_compressed(
        out_file,
        struct_id=np.int32(struct_id),
        body_type=np.array(body_type),
        cluster_bead_indices=all_idx_obj,
        cluster_bead_coords=all_xyz_obj,
        cluster_sizes=sizes,
        kept_mask=np.asarray(kept_mask, dtype=bool),
        all_cluster_centroids=np.asarray(all_centroids, dtype=np.float32),
        kept_cluster_centroids=np.asarray(kept_centroids, dtype=np.float32),
    )
    return out_file


def _pack_ragged_indices(clusters: list) -> tuple:
    """Pack list of 1D index arrays to flat array + offsets."""
    if len(clusters) == 0:
        return np.empty((0,), dtype=np.int32), np.zeros((1,), dtype=np.int64)
    sizes = np.array([len(c) for c in clusters], dtype=np.int64)
    offsets = np.concatenate([[0], np.cumsum(sizes)])
    flat = np.concatenate([np.asarray(c, dtype=np.int32) for c in clusters], axis=0)
    return flat, offsets


def _pack_ragged_coords(coords: np.ndarray, clusters: list) -> tuple:
    """Pack list of cluster xyz matrices to flat xyz array + offsets."""
    if len(clusters) == 0:
        return np.empty((0, 3), dtype=np.float32), np.zeros((1,), dtype=np.int64)
    mats = [coords[c, :].astype(np.float32) for c in clusters]
    sizes = np.array([m.shape[0] for m in mats], dtype=np.int64)
    offsets = np.concatenate([[0], np.cumsum(sizes)])
    flat = np.concatenate(mats, axis=0)
    return flat, offsets


def _save_cluster_dump_h5(struct_id: int,
                          out_h5: str,
                          body_type: str,
                          coords: np.ndarray,
                          all_clusters: list,
                          kept_mask: np.ndarray,
                          all_centroids: np.ndarray,
                          kept_centroids: np.ndarray) -> str:
    """Append/overwrite one structure's all-cluster dump into a single HDF5 file."""
    out_h5 = os.path.abspath(out_h5)
    out_dir = os.path.dirname(out_h5)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Lock file to serialize concurrent writers from parallel workers.
    lock_path = out_h5 + '.lock'
    with open(lock_path, 'w') as lock_fp:
        fcntl.flock(lock_fp.fileno(), fcntl.LOCK_EX)
        with h5py.File(out_h5, 'a') as h5:
            grp_root = h5.require_group('structures')
            sname = str(struct_id)
            if sname in grp_root:
                del grp_root[sname]
            grp = grp_root.create_group(sname)

            idx_flat, idx_offsets = _pack_ragged_indices(all_clusters)
            xyz_flat, xyz_offsets = _pack_ragged_coords(coords, all_clusters)
            sizes = np.array([len(c) for c in all_clusters], dtype=np.int32)

            grp.create_dataset('cluster_bead_indices_flat', data=idx_flat, dtype=idx_flat.dtype)
            grp.create_dataset('cluster_bead_indices_offsets', data=idx_offsets, dtype=idx_offsets.dtype)
            grp.create_dataset('cluster_bead_coords_flat', data=xyz_flat, dtype=xyz_flat.dtype)
            grp.create_dataset('cluster_bead_coords_offsets', data=xyz_offsets, dtype=xyz_offsets.dtype)
            grp.create_dataset('cluster_sizes', data=sizes, dtype=sizes.dtype)
            grp.create_dataset('kept_mask', data=np.asarray(kept_mask, dtype=bool), dtype=bool)
            grp.create_dataset('all_cluster_centroids', data=np.asarray(all_centroids, dtype=np.float32), dtype=np.float32)
            grp.create_dataset('kept_cluster_centroids', data=np.asarray(kept_centroids, dtype=np.float32), dtype=np.float32)

            grp.attrs['struct_id'] = int(struct_id)
            grp.attrs['body_type'] = str(body_type)
            grp.attrs['n_clusters_all'] = int(len(all_clusters))
            grp.attrs['n_clusters_kept'] = int(np.sum(np.asarray(kept_mask, dtype=bool)))

        fcntl.flock(lock_fp.fileno(), fcntl.LOCK_UN)
    return out_h5


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

    all_clusters = _cluster_beads_by_mcl(coords,
                                         radii,
                                         tagged_beads,
                                         dist_thr,
                                         mcl_inflation)

    # First-stage common MCL filtering
    clusters = [c for c in all_clusters if len(c) >= mcl_min_cluster_size]

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

    # Compute centroids for all and for kept clusters
    all_centroids = _clusters_to_centroids(coords, radii, all_clusters)
    kept_centroids = _clusters_to_centroids(coords, radii, clusters)

    # Optional cluster dump output.
    # This dump always includes ALL clusters (before size filtering),
    # plus kept_mask and kept centroids to support downstream selection.
    out_base = params.get('mcl_clusters_out', None)
    if out_base is not None:
        kept_set = set(tuple(c.tolist()) for c in clusters)
        kept_mask = np.array([tuple(c.tolist()) in kept_set for c in all_clusters], dtype=bool)
        if str(out_base).lower().endswith('.h5'):
            _save_cluster_dump_h5(struct_id,
                                  out_base,
                                  body_type,
                                  coords,
                                  all_clusters,
                                  kept_mask,
                                  all_centroids,
                                  kept_centroids)
        else:
            _save_cluster_dump(struct_id,
                               out_base,
                               body_type,
                               coords,
                               all_clusters,
                               kept_mask,
                               all_centroids,
                               kept_centroids)

    return kept_centroids
