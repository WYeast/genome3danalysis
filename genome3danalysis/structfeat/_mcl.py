import numpy as np
import scipy.sparse as sp
from alabtools.utils import Index
from .. import utils
from . import _radial
import os
import h5py
import fcntl
import json
import hashlib
from scipy.spatial import cKDTree

try:
    import markov_clustering as mc
except ImportError:
    mc = None

_TAGGED_BEADS_CACHE = {}
_DOMAIN_CANDIDATE_CACHE = {}


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


def _select_tagged_beads_from_radial_top(struct_id: int,
                                         hss_opt,
                                         params: dict,
                                         index: Index) -> np.ndarray:
    """Fallback selector: choose lowest-radial beads in domain regions.

    Uses:
      - params['auto_regions_from_radial_top'] as fraction in (0,1]
      - params['gap_file'] with 4th-column labels where domain/dom are candidate regions
      - params['_radial_params'] (injected from features.radial) to ensure identical radial mode
    """
    top_frac = float(params['auto_regions_from_radial_top'])
    if not (0 < top_frac <= 1):
        raise ValueError("auto_regions_from_radial_top must be in (0, 1].")

    gap_file = params.get('gap_file', None)
    if gap_file is None:
        raise KeyError("auto_regions_from_radial_top requires gap_file in config.")

    radial_params = params.get('_radial_params', None)
    if radial_params is None:
        raise KeyError("auto_regions_from_radial_top requires features.radial to be configured.")

    # Compute candidate domain beads once and reuse across structures.
    cache_key = (os.path.abspath(str(gap_file)), len(index))
    candidate_beads = _DOMAIN_CANDIDATE_CACHE.get(cache_key, None)
    if candidate_beads is None:
        _, _, _, labels = utils.read_bed(gap_file, val_type=str)
        labels = np.asarray(labels).astype(str)
        domain_tokens = {'domain', 'dom'}
        domain_hap_mask = np.array([x.strip().lower() in domain_tokens for x in labels], dtype=bool)
        domain_mtp_mask = utils.adapt_haploid_to_index(domain_hap_mask, index).astype(bool)
        candidate_beads = np.where(domain_mtp_mask)[0]
        _DOMAIN_CANDIDATE_CACHE[cache_key] = np.asarray(candidate_beads, dtype=int)
    if len(candidate_beads) == 0:
        return np.array([], dtype=int)

    radial = _radial.run(struct_id, hss_opt, radial_params)
    radial_candidates = radial[candidate_beads]
    n_pick = int(np.ceil(len(candidate_beads) * top_frac))
    n_pick = max(1, min(n_pick, len(candidate_beads)))

    # Lower radial means closer to nucleus center for both ellipsoid and
    # experimental radial definitions used in this package.
    order = np.argsort(radial_candidates)
    picked = candidate_beads[order[:n_pick]]
    return np.asarray(picked, dtype=int)


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
    sub_r = radii[tagged_beads].astype(float)
    n = sub_xyz.shape[0]

    tree = cKDTree(sub_xyz)
    r_max = float(np.max(sub_r))
    rows = []
    cols = []

    for i in range(n):
        if dist_threshold is None:
            r_query = float(2.0 * (sub_r[i] + r_max))
            cand = tree.query_ball_point(sub_xyz[i], r=r_query)
            if len(cand) == 0:
                rows.append(i)
                cols.append(i)
                continue
            cand = np.asarray(cand, dtype=int)
            d = np.linalg.norm(sub_xyz[cand] - sub_xyz[i], axis=1)
            keep = d <= (2.0 * (sub_r[i] + sub_r[cand]))
            neigh = cand[keep]
        else:
            neigh = np.asarray(tree.query_ball_point(sub_xyz[i], r=float(dist_threshold)), dtype=int)

        if neigh.size == 0:
            rows.append(i)
            cols.append(i)
        else:
            rows.extend([i] * int(neigh.size))
            cols.extend(neigh.tolist())

    data = np.ones(len(rows), dtype=np.int8)
    mat = sp.csr_matrix((data, (rows, cols)), shape=(n, n), dtype=np.int8)
    mat = mat.maximum(mat.T)
    mat.setdiag(1)
    result = mc.run_mcl(mat, inflation=inflation)
    clusters = mc.get_clusters(result)

    out = [np.array([tagged_beads[i] for i in cluster], dtype=int) for cluster in clusters]
    return out


def _cluster_distinct_chrom_counts(clusters: list, index: Index) -> np.ndarray:
    """Count distinct chromosome strings represented in each cluster.

    Chromosome copy is intentionally ignored: beads from chr1 copy A and
    chr1 copy B both count as chr1.
    """
    if len(clusters) == 0:
        return np.empty((0,), dtype=np.int32)

    counts = []
    for c in clusters:
        chroms = np.asarray(index.chromstr[np.asarray(c, dtype=int)]).astype(str)
        counts.append(len(np.unique(chroms)))
    return np.asarray(counts, dtype=np.int32)


def _filter_clusters_by_min_distinct_chroms(clusters: list,
                                            index: Index,
                                            min_distinct_chroms: int) -> list:
    """Keep clusters represented by at least N distinct chromosome strings."""
    min_distinct_chroms = int(min_distinct_chroms)
    if min_distinct_chroms < 1:
        raise ValueError("mcl_min_distinct_chroms must be >= 1.")
    if min_distinct_chroms == 1 or len(clusters) == 0:
        return clusters

    counts = _cluster_distinct_chrom_counts(clusters, index)
    return [c for c, n_chrom in zip(clusters, counts) if n_chrom >= min_distinct_chroms]


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
                       kept_centroids: np.ndarray,
                       cluster_distinct_chrom_counts: np.ndarray,
                       mcl_min_distinct_chroms: int) -> str:
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
        cluster_distinct_chrom_counts=np.asarray(cluster_distinct_chrom_counts, dtype=np.int32),
        mcl_min_distinct_chroms=np.int32(mcl_min_distinct_chroms),
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
                          kept_centroids: np.ndarray,
                          cluster_distinct_chrom_counts: np.ndarray,
                          mcl_min_distinct_chroms: int) -> str:
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
            grp.create_dataset('cluster_distinct_chrom_counts',
                               data=np.asarray(cluster_distinct_chrom_counts, dtype=np.int32),
                               dtype=np.int32)
            grp.create_dataset('kept_mask', data=np.asarray(kept_mask, dtype=bool), dtype=bool)
            grp.create_dataset('all_cluster_centroids', data=np.asarray(all_centroids, dtype=np.float32), dtype=np.float32)
            grp.create_dataset('kept_cluster_centroids', data=np.asarray(kept_centroids, dtype=np.float32), dtype=np.float32)

            grp.attrs['struct_id'] = int(struct_id)
            grp.attrs['body_type'] = str(body_type)
            grp.attrs['n_clusters_all'] = int(len(all_clusters))
            grp.attrs['n_clusters_kept'] = int(np.sum(np.asarray(kept_mask, dtype=bool)))
            grp.attrs['mcl_min_distinct_chroms'] = int(mcl_min_distinct_chroms)

        fcntl.flock(lock_fp.fileno(), fcntl.LOCK_UN)
    return out_h5


def compute_mcl_cluster_coms(struct_id: int, hss_opt, params: dict, body_type: str) -> np.ndarray:
    index = Index(hss_opt)
    coords = hss_opt['coordinates'][str(struct_id)][:]
    radii = hss_opt['radii'][:]

    has_region_tags = ('regions_file' in params and 'regions_label' in params)
    if has_region_tags:
        tag_bed = params['regions_file']
        tag_name = params['regions_label']
        cache_key = (os.path.abspath(str(tag_bed)), str(tag_name), len(index))
        tagged_beads = _TAGGED_BEADS_CACHE.get(cache_key, None)
        if tagged_beads is None:
            hap_labels = _load_haploid_labels_from_bed(tag_bed, index)
            labels_mtp = utils.adapt_haploid_to_index(hap_labels, index)
            tagged_beads = np.where(labels_mtp == tag_name)[0]
            _TAGGED_BEADS_CACHE[cache_key] = np.asarray(tagged_beads, dtype=int)
    elif 'auto_regions_from_radial_top' in params:
        tagged_beads = _select_tagged_beads_from_radial_top(struct_id, hss_opt, params, index)
    else:
        raise KeyError("MCL requires either regions_file+regions_label or auto_regions_from_radial_top.")

    dist_thr = params.get('mcl_dist_threshold', None)
    mcl_inflation = float(params.get('mcl_inflation', 1.4))
    mcl_min_cluster_size = int(params.get('mcl_min_cluster_size', 2))
    mcl_min_distinct_chroms = int(params.get('mcl_min_distinct_chroms', 1))
    if mcl_min_distinct_chroms < 1:
        raise ValueError("mcl_min_distinct_chroms must be >= 1.")

    # Optional cache: reuse MCL centroids across features (e.g. speckle + speckle_tsa).
    cache_file = params.get('_mcl_cache_h5', None)
    cache_signature = None
    if cache_file is not None:
        sig_obj = {
            'body_type': body_type,
            'regions_file': params.get('regions_file', None),
            'regions_label': params.get('regions_label', None),
            'auto_regions_from_radial_top': params.get('auto_regions_from_radial_top', None),
            'gap_file': params.get('gap_file', None),
            'radial_params': params.get('_radial_params', None),
            'mcl_dist_threshold': params.get('mcl_dist_threshold', None),
            'mcl_inflation': params.get('mcl_inflation', 1.4),
            'mcl_min_cluster_size': params.get('mcl_min_cluster_size', 2),
            'mcl_min_distinct_chroms': params.get('mcl_min_distinct_chroms', 1),
            'min_cluster_size': params.get('min_cluster_size', 5),
            'top_percentage': params.get('top_percentage', 95.0)
        }
        sig_str = json.dumps(sig_obj, sort_keys=True, default=str)
        cache_signature = hashlib.md5(sig_str.encode('utf-8')).hexdigest()
        cache_file = os.path.abspath(str(cache_file))
        if os.path.exists(cache_file):
            lock_path = cache_file + '.lock'
            with open(lock_path, 'a+') as lock_fp:
                fcntl.flock(lock_fp.fileno(), fcntl.LOCK_SH)
                with h5py.File(cache_file, 'r') as h5:
                    path = f"mcl_centroids/{cache_signature}/{struct_id}"
                    if path in h5:
                        out = h5[path][:]
                        fcntl.flock(lock_fp.fileno(), fcntl.LOCK_UN)
                        return out
                fcntl.flock(lock_fp.fileno(), fcntl.LOCK_UN)

    all_clusters = _cluster_beads_by_mcl(coords,
                                         radii,
                                         tagged_beads,
                                         dist_thr,
                                         mcl_inflation)

    # First-stage common MCL filtering
    clusters = [c for c in all_clusters if len(c) >= mcl_min_cluster_size]
    clusters = _filter_clusters_by_min_distinct_chroms(clusters, index, mcl_min_distinct_chroms)

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

    cluster_distinct_chrom_counts = _cluster_distinct_chrom_counts(all_clusters, index)

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
                                  kept_centroids,
                                  cluster_distinct_chrom_counts,
                                  mcl_min_distinct_chroms)
        else:
            _save_cluster_dump(struct_id,
                               out_base,
                               body_type,
                               coords,
                               all_clusters,
                               kept_mask,
                               all_centroids,
                               kept_centroids,
                               cluster_distinct_chrom_counts,
                               mcl_min_distinct_chroms)

    if cache_file is not None and cache_signature is not None:
        cache_dir = os.path.dirname(cache_file)
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)
        lock_path = cache_file + '.lock'
        with open(lock_path, 'w') as lock_fp:
            fcntl.flock(lock_fp.fileno(), fcntl.LOCK_EX)
            with h5py.File(cache_file, 'a') as h5:
                grp = h5.require_group(f"mcl_centroids/{cache_signature}")
                ds = str(struct_id)
                if ds in grp:
                    del grp[ds]
                grp.create_dataset(ds, data=np.asarray(kept_centroids, dtype=np.float32), dtype=np.float32)
            fcntl.flock(lock_fp.fileno(), fcntl.LOCK_UN)

    return kept_centroids
