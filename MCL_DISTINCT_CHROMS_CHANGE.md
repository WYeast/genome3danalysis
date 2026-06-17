# MCL distinct-chromosome cluster filter change

## What changed

Added a new MCL config parameter:

```json
"mcl_min_distinct_chroms": 1
```

This applies to MCL-based nuclear body features:

- `speckle`
- `speckle_tsa`
- `nucleoli`
- `nucleoli_tsa`

The parameter controls how many distinct chromosome strings must be represented in each kept MCL cluster.

## Semantics

- Default: `1`
  - Preserves previous behavior.
- `2`
  - Keeps only MCL clusters containing beads from at least two distinct chromosome strings.
- `N`
  - Keeps only MCL clusters containing beads from at least N distinct chromosome strings.
- Chromosome copy is ignored:
  - `chr1` copy A and `chr1` copy B both count as `chr1`.

## Where implemented

Main implementation is in:

```text
genome3danalysis/structfeat/_mcl.py
```

New helper functions:

```python
_cluster_distinct_chrom_counts(clusters, index)
_filter_clusters_by_min_distinct_chroms(clusters, index, min_distinct_chroms)
```

The code uses:

```python
index.chromstr
```

to count distinct chromosome strings.

## Filter order

Within `compute_mcl_cluster_coms(...)`, cluster filtering now happens in this order:

1. Raw MCL clustering
2. Common size floor:
   ```python
   mcl_min_cluster_size
   ```
3. Common chromosome-diversity filter:
   ```python
   mcl_min_distinct_chroms
   ```
4. Body-specific filters:
   - speckle: `min_cluster_size`
   - nucleoli: `top_percentage`

## Cache behavior

`mcl_min_distinct_chroms` is included in the MCL cache signature. This prevents runs with different distinct-chromosome thresholds from reusing incompatible cached centroids.

## Cluster dump behavior

`mcl_clusters_out` now includes extra metadata.

For NPZ output:

- `cluster_distinct_chrom_counts`
- `mcl_min_distinct_chroms`

For HDF5 output under `/structures/{struct_id}`:

- dataset:
  - `cluster_distinct_chrom_counts`
- attr:
  - `mcl_min_distinct_chroms`

The existing `kept_mask` reflects the final kept clusters after all filters, including this new chromosome-diversity filter.

## Config usage

The parameter is feature-local. Add it only to the nuclear body feature block where the filter should apply.

Example:

```json
{
  "features": {
    "speckle": {
      "use_mcl": true,
      "regions_file": "/path/regions.bed",
      "regions_label": "SON",
      "mcl_min_distinct_chroms": 2
    }
  }
}
```

If both `speckle` and `speckle_tsa` should use the same filtered MCL bodies, put the same parameter in both blocks:

```json
{
  "features": {
    "speckle": {
      "use_mcl": true,
      "regions_file": "/path/regions.bed",
      "regions_label": "SON",
      "mcl_min_distinct_chroms": 2
    },
    "speckle_tsa": {
      "use_mcl": true,
      "regions_file": "/path/regions.bed",
      "regions_label": "SON",
      "mcl_min_distinct_chroms": 2
    }
  }
}
```

## Validation

- `mcl_min_distinct_chroms < 1` raises `ValueError`.
- `mcl_min_distinct_chroms == 1` is a no-op and preserves previous behavior.
