# 3dnalysis
Package for analyzing 3D genome features.

## New configuration options

### 1) ICP supports radius-factor threshold

`icp` can now use either:

- `dist_cutoff`: absolute surface-to-surface cutoff in nm (legacy), or
- `radius_factor`: center-to-center threshold as `N * (r_i + r_j)`.

Example:

```json
{
  "features": {
    "icp": {
      "radius_factor": 1.5
    }
  }
}
```

### 2) speckle/nucleoli with built-in MCL clustering

For `speckle`, `speckle_tsa`, `nucleoli`, `nucleoli_tsa`, you can set `use_mcl=true`
to compute body clusters on the fly from tagged regions.

Common MCL parameters:

- `regions_file`: BED with region labels in 4th column
- `regions_label`: label name to select seed regions (e.g. `SON`, `POL1RE`)
- `mcl_dist_threshold`: optional absolute center distance threshold for graph edges
  (if omitted, uses `2*(r_i+r_j)`)
- `mcl_inflation`: MCL inflation, default `1.4`
- `mcl_min_cluster_size`: MCL cluster floor, default `2`

Body-specific cluster filtering:

- speckle: `min_cluster_size` (default `5`)
- nucleoli: `top_percentage` (default `95.0`), keep clusters with size >= percentile threshold

### 3) Experimental nucleus shape for lamina/lamina_tsa

`lamina` and `lamina_tsa` now support `shape="experimental"`, using per-structure
shape maps from `.bin` or `.mrc`.

Provide one of:

- `exp_shape_file`: template path, e.g. `"/path/lamina/{struct_id}.bin"`
- `exp_shape_dir` + `exp_shape_ext` (default `.bin`)

Optional for `.mrc`:

- `transpose_mode`: `none` (default), `T`, or `swap_xz`
- `voxel_size`: override voxel size `[dx, dy, dz]`
- `exp_shape_threshold`: binarization threshold (default `0.5`)

For experimental mode:

- `lamina` outputs LAD as minimum distance to lamina surface (distance-transform approximation)
- `lamina_tsa` outputs `exp(-tsa_exponent * LAD)`

### 4) gap_file supports boolean or cen/domain-style labels

`gap_file` must be a 4-column BED (no header). The 4th column can be either:

- boolean-like: `True/False`, `1/0`, `yes/no`;
- label-style: `domain`/`dom` are treated as non-gap, all other labels
  (e.g. `cen`, `tel`, `gap`, `contig`) are treated as gap.
