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

`lamina`, `lamina_tsa`, and `radial` now support `shape="experimental"`, using per-structure
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
- `radial` outputs `1 - clip(LAD/R, 0, 1)` per structure, with `R` controlled by `exp_radial_mode`

For shape path resolution, both flat and nested layouts are supported:

- flat: `/path/lamina_bin/{struct_id}.bin`
- nested: `/path/lamina_bin/{struct_id}/lamina/lamina.bin`


`radial` experimental-mode options:

- `exp_radial_mode`: `max_lad` (default, `R=max(LAD)`) or `half_longest_axis` (`R=0.5*longest_axis`)

### 4) gap_file supports boolean or cen/domain-style labels

`gap_file` must be a 4-column BED (no header). The 4th column can be either:

- boolean-like: `True/False`, `1/0`, `yes/no`;
- label-style: `domain`/`dom` are treated as non-gap, all other labels
  (e.g. `cen`, `tel`, `gap`, `contig`) are treated as gap.

### 5) Control output `.sf.h5` filename and folder

By default, `structfeat-run` writes output as:

- `<hss_name with .hss replaced by .sf.h5>`

You can override this in config:

- `sf_name`: explicit output file path (highest priority)
- `output_dir`: output folder used with auto-generated `<basename(hss_name)>.sf.h5`

Example:

```json
{
  "hss_name": "/path/to/input/model.hss",
  "output_dir": ".",
  "features": {
    "radial": {
      "shape": "experimental",
      "exp_shape_file": "/path/to/lamina_bin/{struct_id}.bin"
    }
  }
}
```

or

```json
{
  "hss_name": "/path/to/input/model.hss",
  "sf_name": "./my_custom_output.sf.h5",
  "features": {
    "lamina": {
      "shape": "experimental",
      "exp_shape_file": "/path/to/lamina_bin/{struct_id}.bin"
    }
  }
}
```

## Understanding `.sf.h5` output content

The output file is an HDF5 container with:

- root-level metadata:
  - `configuration` attribute: the full config used for the run.
- root-level groups:
  - `index` (saved from input HSS index)
  - one group per extracted feature, e.g. `radial`, `lamina`, `icp`, ...

For each feature group (e.g. `/radial`), datasets include:

- `matrix`: per-structure values, shape `(n_bead, n_struct)`
- `mean_arr`: haploid mean across structures
- `std_arr`: haploid std across structures
- `mean_arr_lnorm_gwide`: log2-normalized mean (genome-wide)
- `mean_arr_lnorm_cwide`: log2-normalized mean (chromosome-wide)
- `std_arr_lnorm_gwide`: log2-normalized std (genome-wide)
- `std_arr_lnorm_cwide`: log2-normalized std (chromosome-wide)
- optional `association_freq`: present when `contact_threshold` is provided for that feature

### Python example: inspect and read with `h5py`

```python
import h5py

sf_path = "./model.sf.h5"
with h5py.File(sf_path, "r") as h5:
    print("root keys:", list(h5.keys()))
    print("config:\n", h5.attrs["configuration"])

    feature = "radial"
    print(f"datasets in /{feature}:", list(h5[feature].keys()))

    mat = h5[f"{feature}/matrix"][:]  # (n_bead, n_struct)
    mean_arr = h5[f"{feature}/mean_arr"][:]
    print("matrix shape:", mat.shape)
    print("mean shape:", mean_arr.shape)
```

### Python example: read via `SfFile`

```python
from genome3danalysis.structfeat import SfFile

sf = SfFile("./model.sf.h5", "r")

print("features:", sf.feature_list)
for feat in sf.feature_list:
    print(feat, sf.get_feature_dataset_list(feat))

rad_matrix = sf.get_feature("radial", format="matrix")
lam_mean = sf.get_feature("lamina", format="mean_arr")

sf.close()
```
