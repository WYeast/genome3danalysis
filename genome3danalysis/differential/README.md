# genome3danalysis.differential

Mahalanobis-distance-based differential structural feature (SF) analysis for sf.h5 inputs.

## Example A: End-to-end (LPS vs control, 2 reps each)

```python
from genome3danalysis import Sample
from genome3danalysis.differential import (
    compute_differential_sf,
    plot_manhattan,
    plot_volcano,
    plot_sf_tracks,
)

treatment = Sample(
    name="LPS",
    h5_paths=[
        "/path/to/LPS_1_model.sf.h5",
        "/path/to/LPS_2_model.sf.h5",
    ],
    mode="bulk",
)

control = Sample(
    name="unstimulated",
    h5_paths=[
        "/path/to/unstimulated_1_model.sf.h5",
        "/path/to/unstimulated_2_model.sf.h5",
    ],
    mode="bulk",
)

res = compute_differential_sf(
    treatment=treatment,
    control=control,
    feature_name="icp",
    h5_dataset="mean_arr",
    preprocess="raw",
    md_mode="contrast",
    pvalue_threshold=0.05,
    fdr_threshold=None,
    delta_mode="zero",
    cen_file_path="/path/to/domains_100kb.bed",
    output_csv="./results/LPS_vs_unstim_icp.csv",
    output_bed_dir="./results/bed",
    verbose=True,
)

plot_manhattan(res, gate_threshold=0.05, out_path="./results/manhattan")
plot_volcano(
    res,
    gate_threshold=0.05,
    delta_threshold=float(res["delta_cut"].iloc[0]),
    x_mode="delta",  # or "log2fc_mean"
    gene_list=["NFKB1"],
    gene_annotation_file="/u/project/xjzhou/wyeast/wy_data/Sohyeon_data/deep_seq/mcool_files/analysis/hg38.genes.chr.bed",
    show=True,
    out_path="./results/volcano",
)
plot_sf_tracks(
    res,
    treatment=treatment,
    control=control,
    chroms=["chr4", "chr9"],
    feature_name="icp",
    h5_dataset="mean_arr",
    out_path="./results/tracks",
)
```

## Example B: QC only (single sample or list of samples)

```python
from genome3danalysis import Sample
from genome3danalysis.differential import run_qc

samples = [
    Sample("LPS", ["/path/LPS_1.sf.h5", "/path/LPS_2.sf.h5"], mode="bulk"),
    Sample("unstimulated", ["/path/unstim_1.sf.h5", "/path/unstim_2.sf.h5"], mode="bulk"),
]

qc = run_qc(
    samples,
    feature_name="enhd",
    h5_dataset="mean_arr",
    cen_file_path="/path/to/domains_100kb.bed",
    output_dir="./qc",
    verbose=True,
)
```

## Example C: Low-level custom workflow

```python
import numpy as np
from statsmodels.stats.multitest import multipletests

from genome3danalysis import Sample
from genome3danalysis.differential import (
    build_matrix,
    preprocess_matrix,
    mahalanobis_robust,
    compute_pvalues,
    gate_significance,
)

t = Sample("LPS", ["/path/LPS_1.sf.h5", "/path/LPS_2.sf.h5"], mode="bulk")
c = Sample("unstimulated", ["/path/u_1.sf.h5", "/path/u_2.sf.h5"], mode="bulk")

merged, sample_cols = build_matrix(t, c, feature_name="icp", h5_dataset="mean_arr")
X_raw = merged[sample_cols].to_numpy(float)
X = preprocess_matrix(X_raw, mode="log2p1")
_, mu, cov, used = mahalanobis_robust(X, outlier_frac=0.05)

stats = compute_pvalues(
    X,
    mu,
    cov,
    sample_cols,
    treatment_name=t.name,
    control_name=c.name,
    n_treat_reps=2,
    n_ctrl_reps=2,
    md_mode="contrast",
)

_, fdr, _, _ = multipletests(stats["pvals"], method="fdr_bh")
out = merged[["chr", "start", "end"]].copy()
out["pvalue"] = stats["pvals"]
out["fdr"] = fdr
out["delta"] = X[:, [sample_cols.index("LPS_1"), sample_cols.index("LPS_2")]].mean(1) - \
               X[:, [sample_cols.index("unstimulated_1"), sample_cols.index("unstimulated_2")]].mean(1)
out["used_for_cov"] = used
out = gate_significance(out, fdr_threshold=0.05, pvalue_threshold=None, delta_mode="zero")
```
