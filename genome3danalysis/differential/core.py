"""High-level differential structural-feature analysis entry point."""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from genome3danalysis.core import Sample
from genome3danalysis.differential.io import build_matrix, filter_cen_overlap
from genome3danalysis.differential.stats import (
    compute_pairwise_fc,
    compute_pvalues,
    gate_significance,
    mahalanobis_robust,
    preprocess_matrix,
)


def compute_differential_sf(
    treatment: Sample,
    control: Sample,
    feature_name: str,
    *,
    h5_dataset: str = "mean_arr",
    preprocess: str = "log2p1",
    md_mode: str = "contrast",
    robust_outlier_frac: float = 0.05,
    fdr_threshold: float | None = 0.05,
    pvalue_threshold: float | None = None,
    delta_mode: str = "zero",
    delta_sd_multiplier: float = 0.5,
    exclude_chroms: tuple[str, ...] = ("chrX", "chrY", "X", "Y", "chrM", "M", "chrMT", "MT"),
    cen_file_path: str | None = None,
    keep_pairwise_fc: bool = True,
    output_csv: str | None = None,
    output_bed_dir: str | None = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """Run MD-based differential analysis between treatment and control samples.

    Notes
    -----
    ``md2`` in output is the statistic used for p-value computation:
    - full mode: full Mahalanobis distance squared
    - contrast mode: projected z^2

    ``md2_full`` is always reported as the full Mahalanobis distance squared.
    """
    if (fdr_threshold is None) == (pvalue_threshold is None):
        raise ValueError("Exactly one of fdr_threshold or pvalue_threshold must be non-None")

    gate_desc = (
        f"pvalue<{pvalue_threshold}" if pvalue_threshold is not None else f"fdr<{fdr_threshold}"
    )

    if verbose:
        print(f"\n{'='*70}\nMD differential: {treatment.name} vs {control.name}  ({feature_name})")
        print(f"  input_format=h5, h5_dataset={h5_dataset}")
        print(f"  preprocess={preprocess}  md_mode={md_mode}  {gate_desc}  delta_mode={delta_mode}")
        if cen_file_path:
            print(f"  CEN filter: {cen_file_path}")
        print("=" * 70)

    merged, sample_cols = build_matrix(
        treatment=treatment,
        control=control,
        feature_name=feature_name,
        h5_dataset=h5_dataset,
    )
    n_total_in = len(merged)
    raw_df_for_fc = merged.copy() if keep_pairwise_fc else None

    if exclude_chroms:
        keep = ~merged["chr"].isin(set(exclude_chroms))
        merged = merged[keep].reset_index(drop=True)
        if raw_df_for_fc is not None:
            raw_df_for_fc = raw_df_for_fc[keep].reset_index(drop=True)
    n_after_chrom = len(merged)

    n_drop_cen = 0
    if cen_file_path and os.path.exists(cen_file_path):
        merged_new, n_drop_cen = filter_cen_overlap(merged, cen_file_path)
        if n_drop_cen > 0 and raw_df_for_fc is not None:
            keep_keys = set(zip(merged_new["chr"], merged_new["start"]))
            mask_keep = [
                (c, s) in keep_keys
                for c, s in zip(raw_df_for_fc["chr"], raw_df_for_fc["start"])
            ]
            raw_df_for_fc = raw_df_for_fc[mask_keep].reset_index(drop=True)
        merged = merged_new

    if verbose:
        print(
            f"  Bins loaded: {n_total_in}, after chrom filter: {n_after_chrom}, "
            f"after CEN filter: {len(merged)} (dropped {n_drop_cen} CEN bins)"
        )
        print(f"  Sample columns (in order): {sample_cols}")

    sample_mat_raw = merged[sample_cols].values.astype(float)
    finite_mask = np.all(np.isfinite(sample_mat_raw), axis=1)
    if not finite_mask.all():
        if verbose:
            print(f"  Dropping {(~finite_mask).sum()} bins with NaN/Inf")
        merged = merged[finite_mask].reset_index(drop=True)
        if raw_df_for_fc is not None:
            raw_df_for_fc = raw_df_for_fc[finite_mask].reset_index(drop=True)
        sample_mat_raw = merged[sample_cols].values.astype(float)

    sample_mat = preprocess_matrix(sample_mat_raw, mode=preprocess)

    _, mean_used, cov_used, used_mask = mahalanobis_robust(sample_mat, robust_outlier_frac)
    n_treat = len([c for c in sample_cols if c.startswith(f"{treatment.name}_")])
    n_ctrl = len([c for c in sample_cols if c.startswith(f"{control.name}_")])

    pstats = compute_pvalues(
        sample_mat,
        mean_used,
        cov_used,
        sample_cols,
        treatment_name=treatment.name,
        control_name=control.name,
        n_treat_reps=n_treat,
        n_ctrl_reps=n_ctrl,
        md_mode=md_mode,
    )
    md2_full = pstats["md2_full"]
    md2_used = pstats["md2_used"]
    pvals = pstats["pvals"]

    _, fdr_vals, _, _ = multipletests(pvals, method="fdr_bh")

    treat_idx = [i for i, c in enumerate(sample_cols) if c.startswith(f"{treatment.name}_")]
    ctrl_idx = [i for i, c in enumerate(sample_cols) if c.startswith(f"{control.name}_")]
    delta = sample_mat[:, treat_idx].mean(axis=1) - sample_mat[:, ctrl_idx].mean(axis=1)

    out = merged[["chr", "start", "end"]].copy()
    for j, name in enumerate(sample_cols):
        out[f"val_{name}"] = sample_mat[:, j]
    out["md2_full"] = md2_full
    out["md2"] = md2_used
    out["md2_used"] = md2_used
    out["md_mode"] = md_mode
    out["pvalue"] = pvals
    out["fdr"] = fdr_vals
    out["delta"] = delta
    out["used_for_cov"] = used_mask

    out = gate_significance(
        out,
        fdr_threshold=fdr_threshold,
        pvalue_threshold=pvalue_threshold,
        delta_mode=delta_mode,
        delta_sd_multiplier=delta_sd_multiplier,
    )

    if keep_pairwise_fc and raw_df_for_fc is not None:
        fc = compute_pairwise_fc(raw_df_for_fc, treatment.name, control.name)
        for key, val in fc.items():
            out[key] = val
            if key.startswith("fc_"):
                suffix = key.split("_", 1)[1]
                out[f"fdr_{suffix}"] = fdr_vals

    if output_csv is not None:
        os.makedirs(os.path.dirname(os.path.abspath(output_csv)), exist_ok=True)
        out.to_csv(output_csv, index=False)

    if verbose:
        n_sig = int(out["significant"].sum())
        n_up = int((out["regulation"] == "up").sum())
        n_dn = int((out["regulation"] == "down").sum())
        print(f"  Result: {n_sig} sig bins  ({n_up} up, {n_dn} down)")
        print(f"  Sig regions: all={n_sig}, up={n_up}, down={n_dn}")
        print(
            "  delta_cut = "
            f"{float(out['delta_cut'].iloc[0]):.6g}  "
            f"(SD(delta)={np.std(out['delta'].values):.6g}, mean(delta)={np.mean(out['delta'].values):.6g})"
        )
        if output_csv is not None:
            print(f"  Saved CSV: {output_csv}")

    if output_bed_dir is not None:
        os.makedirs(output_bed_dir, exist_ok=True)
        up_p = os.path.join(output_bed_dir, f"{treatment.name}_vs_{control.name}_{feature_name}_up.bed")
        dn_p = os.path.join(output_bed_dir, f"{treatment.name}_vs_{control.name}_{feature_name}_down.bed")
        out[out["regulation"] == "up"][["chr", "start", "end"]].to_csv(
            up_p, sep="\t", header=False, index=False
        )
        out[out["regulation"] == "down"][["chr", "start", "end"]].to_csv(
            dn_p, sep="\t", header=False, index=False
        )
        if verbose:
            print(f"  Saved BED: {up_p}\n  Saved BED: {dn_p}")

    return out
