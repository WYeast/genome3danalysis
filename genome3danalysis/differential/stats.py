"""Statistical methods for differential structural-feature analysis."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import chi2, rankdata


def preprocess_matrix(X: np.ndarray, mode: str) -> np.ndarray:
    """Per-column preprocessing of a ``(n_bins, n_samples)`` matrix."""
    X = np.asarray(X, dtype=float).copy()
    n_bins, n_samples = X.shape

    if mode == "raw":
        return X
    if mode == "minmax":
        out = np.empty_like(X)
        for j in range(n_samples):
            col = X[:, j]
            cmin, cmax = np.nanmin(col), np.nanmax(col)
            out[:, j] = 0.0 if cmax - cmin == 0 else (col - cmin) / (cmax - cmin)
        return out
    if mode == "log2p1":
        return np.log2(X + 1.0)
    if mode == "log2":
        positives = X[X > 0]
        if positives.size == 0:
            raise ValueError("log2 mode but matrix has no positive values")
        eps = positives.min() / 2.0
        return np.log2(X + eps)
    if mode == "quantile":
        ranks = np.empty_like(X)
        sorted_vals = np.empty_like(X)
        for j in range(n_samples):
            ranks[:, j] = rankdata(X[:, j], method="average")
            sorted_vals[:, j] = np.sort(X[:, j])
        target = sorted_vals.mean(axis=1)
        out = np.empty_like(X)
        for j in range(n_samples):
            r = ranks[:, j]
            r_lo = np.clip(np.floor(r).astype(int) - 1, 0, n_bins - 1)
            r_hi = np.clip(np.ceil(r).astype(int) - 1, 0, n_bins - 1)
            frac = r - np.floor(r)
            out[:, j] = target[r_lo] * (1 - frac) + target[r_hi] * frac
        return out

    raise ValueError(f"Unknown preprocess mode: {mode!r}")


def mahalanobis_robust(
    X: np.ndarray,
    outlier_frac: float = 0.05,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Two-pass robust Mahalanobis distance."""
    n_bins, n_samples = X.shape
    mean_all = X.mean(axis=0)
    cov_all = np.cov(X, rowvar=False)
    cov_inv = np.linalg.pinv(cov_all)
    diff = X - mean_all
    md2_init = np.einsum("ij,jk,ik->i", diff, cov_inv, diff)

    if outlier_frac <= 0:
        return md2_init, mean_all, cov_all, np.ones(n_bins, dtype=bool)

    cutoff = np.quantile(md2_init, 1.0 - outlier_frac)
    keep_mask = md2_init <= cutoff
    if keep_mask.sum() < n_samples + 5:
        return md2_init, mean_all, cov_all, np.ones(n_bins, dtype=bool)

    Xk = X[keep_mask]
    mean_final = Xk.mean(axis=0)
    cov_final = np.cov(Xk, rowvar=False)
    cov_inv_final = np.linalg.pinv(cov_final)
    diff_final = X - mean_final
    md2_final = np.einsum("ij,jk,ik->i", diff_final, cov_inv_final, diff_final)
    return md2_final, mean_final, cov_final, keep_mask


def compute_pvalues(
    X: np.ndarray,
    mean: np.ndarray,
    cov: np.ndarray,
    sample_cols: list[str],
    treatment_name: str,
    control_name: str,
    n_treat_reps: int,
    n_ctrl_reps: int,
    md_mode: str = "contrast",
) -> dict[str, np.ndarray]:
    """Compute full/contrast MD statistics and p-values."""
    diff = X - mean
    cov_inv = np.linalg.pinv(cov)
    md2_full = np.einsum("ij,jk,ik->i", diff, cov_inv, diff)
    n_samples = X.shape[1]

    if md_mode == "full":
        md2_used = md2_full
        pvals = chi2.sf(md2_used, df=n_samples)
    elif md_mode == "contrast":
        treat_keys = [f"{treatment_name}_{i+1}" for i in range(n_treat_reps)]
        ctrl_keys = [f"{control_name}_{i+1}" for i in range(n_ctrl_reps)]
        treat_idx = [sample_cols.index(k) for k in treat_keys if k in sample_cols]
        ctrl_idx = [sample_cols.index(k) for k in ctrl_keys if k in sample_cols]
        if len(treat_idx) == 0 or len(ctrl_idx) == 0:
            raise ValueError("Could not find treatment/control columns for contrast mode")

        contrast_vec = np.zeros(n_samples)
        contrast_vec[treat_idx] = 1.0 / len(treat_idx)
        contrast_vec[ctrl_idx] = -1.0 / len(ctrl_idx)
        var_proj = float(contrast_vec @ cov @ contrast_vec)
        proj = diff @ contrast_vec
        md2_used = (proj**2) / max(var_proj, 1e-30)
        pvals = chi2.sf(md2_used, df=1)
    else:
        raise ValueError(f"Unknown md_mode: {md_mode}")

    pvals = np.clip(pvals, 1e-300, 1.0)
    return {"md2_full": md2_full, "md2_used": md2_used, "pvals": pvals}


def compute_pairwise_fc(
    raw_df: pd.DataFrame,
    treatment_name: str,
    control_name: str,
    eps: float = 1.0,
) -> dict[str, np.ndarray]:
    """Compute index-matched pairwise fold-changes for available replicates."""
    t_cols = [c for c in raw_df.columns if c.startswith(f"{treatment_name}_")]
    c_cols = [c for c in raw_df.columns if c.startswith(f"{control_name}_")]

    def _rep_idx(col: str) -> int:
        try:
            return int(col.rsplit("_", 1)[1])
        except ValueError:
            return 10**9

    t_cols = sorted(t_cols, key=_rep_idx)
    c_cols = sorted(c_cols, key=_rep_idx)
    k = min(len(t_cols), len(c_cols))
    if k == 0:
        return {}

    out: dict[str, np.ndarray] = {}
    for i in range(k):
        pair_id = i + 1
        t = raw_df[t_cols[i]].values + eps
        c = raw_df[c_cols[i]].values + eps
        fc = t / c
        out[f"fc_{pair_id}{pair_id}"] = fc
        out[f"log2fc_{pair_id}{pair_id}"] = np.log2(fc)
    return out


def gate_significance(
    df: pd.DataFrame,
    *,
    fdr_threshold: float | None = None,
    pvalue_threshold: float | None = None,
    delta_mode: str = "zero",
    delta_sd_multiplier: float = 0.5,
) -> pd.DataFrame:
    """Apply significance gating and label regulation direction."""
    if (fdr_threshold is None) == (pvalue_threshold is None):
        raise ValueError("Exactly one of fdr_threshold or pvalue_threshold must be non-None")

    out = df.copy()
    delta = out["delta"].values

    if delta_mode == "zero":
        delta_cut = 0.0
    elif delta_mode == "sd":
        delta_cut = float(delta_sd_multiplier * np.std(delta))
    else:
        raise ValueError(f"Unknown delta_mode: {delta_mode}")

    if pvalue_threshold is not None:
        stat_pass = out["pvalue"].values < pvalue_threshold
    else:
        stat_pass = out["fdr"].values < float(fdr_threshold)

    delta_pass = np.abs(delta) > delta_cut
    significant = stat_pass & delta_pass
    regulation = np.where(
        significant & (delta > 0),
        "up",
        np.where(significant & (delta < 0), "down", "ns"),
    )

    out["delta_cut"] = delta_cut
    out["significant"] = significant
    out["regulation"] = regulation
    return out
