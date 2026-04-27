"""Quality-control utilities for differential structural-feature workflows."""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

from genome3danalysis.core import Sample
from genome3danalysis.differential.io import filter_cen_overlap
from genome3danalysis.differential.plotting import plot_replicate_scatter, plot_sample_correlation


def _rep_stats(df: pd.DataFrame) -> dict[str, float]:
    val = df.iloc[:, 3].values.astype(float)
    return {
        "mean": float(np.nanmean(val)),
        "median": float(np.nanmedian(val)),
        "iqr": float(np.nanpercentile(val, 75) - np.nanpercentile(val, 25)),
        "nan_count": int(np.isnan(val).sum()),
        "n_bins": int(len(val)),
    }


def _plot_hist(rep_data: dict[str, pd.DataFrame], out_path: str | None = None) -> None:
    import matplotlib.pyplot as plt

    n = len(rep_data)
    fig, axes = plt.subplots(n, 1, figsize=(8, max(2.5 * n, 3)))
    if n == 1:
        axes = [axes]
    for ax, (name, df) in zip(axes, rep_data.items()):
        vals = df.iloc[:, 3].values.astype(float)
        ax.hist(vals[np.isfinite(vals)], bins=80, color="gray", alpha=0.8)
        ax.set_title(name)
        ax.set_xlabel("value")
        ax.set_ylabel("count")
    fig.tight_layout()
    if out_path:
        os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
        base = out_path.rsplit(".", 1)[0]
        fig.savefig(base + ".png", dpi=200, bbox_inches="tight")
        fig.savefig(base + ".pdf", dpi=200, bbox_inches="tight")
    plt.show()


def run_qc(
    samples: Sample | list[Sample],
    feature_name: str,
    *,
    h5_dataset: str = "mean_arr",
    cen_file_path: str | None = None,
    exclude_chroms: tuple[str, ...] = ("chrX", "chrY", "X", "Y", "chrM", "M", "chrMT", "MT"),
    output_dir: str | None = None,
    colors: dict | None = None,
    verbose: bool = True,
) -> dict:
    """Run QC for one or multiple samples and generate QC plots/metrics."""
    sample_list = [samples] if isinstance(samples, Sample) else list(samples)
    metrics: dict = {"samples": {}, "inter_sample_corr": None}

    per_sample_mean: dict[str, np.ndarray] = {}

    for sample in sample_list:
        rep_data = sample.load(feature_name=feature_name, h5_dataset=h5_dataset, value_col_name=feature_name)

        filtered: dict[str, pd.DataFrame] = {}
        rep_stats: dict[str, dict] = {}
        rep_names = list(rep_data.keys())
        corr_within = pd.DataFrame(index=rep_names, columns=rep_names, dtype=float)

        for rep_name, df in rep_data.items():
            work = df.copy()
            if exclude_chroms:
                work = work[~work["chr"].isin(set(exclude_chroms))].reset_index(drop=True)
            work, _ = filter_cen_overlap(work, cen_file_path)
            filtered[rep_name] = work
            rep_stats[rep_name] = _rep_stats(work)

        for i, a in enumerate(rep_names):
            for j, b in enumerate(rep_names):
                m = filtered[a].merge(filtered[b], on=["chr", "start", "end"], suffixes=("_a", "_b"))
                va = m.iloc[:, 3].values
                vb = m.iloc[:, 4].values
                corr_within.iloc[i, j] = np.corrcoef(va, vb)[0, 1] if len(m) > 1 else np.nan

        if verbose:
            print(f"\nQC sample: {sample.name}")
            for rep_name in rep_names:
                st = rep_stats[rep_name]
                print(
                    f"  {rep_name}: mean={st['mean']:.6g}, median={st['median']:.6g}, "
                    f"IQR={st['iqr']:.6g}, NaN={st['nan_count']}"
                )
            print("  replicate Pearson r:")
            print(corr_within.to_string())

        out_prefix = None
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            out_prefix = os.path.join(output_dir, sample.name)

        plot_sample_correlation(
            filtered,
            colors=colors,
            out_path=f"{out_prefix}_correlation" if out_prefix else None,
            title=f"{sample.name} replicate correlation",
        )
        plot_replicate_scatter(
            filtered,
            sample_name=sample.name,
            colors=colors,
            out_path=f"{out_prefix}_scatter" if out_prefix else None,
        )
        _plot_hist(filtered, out_path=f"{out_prefix}_hist" if out_prefix else None)

        merged = None
        for k, d in filtered.items():
            renamed = d.rename(columns={feature_name: k})
            merged = renamed if merged is None else merged.merge(renamed, on=["chr", "start", "end"], how="inner")
        per_sample_mean[sample.name] = merged[rep_names].mean(axis=1).values

        metrics["samples"][sample.name] = {
            "rep_stats": rep_stats,
            "replicate_corr": corr_within,
            "n_replicates": len(rep_names),
        }

    if len(sample_list) > 1:
        names = [s.name for s in sample_list]
        inter = pd.DataFrame(index=names, columns=names, dtype=float)
        for i, a in enumerate(names):
            for j, b in enumerate(names):
                va, vb = per_sample_mean[a], per_sample_mean[b]
                n = min(len(va), len(vb))
                inter.iloc[i, j] = np.corrcoef(va[:n], vb[:n])[0, 1] if n > 1 else np.nan
        metrics["inter_sample_corr"] = inter

        if verbose:
            print("\nInter-sample Pearson r (sample means):")
            print(inter.to_string())

    return metrics
