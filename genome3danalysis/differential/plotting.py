"""Plotting functions for differential structural-feature analysis."""

from __future__ import annotations

import itertools
import os
from math import ceil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2

from genome3danalysis.core import Sample

DEFAULT_COLORS = {
    "treatment": "seagreen",
    "control": "dimgray",
    "up": "red",
    "down": "blue",
    "cen": "gray",
    "ns": "lightgray",
}


def _resolve_colors(user_colors: dict | None) -> dict[str, str]:
    out = dict(DEFAULT_COLORS)
    if user_colors:
        out.update(user_colors)
    return out


def _save_show(fig: plt.Figure, out_path: str | None) -> None:
    if out_path:
        root, ext = os.path.splitext(out_path)
        base = root if ext.lower() in {".png", ".pdf"} else out_path
        os.makedirs(os.path.dirname(os.path.abspath(base)), exist_ok=True)
        fig.savefig(base + ".png", dpi=200, bbox_inches="tight")
        fig.savefig(base + ".pdf", dpi=200, bbox_inches="tight")
    plt.show()


def _load_cen(cen_file_path: str | None) -> pd.DataFrame:
    if not cen_file_path or not os.path.exists(cen_file_path):
        return pd.DataFrame(columns=["chr", "start", "end", "domain"])
    cen = pd.read_csv(cen_file_path, sep=r"\s+", header=None, names=["chr", "start", "end", "domain"])
    return cen[cen["domain"] == "cen"].copy()


def plot_sf_tracks(
    df: pd.DataFrame,
    treatment: Sample,
    chroms,
    feature_name: str,
    control: Sample | None = None,
    value_col: str | None = None,
    ylim: tuple[float, float] | None = None,
    cen_file_path: str | None = None,
    colors: dict | None = None,
    out_path: str | None = None,
    title: str | None = None,
    h5_dataset: str = "mean_arr",
):
    c = _resolve_colors(colors)
    value_col = feature_name if value_col is None else value_col

    rep_t = treatment.load(feature_name=feature_name, h5_dataset=h5_dataset, value_col_name=value_col)
    rep_c = (
        control.load(feature_name=feature_name, h5_dataset=h5_dataset, value_col_name=value_col)
        if control is not None
        else {}
    )
    rep_all = {**rep_t, **rep_c}
    keys = list(rep_all.keys())

    n_rows = len(chroms)
    fig, axes = plt.subplots(n_rows, 1, figsize=(16, max(3.2 * n_rows, 4)), sharex=False)
    if n_rows == 1:
        axes = [axes]

    linestyles = ["-", "--", ":", "-."]
    cen = _load_cen(cen_file_path)

    for ax, chrom in zip(axes, chroms):
        for i, key in enumerate(keys):
            d = rep_all[key]
            sub = d[d["chr"] == chrom]
            if sub.empty:
                continue
            xv = (sub["start"].values + sub["end"].values) / 2 / 1e6
            cond_color = c["treatment"] if key.startswith(f"{treatment.name}_") else c["control"]
            ax.plot(
                xv,
                sub[value_col].values,
                color=cond_color,
                linestyle=linestyles[i % len(linestyles)],
                linewidth=1.2,
                alpha=0.85,
                label=key,
            )

        cen_chr = cen[cen["chr"] == chrom]
        for _, row in cen_chr.iterrows():
            ax.axvspan(row["start"] / 1e6, row["end"] / 1e6, color=c["cen"], alpha=0.25, zorder=0)

        if control is not None and "significant" in df.columns:
            sig_chr = df[(df["chr"] == chrom) & (df["significant"])].copy()
            for _, srow in sig_chr.iterrows():
                cc = c["up"] if srow["regulation"] == "up" else c["down"]
                ax.axvspan(srow["start"] / 1e6, srow["end"] / 1e6, color=cc, alpha=0.18, zorder=1)

        if ylim is not None:
            ax.set_ylim(*ylim)
        ax.grid(True, alpha=0.3)
        ax.set_title(chrom, loc="left")
        ax.set_ylabel(value_col)
        ax.set_xlabel("Position (Mb)")

    if title:
        fig.suptitle(title)
    fig.tight_layout()
    _save_show(fig, out_path)
    return fig


def plot_diff_tracks(df, chroms, value="delta", cen_file_path=None, colors=None, out_path=None, title=None):
    c = _resolve_colors(colors)
    fig, axes = plt.subplots(len(chroms), 1, figsize=(16, max(3.2 * len(chroms), 4)))
    if len(chroms) == 1:
        axes = [axes]
    cen = _load_cen(cen_file_path)

    for ax, chrom in zip(axes, chroms):
        sub = df[df["chr"] == chrom].copy()
        if sub.empty:
            continue
        x = (sub["start"].values + sub["end"].values) / 2 / 1e6
        y = sub[value].values
        ax.plot(x, y, color="black", linewidth=1.0)
        ax.axhline(0, color="black", linestyle="--", linewidth=0.9)

        sig = sub[sub.get("significant", False)] if "significant" in sub.columns else pd.DataFrame()
        if not sig.empty:
            xu = (sig["start"] + sig["end"]) / 2 / 1e6
            yu = sig[value]
            cols = np.where(sig["regulation"].values == "up", c["up"], c["down"])
            ax.scatter(xu, yu, c=cols, s=10)

        cen_chr = cen[cen["chr"] == chrom]
        for _, row in cen_chr.iterrows():
            ax.axvspan(row["start"] / 1e6, row["end"] / 1e6, color=c["cen"], alpha=0.25, zorder=0)

        ax.set_title(chrom, loc="left")
        ax.set_ylabel(value)
        ax.set_xlabel("Position (Mb)")

    if title:
        fig.suptitle(title)
    fig.tight_layout()
    _save_show(fig, out_path)
    return fig


def plot_manhattan(df, gate_threshold=None, colors=None, out_path=None, title=None):
    c = _resolve_colors(colors)
    data = df.copy()
    chroms = list(dict.fromkeys(data["chr"].tolist()))
    x_all = np.zeros(len(data))
    offset = 0
    boundaries = []

    for ch in chroms:
        idx = data["chr"] == ch
        starts = data.loc[idx, "start"].values
        x_all[idx.values] = starts + offset
        if starts.size:
            offset = int(x_all[idx.values].max()) + int(data.loc[idx, "end"].max() - data.loc[idx, "start"].min())
        boundaries.append(offset)

    y = -np.log10(np.clip(data["pvalue"].values.astype(float), 1e-300, 1.0))
    cols = np.where(data.get("regulation", "ns") == "up", c["up"], np.where(data.get("regulation", "ns") == "down", c["down"], c["ns"]))

    fig, ax = plt.subplots(figsize=(18, 5))
    ax.scatter(x_all, y, c=cols, s=8, alpha=0.8)
    for b in boundaries:
        ax.axvline(b, color="lightgray", linewidth=0.5)
    if gate_threshold is not None:
        ax.axhline(-np.log10(max(gate_threshold, 1e-300)), color="black", linestyle="--")
    ax.set_ylabel("-log10(pvalue)")
    ax.set_xlabel("Genome position")
    if title:
        ax.set_title(title)
    fig.tight_layout()
    _save_show(fig, out_path)
    return fig


def plot_volcano(
    df,
    gate_threshold=None,
    delta_threshold=None,
    colors=None,
    out_path=None,
    title=None,
    x_mode="delta",
):
    """Plot volcano scatter for differential regions.

    Parameters
    ----------
    df : pandas.DataFrame
        Differential result table.
    gate_threshold : float | None, optional
        P-value (or FDR) threshold to draw as horizontal dashed line.
    delta_threshold : float | None, optional
        Absolute x-threshold to draw as vertical dashed lines.
    colors : dict | None, optional
        Color overrides.
    out_path : str | None, optional
        Base output path; when set saves PNG/PDF and shows.
    title : str | None, optional
        Plot title.
    x_mode : str, optional
        X-axis mode:
        - ``"delta"``: use ``df["delta"]`` directly.
        - ``"log2fc_mean"``: use mean of available ``log2fc_*`` columns,
          e.g. ``(log2fc_11 + log2fc_22)/2``.
    """
    c = _resolve_colors(colors)

    if x_mode == "delta":
        x = df["delta"].values.astype(float)
        x_label = "delta"
    elif x_mode == "log2fc_mean":
        fc_cols = sorted([col for col in df.columns if col.startswith("log2fc_")])
        if not fc_cols:
            raise ValueError(
                "x_mode='log2fc_mean' requires at least one 'log2fc_*' column in df"
            )
        x = df[fc_cols].mean(axis=1).values.astype(float)
        x_label = "mean log2 fold change"
    else:
        raise ValueError("x_mode must be one of {'delta', 'log2fc_mean'}")

    y = -np.log10(np.clip(df["pvalue"].values.astype(float), 1e-300, 1.0))
    reg = df.get("regulation", pd.Series(["ns"] * len(df))).values
    cols = np.where(reg == "up", c["up"], np.where(reg == "down", c["down"], c["ns"]))

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(x, y, c=cols, s=10, alpha=0.75)
    if delta_threshold is not None:
        ax.axvline(delta_threshold, color="black", linestyle="--")
        ax.axvline(-delta_threshold, color="black", linestyle="--")
    if gate_threshold is not None:
        ax.axhline(-np.log10(max(gate_threshold, 1e-300)), color="black", linestyle="--")
    ax.set_xlabel(x_label)
    ax.set_ylabel("-log10(pvalue)")
    if title:
        ax.set_title(title)
    fig.tight_layout()
    _save_show(fig, out_path)
    return fig


def plot_qq(df, md_mode="contrast", colors=None, out_path=None, title=None):
    _ = _resolve_colors(colors)
    obs = np.sort(df["md2_used"].values.astype(float))
    n = len(obs)
    probs = (np.arange(1, n + 1) - 0.5) / n

    if md_mode == "contrast":
        df_chi = 1
    else:
        n_samples = len([c for c in df.columns if c.startswith("val_")])
        df_chi = max(n_samples, 1)

    theo = chi2.ppf(probs, df=df_chi)
    lam = float(np.median(obs) / max(np.median(theo), 1e-30))

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(theo, obs, s=10, alpha=0.8)
    lim = max(float(np.nanmax(theo)), float(np.nanmax(obs)))
    ax.plot([0, lim], [0, lim], "k--", linewidth=1)
    ax.text(0.03, 0.95, f"lambda={lam:.3f}", transform=ax.transAxes, va="top")
    ax.set_xlabel("Theoretical chi-square quantile")
    ax.set_ylabel("Observed md2_used quantile")
    if title:
        ax.set_title(title)
    fig.tight_layout()
    _save_show(fig, out_path)
    return fig


def plot_replicate_scatter(rep_data: dict, sample_name: str, colors=None, out_path=None):
    c = _resolve_colors(colors)
    keys = list(rep_data.keys())
    pairs = list(itertools.combinations(keys, 2))
    n = len(pairs)
    if n == 0:
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.set_title(f"{sample_name}: only one replicate")
        _save_show(fig, out_path)
        return fig

    ncols = min(3, n)
    nrows = ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = np.atleast_1d(axes).ravel()

    for ax, (k1, k2) in zip(axes, pairs):
        d1, d2 = rep_data[k1], rep_data[k2]
        m = d1.merge(d2, on=["chr", "start", "end"], suffixes=("_1", "_2"))
        v1 = m.iloc[:, 3].values
        v2 = m.iloc[:, 4].values
        r = np.corrcoef(v1, v2)[0, 1] if len(v1) > 1 else np.nan
        ax.scatter(v1, v2, s=8, alpha=0.5, color=c["treatment"])
        lo = np.nanmin([v1.min(initial=0), v2.min(initial=0)])
        hi = np.nanmax([v1.max(initial=1), v2.max(initial=1)])
        ax.plot([lo, hi], [lo, hi], "k--", linewidth=1)
        ax.set_title(f"{k1} vs {k2} (r={r:.3f})")
        ax.set_xlabel(k1)
        ax.set_ylabel(k2)

    for ax in axes[n:]:
        ax.axis("off")

    fig.suptitle(f"Replicate scatter: {sample_name}")
    fig.tight_layout()
    _save_show(fig, out_path)
    return fig


def plot_md2_vs_delta(df, colors=None, out_path=None, title=None):
    c = _resolve_colors(colors)
    x = np.abs(df["delta"].values)
    y = df["md2_used"].values
    reg = df.get("regulation", pd.Series(["ns"] * len(df))).values
    cols = np.where(reg == "up", c["up"], np.where(reg == "down", c["down"], c["ns"]))

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(x, y, c=cols, s=10, alpha=0.7)
    if "delta_cut" in df.columns and len(df):
        ax.axvline(float(df["delta_cut"].iloc[0]), color="black", linestyle="--")
    ax.set_xlabel("|delta|")
    ax.set_ylabel("md2_used")
    if title:
        ax.set_title(title)
    fig.tight_layout()
    _save_show(fig, out_path)
    return fig


def plot_sample_correlation(rep_data: dict, colors=None, out_path=None, title=None):
    _ = _resolve_colors(colors)
    keys = list(rep_data.keys())
    if not keys:
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.set_title("No replicate data")
        _save_show(fig, out_path)
        return fig

    merged = None
    for k in keys:
        d = rep_data[k].copy().rename(columns={rep_data[k].columns[3]: k})
        merged = d if merged is None else merged.merge(d, on=["chr", "start", "end"], how="inner")
    mat = merged[keys].corr(method="pearson")

    off_diag = mat.values[~np.eye(mat.shape[0], dtype=bool)]
    vmin = float(np.nanmin(off_diag)) if off_diag.size else -1.0

    fig, ax = plt.subplots(figsize=(1.1 * len(keys) + 3, 1.0 * len(keys) + 2))
    im = ax.imshow(mat.values, cmap="coolwarm", vmin=vmin, vmax=1.0)
    ax.set_xticks(np.arange(len(keys)))
    ax.set_yticks(np.arange(len(keys)))
    ax.set_xticklabels(keys, rotation=45, ha="right")
    ax.set_yticklabels(keys)
    for i in range(len(keys)):
        for j in range(len(keys)):
            ax.text(j, i, f"{mat.values[i, j]:.2f}", ha="center", va="center", fontsize=9)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    if title:
        ax.set_title(title)
    fig.tight_layout()
    _save_show(fig, out_path)
    return fig


def plot_delta_histogram(df, bins=80, colors=None, out_path=None, title=None):
    c = _resolve_colors(colors)
    d = df["delta"].values
    up = df[df.get("regulation", "ns") == "up"]["delta"].values
    down = df[df.get("regulation", "ns") == "down"]["delta"].values

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(d, bins=bins, color=c["ns"], alpha=0.8, label="all")
    if len(up):
        ax.hist(up, bins=bins, color=c["up"], alpha=0.6, label="up")
    if len(down):
        ax.hist(down, bins=bins, color=c["down"], alpha=0.6, label="down")

    if "delta_cut" in df.columns and len(df):
        dc = float(df["delta_cut"].iloc[0])
        ax.axvline(dc, color="black", linestyle="--")
        ax.axvline(-dc, color="black", linestyle="--")

    ax.set_xlabel("delta")
    ax.set_ylabel("count")
    ax.legend()
    if title:
        ax.set_title(title)
    fig.tight_layout()
    _save_show(fig, out_path)
    return fig
