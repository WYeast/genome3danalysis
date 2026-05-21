"""IO utilities for differential structural-feature analysis."""

from __future__ import annotations

import os

import h5py
import numpy as np
import pandas as pd

from genome3danalysis.core import Sample


def load_sample_h5(
    file_path: str,
    feature_name: str,
    h5_dataset: str = "mean_arr",
    value_col_name: str = "value",
) -> pd.DataFrame:
    """Read a sample's per-bin haploid feature value from sf.h5.

    Parameters
    ----------
    file_path
        Path to sf.h5.
    feature_name
        Group inside h5, e.g. ``icp``, ``enhd``, ``rg``, ``lamina``.
    h5_dataset
        Dataset under the feature group, e.g. ``mean_arr``.
    value_col_name
        Output value-column name.

    Returns
    -------
    pandas.DataFrame
        Columns ``['chr', 'start', 'end', value_col_name]``.
    """
    with h5py.File(file_path, "r") as h5:
        if feature_name not in h5:
            raise KeyError(
                f"feature {feature_name!r} not in {file_path}; "
                f"available: {list(h5.keys())}"
            )
        feat_group = h5[feature_name]
        if h5_dataset not in feat_group:
            raise KeyError(
                f"dataset {h5_dataset!r} not in /{feature_name}; "
                f"available: {list(feat_group.keys())}"
            )
        values = feat_group[h5_dataset][:]

        if values.ndim != 1:
            raise ValueError(
                f"Expected 1D /{feature_name}/{h5_dataset} for bulk loading; got shape={values.shape}"
            )

        if "index" not in h5:
            raise KeyError(f"/index group missing in {file_path}")
        idx = h5["index"]
        for k in ["copy", "chromstr", "start", "end"]:
            if k not in idx:
                raise KeyError(f"/index/{k} missing in {file_path}")

        copy_arr = idx["copy"][:]
        first_copy = copy_arr.min()
        haploid_mask = copy_arr == first_copy

        chromstr = idx["chromstr"][:]
        if chromstr.dtype.kind == "S":
            chromstr = np.array([s.decode("utf-8") for s in chromstr])
        starts = idx["start"][:].astype(int)
        ends = idx["end"][:].astype(int)

        chr_h = chromstr[haploid_mask]
        start_h = starts[haploid_mask]
        end_h = ends[haploid_mask]

        if len(chr_h) != len(values):
            raise ValueError(
                f"Length mismatch in {file_path}: "
                f"haploid coords (copy=={first_copy}) = {len(chr_h)}, "
                f"but /{feature_name}/{h5_dataset} has length = {len(values)}"
            )

    return pd.DataFrame(
        {
            "chr": chr_h,
            "start": start_h,
            "end": end_h,
            value_col_name: values,
        }
    )


def build_matrix(
    treatment: Sample,
    control: Sample,
    feature_name: str,
    h5_dataset: str = "mean_arr",
) -> tuple[pd.DataFrame, list[str]]:
    """Load treatment/control sample data and merge on genomic bins."""
    merged: pd.DataFrame | None = None
    sample_cols: list[str] = []

    all_data: dict[str, pd.DataFrame] = {}
    for sample in [control, treatment]:
        rep_data = sample.load(feature_name=feature_name, h5_dataset=h5_dataset, value_col_name=feature_name)
        for key, df in rep_data.items():
            sample_cols.append(key)
            all_data[key] = df.rename(columns={feature_name: key})

    for key in sample_cols:
        df = all_data[key]
        merged = df if merged is None else merged.merge(df, on=["chr", "start", "end"], how="inner")

    if merged is None:
        raise ValueError("No sample data could be loaded")
    return merged, sample_cols


def filter_cen_overlap(df: pd.DataFrame, cen_file_path: str | None) -> tuple[pd.DataFrame, int]:
    """Return ``(filtered_df, n_dropped)`` by removing CEN-overlapping bins."""
    if not cen_file_path or not os.path.exists(cen_file_path):
        return df, 0

    cen = pd.read_csv(cen_file_path, sep=r"\s+", header=None, names=["chr", "start", "end", "domain"])
    cen = cen[cen["domain"] == "cen"].copy()
    if cen.empty:
        return df, 0

    cen_by_chr = {c: sub for c, sub in cen.groupby("chr")}

    def _hits(row: pd.Series) -> bool:
        sub = cen_by_chr.get(row["chr"])
        if sub is None:
            return False
        for _, c in sub.iterrows():
            if not (row["end"] <= c["start"] or c["end"] <= row["start"]):
                return True
        return False

    mask = df.apply(_hits, axis=1)
    n_drop = int(mask.sum())
    return df[~mask].reset_index(drop=True), n_drop
