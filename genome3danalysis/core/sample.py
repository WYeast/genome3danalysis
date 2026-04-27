"""Sample abstraction for structural-feature differential analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import h5py
import numpy as np
import pandas as pd


@dataclass
class Sample:
    """A single biological sample.

    Bulk mode
    ---------
    ``h5_paths`` contains one sf.h5 path per replicate.

    Single-cell mode
    ----------------
    ``h5_paths`` may contain one sf.h5 with multiple structures/cells. If the
    requested dataset is 2D, each row is treated as one cell. If the dataset
    is 1D, each file is treated as one cell.

    Parameters
    ----------
    name
        Display name, used as label prefix and in plots.
    h5_paths
        Paths to sf.h5 file(s).
    mode
        ``'bulk'`` or ``'single_cell'``.
    """

    name: str
    h5_paths: list[str]
    mode: Literal["bulk", "single_cell"] = "bulk"

    def load(
        self,
        feature_name: str,
        h5_dataset: str = "mean_arr",
        value_col_name: str | None = None,
    ) -> dict[str, pd.DataFrame]:
        """Load one feature from all configured h5 files.

        Parameters
        ----------
        feature_name
            H5 feature group (e.g., ``icp``, ``enhd``).
        h5_dataset
            Dataset name under ``/{feature_name}``.
        value_col_name
            Output value-column name. Defaults to ``feature_name``.

        Returns
        -------
        dict[str, pandas.DataFrame]
            In bulk mode: ``{f'{name}_{i+1}': df}``.
            In single-cell mode: ``{f'{name}_cell_{k:04d}': df}``.
            Each dataframe has columns ``['chr', 'start', 'end', value_col_name]``.
        """
        value_name = feature_name if value_col_name is None else value_col_name

        if self.mode == "bulk":
            out: dict[str, pd.DataFrame] = {}
            for i, fp in enumerate(self.h5_paths, start=1):
                out[f"{self.name}_{i}"] = _load_sample_h5(
                    fp,
                    feature_name=feature_name,
                    h5_dataset=h5_dataset,
                    value_col_name=value_name,
                )
            return out

        if self.mode != "single_cell":
            raise ValueError(f"Unknown mode: {self.mode!r}")

        out_sc: dict[str, pd.DataFrame] = {}
        cell_idx = 1
        for fp in self.h5_paths:
            with h5py.File(fp, "r") as h5:
                values, chr_h, start_h, end_h = _read_h5_arrays(h5, feature_name, h5_dataset)

            if values.ndim == 1:
                out_sc[f"{self.name}_cell_{cell_idx:04d}"] = pd.DataFrame(
                    {"chr": chr_h, "start": start_h, "end": end_h, value_name: values}
                )
                cell_idx += 1
                continue

            if values.ndim == 2:
                n_cells, n_bins = values.shape
                if n_bins != len(chr_h):
                    raise ValueError(
                        f"Length mismatch in {fp}: coords={len(chr_h)} but dataset bins={n_bins}"
                    )
                for r in range(n_cells):
                    out_sc[f"{self.name}_cell_{cell_idx:04d}"] = pd.DataFrame(
                        {
                            "chr": chr_h,
                            "start": start_h,
                            "end": end_h,
                            value_name: values[r, :],
                        }
                    )
                    cell_idx += 1
                continue

            raise NotImplementedError(
                f"TODO: unsupported dataset ndim={values.ndim} for single_cell in {fp}"
            )

        return out_sc


def _read_h5_arrays(
    h5: h5py.File,
    feature_name: str,
    h5_dataset: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if feature_name not in h5:
        raise KeyError(f"feature {feature_name!r} not in file; available: {list(h5.keys())}")
    feat_group = h5[feature_name]
    if h5_dataset not in feat_group:
        raise KeyError(
            f"dataset {h5_dataset!r} not in /{feature_name}; available: {list(feat_group.keys())}"
        )
    values = feat_group[h5_dataset][:]

    if "index" not in h5:
        raise KeyError("/index group missing")
    idx = h5["index"]
    for k in ["copy", "chromstr", "start", "end"]:
        if k not in idx:
            raise KeyError(f"/index/{k} missing")

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
    return values, chr_h, start_h, end_h


def _load_sample_h5(
    file_path: str,
    feature_name: str,
    h5_dataset: str = "mean_arr",
    value_col_name: str = "value",
) -> pd.DataFrame:
    """Read a sample's per-bin haploid feature value from sf.h5."""
    with h5py.File(file_path, "r") as h5:
        values, chr_h, start_h, end_h = _read_h5_arrays(h5, feature_name, h5_dataset)

        if values.ndim != 1:
            raise ValueError(
                f"Expected 1D dataset for bulk loading in {file_path}, got shape={values.shape}"
            )

        if len(chr_h) != len(values):
            raise ValueError(
                f"Length mismatch in {file_path}: "
                f"haploid coords (copy=={h5['index']['copy'][:].min()}) = {len(chr_h)}, "
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
