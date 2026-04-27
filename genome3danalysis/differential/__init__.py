"""Differential structural-feature analysis public API."""

from .core import compute_differential_sf
from .io import build_matrix, filter_cen_overlap, load_sample_h5
from .plotting import (
    DEFAULT_COLORS,
    plot_delta_histogram,
    plot_diff_tracks,
    plot_manhattan,
    plot_md2_vs_delta,
    plot_qq,
    plot_replicate_scatter,
    plot_sample_correlation,
    plot_sf_tracks,
    plot_volcano,
)
from .qc import run_qc
from .stats import (
    compute_pairwise_fc,
    compute_pvalues,
    gate_significance,
    mahalanobis_robust,
    preprocess_matrix,
)

__all__ = [
    "compute_differential_sf",
    "run_qc",
    "load_sample_h5",
    "build_matrix",
    "filter_cen_overlap",
    "preprocess_matrix",
    "mahalanobis_robust",
    "compute_pvalues",
    "compute_pairwise_fc",
    "gate_significance",
    "DEFAULT_COLORS",
    "plot_sf_tracks",
    "plot_diff_tracks",
    "plot_manhattan",
    "plot_volcano",
    "plot_qq",
    "plot_replicate_scatter",
    "plot_md2_vs_delta",
    "plot_sample_correlation",
    "plot_delta_histogram",
]
