"""Fit FMBS layer-wise likelihood functions from similarity distributions."""

from __future__ import annotations

import json
import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from .constants import LAYERS
from .data import DataValidationError, validate_similarity_file


def exponential_likelihood(similarity: np.ndarray, coefficient: float) -> np.ndarray:
    """Evaluate the one-parameter likelihood curve used by the original FMBS code."""
    return np.exp(coefficient * similarity)


@dataclass(frozen=True)
class LayerParameters:
    coefficient: float
    shared_counts: tuple[int, ...]
    nonshared_counts: tuple[int, ...]
    likelihood_ratios: tuple[float | None, ...]

    def to_dict(self) -> dict:
        return {
            "coefficient": self.coefficient,
            "shared_counts": list(self.shared_counts),
            "nonshared_counts": list(self.nonshared_counts),
            "likelihood_ratios": list(self.likelihood_ratios),
        }


@dataclass(frozen=True)
class FMBSParameters:
    """Serializable parameters for all 25 Chemical Checker spaces."""

    bin_edges: tuple[float, ...]
    layers: Mapping[str, LayerParameters]
    shared_rows: int
    nonshared_rows: int
    created_at: str
    format_version: int = 1

    def coefficients(self, layers: Sequence[str] = LAYERS) -> np.ndarray:
        missing = [layer for layer in layers if layer not in self.layers]
        if missing:
            raise DataValidationError(
                "Parameter file is missing layers: " + ", ".join(missing)
            )
        return np.asarray([self.layers[layer].coefficient for layer in layers])

    def to_dict(self) -> dict:
        return {
            "format_version": self.format_version,
            "created_at": self.created_at,
            "method": "layer-wise likelihood ratio with exp(a * similarity) fit",
            "bin_edges": list(self.bin_edges),
            "shared_rows": self.shared_rows,
            "nonshared_rows": self.nonshared_rows,
            "layers": {layer: params.to_dict() for layer, params in self.layers.items()},
        }

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(f".{path.name}.tmp")
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(self.to_dict(), handle, indent=2, ensure_ascii=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)

    @classmethod
    def load(cls, path: Path) -> FMBSParameters:
        with path.open(encoding="utf-8") as handle:
            payload = json.load(handle)
        if payload.get("format_version") != 1:
            raise DataValidationError(
                f"Unsupported FMBS parameter format: {payload.get('format_version')!r}"
            )
        layers = {
            layer: LayerParameters(
                coefficient=float(values["coefficient"]),
                shared_counts=tuple(int(value) for value in values["shared_counts"]),
                nonshared_counts=tuple(
                    int(value) for value in values["nonshared_counts"]
                ),
                likelihood_ratios=tuple(
                    None if value is None else float(value)
                    for value in values["likelihood_ratios"]
                ),
            )
            for layer, values in payload["layers"].items()
        }
        return cls(
            bin_edges=tuple(float(value) for value in payload["bin_edges"]),
            layers=layers,
            shared_rows=int(payload["shared_rows"]),
            nonshared_rows=int(payload["nonshared_rows"]),
            created_at=str(payload["created_at"]),
            format_version=1,
        )


def _histograms(
    path: Path,
    layers: Sequence[str],
    bin_edges: np.ndarray,
    chunksize: int,
) -> tuple[dict[str, np.ndarray], dict[str, int], int]:
    validate_similarity_file(path)
    counts = {layer: np.zeros(len(bin_edges) - 1, dtype=np.int64) for layer in layers}
    totals = {layer: 0 for layer in layers}
    row_count = 0

    for chunk in pd.read_csv(path, usecols=list(layers), chunksize=chunksize):
        row_count += len(chunk)
        for layer in layers:
            values = pd.to_numeric(chunk[layer], errors="coerce").to_numpy(dtype=float)
            finite = values[np.isfinite(values)]
            totals[layer] += len(finite)
            counts[layer] += np.histogram(finite, bins=bin_edges)[0]

    if row_count == 0:
        raise DataValidationError(f"{path} contains no data rows.")
    return counts, totals, row_count


def fit_parameters(
    shared_path: Path,
    nonshared_path: Path,
    *,
    layers: Sequence[str] = LAYERS,
    bins: int = 19,
    chunksize: int = 100_000,
) -> FMBSParameters:
    """Fit FMBS coefficients while reading large CSVs in bounded-memory chunks."""
    if bins < 2:
        raise ValueError("bins must be at least 2")
    bin_edges = np.linspace(-1.0, 1.0, bins + 1)
    shared_counts, shared_totals, shared_rows = _histograms(
        shared_path, layers, bin_edges, chunksize
    )
    nonshared_counts, nonshared_totals, nonshared_rows = _histograms(
        nonshared_path, layers, bin_edges, chunksize
    )

    x_values = bin_edges[:-1]
    fitted: dict[str, LayerParameters] = {}
    for layer in layers:
        if shared_totals[layer] == 0 or nonshared_totals[layer] == 0:
            raise DataValidationError(f"{layer} contains no finite similarity values.")
        shared_probability = shared_counts[layer] / shared_totals[layer]
        nonshared_probability = nonshared_counts[layer] / nonshared_totals[layer]
        with np.errstate(divide="ignore", invalid="ignore"):
            ratios = shared_probability / nonshared_probability
        valid = np.isfinite(ratios) & (nonshared_counts[layer] > 0)
        if valid.sum() < 2:
            raise DataValidationError(
                f"{layer} has fewer than two usable histogram bins for fitting."
            )
        coefficient, _ = curve_fit(
            exponential_likelihood,
            x_values[valid],
            ratios[valid],
            p0=(1.0,),
            maxfev=20_000,
        )
        serializable_ratios = tuple(
            float(value) if np.isfinite(value) else None for value in ratios
        )
        fitted[layer] = LayerParameters(
            coefficient=float(coefficient[0]),
            shared_counts=tuple(int(value) for value in shared_counts[layer]),
            nonshared_counts=tuple(int(value) for value in nonshared_counts[layer]),
            likelihood_ratios=serializable_ratios,
        )

    return FMBSParameters(
        bin_edges=tuple(float(value) for value in bin_edges),
        layers=fitted,
        shared_rows=shared_rows,
        nonshared_rows=nonshared_rows,
        created_at=datetime.now(timezone.utc).isoformat(),
    )
