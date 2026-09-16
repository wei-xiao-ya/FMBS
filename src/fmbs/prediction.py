"""Target prediction from Chemical Checker bioactivity signatures."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Callable, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from .constants import LAYERS
from .training import FMBSParameters

SignatureProvider = Callable[[str, Sequence[str]], np.ndarray]


def signaturizer_global_provider(smiles: Sequence[str]) -> dict[str, np.ndarray]:
    """Generate all 25 layers in one batched Signaturizer call."""
    try:
        from signaturizer import Signaturizer
    except ImportError as exc:
        raise RuntimeError(
            "Prediction requires the optional dependencies. Install them with "
            "`pip install -e '.[predict]'`."
        ) from exc

    result = Signaturizer(list(LAYERS), applicability=False).predict(list(smiles))
    signatures = np.asarray(result.signature, dtype=float)
    if signatures.ndim == 1:
        signatures = signatures.reshape(1, -1)
    if signatures.shape[0] != len(smiles):
        raise RuntimeError(
            f"Signaturizer returned {signatures.shape[0]} rows for "
            f"{len(smiles)} SMILES values."
        )
    expected_columns = len(LAYERS) * 128
    if signatures.shape[1] != expected_columns:
        raise RuntimeError(
            f"Signaturizer returned {signatures.shape[1]} columns; "
            f"expected {expected_columns}."
        )
    return {
        layer: signatures[:, index * 128 : (index + 1) * 128]
        for index, layer in enumerate(LAYERS)
    }


def _unique(values: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(values))


def _normalize_rows(values: np.ndarray, layer: str) -> np.ndarray:
    if values.ndim != 2 or values.shape[1] < 2:
        raise ValueError(f"{layer} signatures must be a two-dimensional matrix.")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{layer} contains a non-finite signature value.")
    centered = values - values.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(centered, axis=1, keepdims=True)
    if np.any(norms == 0):
        raise ValueError(f"{layer} contains a constant signature with undefined correlation.")
    return centered / norms


def predict_targets(
    query_smiles: Sequence[str],
    drug_targets: Mapping[str, Sequence[str]],
    parameters: FMBSParameters,
    *,
    threshold: float = 1_000.0,
    top_k: int = 1,
    signature_provider: SignatureProvider | None = None,
) -> pd.DataFrame:
    """Predict ranked UniProt targets for each query SMILES value."""
    if threshold <= 0:
        raise ValueError("threshold must be greater than zero")
    if top_k < 1:
        raise ValueError("top_k must be at least one")

    queries = list(query_smiles)
    references = list(drug_targets)
    if not queries or not references:
        raise ValueError("At least one query and one reference drug are required.")

    all_smiles = _unique([*references, *queries])
    positions = {smiles: index for index, smiles in enumerate(all_smiles)}
    query_indices = np.asarray([positions[smiles] for smiles in queries])
    reference_indices = np.asarray([positions[smiles] for smiles in references])

    log_scores = np.zeros((len(queries), len(references)), dtype=float)
    coefficients = parameters.coefficients(LAYERS)
    global_signatures = (
        signaturizer_global_provider(all_smiles) if signature_provider is None else None
    )
    for layer, coefficient in zip(LAYERS, coefficients):
        signatures = (
            global_signatures[layer]
            if global_signatures is not None
            else signature_provider(layer, all_smiles)
        )
        normalized = _normalize_rows(np.asarray(signatures, dtype=float), layer)
        correlations = normalized[query_indices] @ normalized[reference_indices].T
        log_scores += coefficient * np.clip(correlations, -1.0, 1.0)

    log_threshold = float(np.log(threshold))
    records: list[dict] = []
    for query_index, query in enumerate(queries):
        eligible = log_scores[query_index] > log_threshold
        for reference_index, reference in enumerate(references):
            if reference == query:
                eligible[reference_index] = False

        target_values: dict[str, list[float]] = defaultdict(list)
        eligible_indices = np.flatnonzero(eligible)
        for reference_index in eligible_indices:
            score = float(np.exp(min(log_scores[query_index, reference_index], 709.0)))
            for target in drug_targets[references[reference_index]]:
                target_values[target].append(score)

        ranked = sorted(
            (
                (target, float(np.mean(scores)), len(scores))
                for target, scores in target_values.items()
            ),
            key=lambda item: (-item[1], item[0]),
        )[:top_k]

        if not ranked:
            records.append(
                {
                    "smiles": query,
                    "rank": None,
                    "target": None,
                    "score": None,
                    "supporting_references": 0,
                }
            )
            continue
        for rank, (target, score, support) in enumerate(ranked, start=1):
            records.append(
                {
                    "smiles": query,
                    "rank": rank,
                    "target": target,
                    "score": score,
                    "supporting_references": support,
                }
            )

    return pd.DataFrame.from_records(records)
