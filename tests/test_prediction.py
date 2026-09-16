from datetime import datetime, timezone

import numpy as np

from fmbs.constants import LAYERS
from fmbs.prediction import predict_targets
from fmbs.training import FMBSParameters, LayerParameters


def _parameters() -> FMBSParameters:
    empty_counts = (1, 1)
    layers = {
        layer: LayerParameters(
            coefficient=10.0 if layer == "A1" else 0.0,
            shared_counts=empty_counts,
            nonshared_counts=empty_counts,
            likelihood_ratios=(1.0, 1.0),
        )
        for layer in LAYERS
    }
    return FMBSParameters(
        bin_edges=(-1.0, 0.0, 1.0),
        layers=layers,
        shared_rows=2,
        nonshared_rows=2,
        created_at=datetime.now(timezone.utc).isoformat(),
    )


def _signature_provider(layer: str, smiles: list[str]) -> np.ndarray:
    vectors = {
        "query": [1.0, 2.0, 3.0],
        "similar": [2.0, 4.0, 6.0],
        "opposite": [3.0, 2.0, 1.0],
    }
    return np.asarray([vectors[value] for value in smiles])


def test_predict_targets_ranks_supported_target() -> None:
    result = predict_targets(
        ["query"],
        {"similar": ("P12345",), "opposite": ("P99999",)},
        _parameters(),
        threshold=1_000.0,
        signature_provider=_signature_provider,
    )

    assert result.loc[0, "target"] == "P12345"
    assert result.loc[0, "rank"] == 1
    assert result.loc[0, "supporting_references"] == 1
    assert result.loc[0, "score"] > 1_000


def test_predict_targets_excludes_identical_reference() -> None:
    result = predict_targets(
        ["similar"],
        {"similar": ("P12345",), "opposite": ("P99999",)},
        _parameters(),
        threshold=1_000.0,
        signature_provider=_signature_provider,
    )

    assert result.loc[0, "target"] is None
    assert result.loc[0, "supporting_references"] == 0
