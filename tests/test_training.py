from pathlib import Path

import numpy as np
import pandas as pd

from fmbs.constants import LAYERS
from fmbs.training import FMBSParameters, fit_parameters


def _similarity_frame(values: list[float]) -> pd.DataFrame:
    frame = pd.DataFrame(
        {
            "drug 1": [f"left-{index}" for index in range(len(values))],
            "drug 2": [f"right-{index}" for index in range(len(values))],
        }
    )
    for layer in LAYERS:
        frame[layer] = values
    return frame


def test_fit_and_json_round_trip(tmp_path: Path) -> None:
    shared_values = [-0.75] + [-0.25] * 2 + [0.25] * 3 + [0.75] * 5
    nonshared_values = [-0.75] * 5 + [-0.25] * 3 + [0.25] * 2 + [0.75]
    shared = tmp_path / "shared.csv"
    nonshared = tmp_path / "nonshared.csv"
    model_path = tmp_path / "model.json"
    _similarity_frame(shared_values).to_csv(shared, index=False)
    _similarity_frame(nonshared_values).to_csv(nonshared, index=False)

    parameters = fit_parameters(shared, nonshared, bins=4, chunksize=3)
    parameters.save(model_path)
    restored = FMBSParameters.load(model_path)

    assert parameters.shared_rows == len(shared_values)
    assert parameters.nonshared_rows == len(nonshared_values)
    assert set(parameters.layers) == set(LAYERS)
    assert np.all(parameters.coefficients() > 0)
    np.testing.assert_allclose(restored.coefficients(), parameters.coefficients())
    assert restored.bin_edges == parameters.bin_edges
