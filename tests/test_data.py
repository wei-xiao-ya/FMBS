from pathlib import Path

import pandas as pd
import pytest

from fmbs.data import DataValidationError, load_drug_targets, load_query_smiles


def test_load_drug_targets_accepts_original_uniport_spelling(tmp_path: Path) -> None:
    path = tmp_path / "targets.csv"
    pd.DataFrame(
        {
            "smiles": ["CC", "CC", "CCC"],
            "drug": ["D1", "D1", "D2"],
            "uniport": ["P1", "P1", "P2"],
        }
    ).to_csv(path, index=False)

    assert load_drug_targets(path) == {"CC": ("P1",), "CCC": ("P2",)}


def test_load_query_smiles_requires_smiles_column(tmp_path: Path) -> None:
    path = tmp_path / "queries.csv"
    pd.DataFrame({"molecule": ["CC"]}).to_csv(path, index=False)

    with pytest.raises(DataValidationError, match="smiles"):
        load_query_smiles(path)
