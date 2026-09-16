"""Input validation and tabular data helpers."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from .constants import LAYERS, PAIR_COLUMNS


class DataValidationError(ValueError):
    """Raised when an FMBS input file does not match the documented schema."""


def _column_lookup(columns: Sequence[str]) -> dict[str, str]:
    return {str(column).strip().lower(): str(column) for column in columns}


def validate_similarity_file(path: Path) -> None:
    """Validate the header of a shared or non-shared similarity file."""
    header = pd.read_csv(path, nrows=0)
    missing = [column for column in (*PAIR_COLUMNS, *LAYERS) if column not in header.columns]
    if missing:
        raise DataValidationError(f"{path} is missing columns: {', '.join(missing)}")


def load_drug_targets(path: Path) -> dict[str, tuple[str, ...]]:
    """Load a SMILES-to-UniProt mapping from the bundled reference table."""
    frame = pd.read_csv(path, dtype=str)
    lookup = _column_lookup(frame.columns)
    smiles_column = lookup.get("smiles")
    # Keep support for the misspelling in the original public data file.
    target_column = lookup.get("uniprot") or lookup.get("uniport") or lookup.get("target")
    if smiles_column is None or target_column is None:
        raise DataValidationError(
            f"{path} must contain 'smiles' and 'uniprot' columns "
            "('uniport' is accepted for backward compatibility)."
        )

    mapping: dict[str, list[str]] = defaultdict(list)
    subset = frame[[smiles_column, target_column]].dropna()
    for smiles, target in subset.itertuples(index=False, name=None):
        smiles = str(smiles).strip()
        target = str(target).strip()
        if smiles and target and target not in mapping[smiles]:
            mapping[smiles].append(target)

    if not mapping:
        raise DataValidationError(f"{path} contains no usable SMILES-to-target rows.")
    return {smiles: tuple(targets) for smiles, targets in mapping.items()}


def load_query_smiles(path: Path) -> list[str]:
    """Load non-empty SMILES values from a CSV file."""
    frame = pd.read_csv(path, dtype=str)
    lookup = _column_lookup(frame.columns)
    smiles_column = lookup.get("smiles")
    if smiles_column is None:
        raise DataValidationError(f"{path} must contain a 'smiles' column.")
    smiles = [value.strip() for value in frame[smiles_column].dropna() if value.strip()]
    if not smiles:
        raise DataValidationError(f"{path} contains no non-empty SMILES values.")
    return smiles
