"""Shared constants for the FMBS pipeline."""

from pathlib import Path

LAYERS = tuple(f"{group}{index}" for group in "ABCDE" for index in range(1, 6))
PAIR_COLUMNS = ("drug 1", "drug 2")

DEFAULT_SHARED_PATH = Path("data/shared_all.csv")
DEFAULT_NONSHARED_PATH = Path("data/nonshared_all.csv.gz")
DEFAULT_TARGETS_PATH = Path("data/drug_targets.csv")
DEFAULT_MODEL_PATH = Path("data/model/fmbs_model.json")
