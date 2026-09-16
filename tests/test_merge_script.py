import subprocess
import sys
from pathlib import Path

import pandas as pd

from fmbs.constants import LAYERS


def test_merge_similarity_files(tmp_path: Path) -> None:
    input_dir = tmp_path / "layers"
    input_dir.mkdir()
    pairs = pd.DataFrame({"drug 1": ["CC", "CCC"], "drug 2": ["CO", "CCO"]})
    for index, layer in enumerate(LAYERS):
        frame = pairs.assign(**{"Pearson Correlation": [index / 100, index / 200]})
        frame.to_csv(input_dir / f"nonshared_similarity{layer}.csv", index=False)

    output = tmp_path / "nonshared_all.csv"
    script = Path(__file__).parents[1] / "scripts" / "merge_similarity_files.py"
    subprocess.run(
        [sys.executable, str(script), str(input_dir), str(output)],
        check=True,
        capture_output=True,
        text=True,
    )

    merged = pd.read_csv(output)
    assert list(merged.columns) == ["drug 1", "drug 2", *LAYERS]
    assert len(merged) == 2
    assert merged.loc[0, "A1"] == 0.0
    assert merged.loc[0, "E5"] == 0.24
