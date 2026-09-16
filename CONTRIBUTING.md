# Contributing

Contributions that improve reproducibility, validation, documentation, or
performance are welcome.

## Development setup

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
ruff check .
pytest
```

Install `.[predict]` only when testing the Signaturizer integration; it includes
TensorFlow and is substantially larger than the core development environment.

## Pull requests

- Keep scientific-method changes separate from refactoring changes.
- Add or update tests for behavioral changes.
- Run `ruff check .` and `pytest` before opening a pull request.
- Do not commit generated outputs or the uncompressed non-shared table.
- Use Git LFS for `data/nonshared_all.csv.gz`.

Changes to bundled data should include updated row counts in `data/README.md`,
updated checksums in `data/SHA256SUMS`, and a regenerated model parameter file.
