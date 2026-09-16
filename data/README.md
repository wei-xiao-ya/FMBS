# Data files

| File | Purpose | Rows |
| --- | --- | ---: |
| `drug_targets.csv` | Reference SMILES, DrugBank identifiers, and UniProt targets | 3,199 |
| `shared_all.csv` | Drug pairs sharing at least one target, with A1–E5 similarities | 11,740 |
| `nonshared_all.csv.gz` | Drug pairs not sharing a target, with A1–E5 similarities | 648,935 |
| `model/fmbs_model.json` | Fitted layer coefficients and histogram provenance | generated |
| `examples/queries.csv` | Example prediction input | 2 |

The similarity tables use the columns `drug 1`, `drug 2`, followed by the 25
Chemical Checker spaces `A1` through `E5`. `shared_all.csv` also contains an
`ECFP4` column; FMBS does not use that extra column.

The uncompressed `nonshared_all.csv` is retained as a local working file and is
ignored by Git. Its compressed copy is managed by Git LFS because it exceeds
GitHub's normal per-file size limit.

The target column in `drug_targets.csv` is named `uniport` in the original
dataset. The loader accepts that historical spelling as well as `uniprot`.

From this directory, verify downloaded artifacts with:

```bash
shasum -a 256 -c SHA256SUMS
```
