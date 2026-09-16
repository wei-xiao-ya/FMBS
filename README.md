# FMBS

FMBS (Fused Multi-level Bioactivity Similarity) is a research pipeline for
predicting protein targets from a molecule's SMILES representation. It combines
Pearson similarities from all 25 Chemical Checker Signaturizer spaces (A1–E5)
with likelihood ratios learned from drug pairs that do or do not share a known
target.

> **Research status:** this repository is an early research implementation. Its
> predictions are not validated for clinical use and must not be used for
> diagnosis or treatment decisions.

## How it works

1. Split known drug pairs into shared-target (ST) and non-shared-target
   (non-ST) groups.
2. For every Chemical Checker space, estimate the similarity distributions
   `P(s_i | ST)` and `P(s_i | non-ST)` with 19 bins over `[-1, 1]`.
3. Calculate the per-layer likelihood ratio and fit
   `L_i(s_i) = exp(a_i × s_i)`.
4. Generate A1–E5 signatures for a query and every reference drug, then compute
   their Pearson similarities.
5. Fuse the 25 fitted likelihoods as `L = product(L_i)`. Reference drugs above
   the threshold vote for their known UniProt targets; targets are ranked by
   their mean supporting score.

The refactored implementation preserves this model while fixing the missing
inputs and undefined variables in the original notebook-style script. It also
generates all 25 signature spaces in one batch and uses vectorized correlations.

## Repository layout

```text
FMBS/
├── data/
│   ├── drug_targets.csv          # reference SMILES → DrugBank → UniProt table
│   ├── shared_all.csv            # 11,740 shared-target pairs
│   ├── nonshared_all.csv.gz      # 648,935 non-shared-target pairs (Git LFS)
│   ├── examples/queries.csv      # example input
│   └── model/fmbs_model.json     # fitted likelihood parameters
├── scripts/
│   └── merge_similarity_files.py
├── src/fmbs/                     # reusable package and CLI
├── tests/                        # unit tests
├── environment.yml              # recommended Conda environment
└── pyproject.toml                # package metadata and dependencies
```

See [`data/README.md`](data/README.md) for the schemas and row counts.

## Installation

The Signaturizer 1.1.16 dependency pins TensorFlow 2.15.1, so use Python 3.9,
3.10, or 3.11. Conda is recommended because it provides RDKit consistently.

```bash
git lfs install
git clone https://github.com/wei-xiao-ya/FMBS.git
cd FMBS
conda env create -f environment.yml
conda activate fmbs
```

For development or model fitting without Signaturizer:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

To add prediction dependencies to an existing Python 3.9–3.11 environment:

```bash
python -m pip install -e ".[predict]"
```

The first Signaturizer run may download its pretrained model resources and can
take longer than subsequent runs.

On Apple Silicon, use a native `arm64` Python distribution (for example,
Miniforge). An Intel `x86_64` Python running through Rosetta can install the
packages but may crash when TensorFlow starts. Check with:

```bash
python -c "import platform; print(platform.machine())"
```

The output should be `arm64` on an Apple Silicon Mac.

## Validate the data

```bash
fmbs validate-data
```

This verifies that both similarity tables contain the two drug columns and all
25 A1–E5 similarity columns, and that the target mapping is readable.

## Fit the FMBS parameters

A fitted parameter file is included for prediction. To reproduce it from the
training tables:

```bash
fmbs train \
  --shared data/shared_all.csv \
  --nonshared data/nonshared_all.csv.gz \
  --output data/model/fmbs_model.json
```

Large files are processed in chunks, so fitting does not require loading the
complete 400 MB uncompressed table into memory.

## Predict targets

Create a CSV with one required column:

```csv
smiles
CC(=O)OC1=CC=CC=C1C(=O)O
```

Run:

```bash
fmbs predict data/examples/queries.csv outputs/predictions.csv
```

Useful options:

```bash
fmbs predict INPUT.csv OUTPUT.csv --threshold 1000 --top-k 5
```

The output columns are:

| Column | Meaning |
| --- | --- |
| `smiles` | Query molecule |
| `rank` | Target rank for the query |
| `target` | Predicted UniProt accession |
| `score` | Mean fused likelihood among supporting reference drugs |
| `supporting_references` | Number of reference drugs contributing to the target |

If no reference exceeds the threshold, the query is retained with empty target
and score fields.

## Rebuild the non-shared table

The source data were originally stored as 25 files named
`nonshared_similarityA1.csv` through `nonshared_similarityE5.csv`. Rebuild the
validated combined table with:

```bash
python scripts/merge_similarity_files.py PATH/TO/nonshared data/nonshared_all.csv
```

The script streams all inputs, confirms that every row contains the same drug
pair in every layer, validates numeric similarities, and writes the output
atomically.

## Development

```bash
ruff check .
pytest
```

GitHub Actions runs both commands on Python 3.9, 3.10, and 3.11.
See [`CONTRIBUTING.md`](CONTRIBUTING.md) before proposing changes.

## References

- Duran-Frigola M, et al. *Extending the small-molecule similarity principle to
  all levels of biology with the Chemical Checker.* Nature Biotechnology
  (2020). <https://doi.org/10.1038/s41587-020-0502-7>
- Bertoni M, et al. *Bioactivity descriptors for uncharacterized chemical
  compounds.* Nature Communications (2021).
  <https://doi.org/10.1038/s41467-021-24150-4>
- Signaturizer source code: <https://github.com/sbnb-irb/signaturizer>

If you use FMBS in published work, please cite the Chemical Checker and
Signaturizer papers above, together with the corresponding FMBS publication or
repository version.

## License

No software or data license has been declared yet. Add explicit licenses before
redistributing the code or datasets from this repository.
