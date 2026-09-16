"""Command-line interface for FMBS."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

from .constants import (
    DEFAULT_MODEL_PATH,
    DEFAULT_NONSHARED_PATH,
    DEFAULT_SHARED_PATH,
    DEFAULT_TARGETS_PATH,
)
from .data import load_drug_targets, load_query_smiles, validate_similarity_file
from .prediction import predict_targets
from .training import FMBSParameters, fit_parameters


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="fmbs",
        description="Fused multi-level bioactivity similarity target prediction.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    train = subparsers.add_parser("train", help="Fit the 25 likelihood functions.")
    train.add_argument("--shared", type=Path, default=DEFAULT_SHARED_PATH)
    train.add_argument("--nonshared", type=Path, default=DEFAULT_NONSHARED_PATH)
    train.add_argument("--output", type=Path, default=DEFAULT_MODEL_PATH)
    train.add_argument("--bins", type=int, default=19)
    train.add_argument("--chunksize", type=int, default=100_000)

    predict = subparsers.add_parser("predict", help="Predict targets for query SMILES.")
    predict.add_argument("input", type=Path, help="CSV containing a smiles column.")
    predict.add_argument("output", type=Path, help="Destination prediction CSV.")
    predict.add_argument("--targets", type=Path, default=DEFAULT_TARGETS_PATH)
    predict.add_argument("--model", type=Path, default=DEFAULT_MODEL_PATH)
    predict.add_argument("--threshold", type=float, default=1_000.0)
    predict.add_argument("--top-k", type=int, default=1)

    validate = subparsers.add_parser("validate-data", help="Validate bundled data files.")
    validate.add_argument("--shared", type=Path, default=DEFAULT_SHARED_PATH)
    validate.add_argument("--nonshared", type=Path, default=DEFAULT_NONSHARED_PATH)
    validate.add_argument("--targets", type=Path, default=DEFAULT_TARGETS_PATH)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "train":
        parameters = fit_parameters(
            args.shared,
            args.nonshared,
            bins=args.bins,
            chunksize=args.chunksize,
        )
        parameters.save(args.output)
        print(
            f"Fitted 25 layers from {parameters.shared_rows:,} shared and "
            f"{parameters.nonshared_rows:,} non-shared pairs."
        )
        print(f"Saved model parameters to {args.output}")
        return 0

    if args.command == "predict":
        queries = load_query_smiles(args.input)
        targets = load_drug_targets(args.targets)
        parameters = FMBSParameters.load(args.model)
        predictions = predict_targets(
            queries,
            targets,
            parameters,
            threshold=args.threshold,
            top_k=args.top_k,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        predictions.to_csv(args.output, index=False)
        print(f"Wrote {len(predictions):,} prediction rows to {args.output}")
        return 0

    if args.command == "validate-data":
        validate_similarity_file(args.shared)
        validate_similarity_file(args.nonshared)
        targets = load_drug_targets(args.targets)
        print(f"Similarity schemas are valid; loaded {len(targets):,} reference drugs.")
        return 0

    raise AssertionError(f"Unhandled command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
