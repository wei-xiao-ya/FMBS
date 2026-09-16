#!/usr/bin/env python3
"""Merge the 25 per-layer non-shared similarity CSV files.

The input files must be named ``nonshared_similarityA1.csv`` through
``nonshared_similarityE5.csv``.  Each file is expected to contain the same
ordered drug pairs and one similarity column.  The script validates those
assumptions while streaming, so it does not load the multi-gigabyte input set
into memory.
"""

from __future__ import annotations

import argparse
import csv
import io
import math
import os
from contextlib import ExitStack
from itertools import zip_longest
from pathlib import Path

LAYERS = tuple(f"{group}{index}" for group in "ABCDE" for index in range(1, 6))
PAIR_COLUMNS = ("drug 1", "drug 2")


def merge_similarity_files(input_dir: Path, output: Path) -> int:
    """Merge and validate layer files, returning the number of data rows."""
    input_paths = [input_dir / f"nonshared_similarity{layer}.csv" for layer in LAYERS]
    missing = [str(path) for path in input_paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing input files:\n" + "\n".join(missing))

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp")

    try:
        with ExitStack() as stack:
            handles = [stack.enter_context(path.open("rb")) for path in input_paths]

            for layer, handle in zip(LAYERS, handles):
                header_bytes = handle.readline()
                if not header_bytes:
                    raise ValueError(f"{layer}: input file is empty")
                header_text = header_bytes.decode("utf-8-sig").rstrip("\r\n")
                header = next(csv.reader(io.StringIO(header_text)))
                if len(header) < 3 or tuple(header[:2]) != PAIR_COLUMNS:
                    raise ValueError(
                        f"{layer}: expected first columns {PAIR_COLUMNS}, got {header!r}"
                    )

            output_handle = stack.enter_context(temporary.open("wb"))
            output_handle.write(
                (",".join((*PAIR_COLUMNS, *LAYERS)) + "\n").encode("utf-8")
            )

            row_count = 0
            for row_number, lines in enumerate(
                zip_longest(*(iter(handle.readline, b"") for handle in handles)), start=2
            ):
                if any(line is None for line in lines):
                    lengths = ["ended" if line is None else "present" for line in lines]
                    raise ValueError(
                        f"Layer files have different lengths at CSV row {row_number}: "
                        f"{dict(zip(LAYERS, lengths))}"
                    )

                pair_prefix: bytes | None = None
                values: list[bytes] = []
                for layer, line in zip(LAYERS, lines):
                    content = line.rstrip(b"\r\n")
                    pair, separator, value = content.rpartition(b",")
                    if not separator or not pair or not value:
                        raise ValueError(f"{layer}: incomplete data at CSV row {row_number}")
                    if pair_prefix is None:
                        pair_prefix = pair
                    elif pair != pair_prefix:
                        raise ValueError(
                            f"{layer}: drug-pair mismatch at CSV row {row_number}"
                        )
                    try:
                        similarity = float(value)
                    except ValueError as exc:
                        raise ValueError(
                            f"{layer}: invalid similarity {value!r} at CSV row {row_number}"
                        ) from exc
                    if not math.isfinite(similarity):
                        raise ValueError(
                            f"{layer}: non-finite similarity {value!r} at CSV row "
                            f"{row_number}"
                        )
                    values.append(value)

                output_handle.write(pair_prefix + b"," + b",".join(values) + b"\n")
                row_count += 1
                if row_count % 100_000 == 0:
                    print(f"Merged {row_count:,} rows", flush=True)

            output_handle.flush()
            os.fsync(output_handle.fileno())

        os.replace(temporary, output)
        return row_count
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Merge 25 non-shared similarity files into one validated CSV."
    )
    parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing nonshared_similarityA1.csv through E5.csv.",
    )
    parser.add_argument("output", type=Path, help="Destination nonshared_all.csv path.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    rows = merge_similarity_files(args.input_dir, args.output)
    print(f"Wrote {rows:,} rows to {args.output}")


if __name__ == "__main__":
    main()
