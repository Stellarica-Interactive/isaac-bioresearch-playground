"""Regenerate the committed normalized datasets from the raw published files.

This is the only thing that writes ``worm/data/normalized/``. Run it after
changing an importer, then review the ``git diff``: because the on-disk format is
sorted CSV, a parser change shows up as "these three edges changed weight" rather
than as an unreadable blob.

Only *anatomy* is written. Role, neurotransmitter and class annotations are
applied at load time, so a corrected annotation table never dirties these files.
"""

from __future__ import annotations

import argparse
import sys

from common.data.io import save_connectome
from common.data.registry import list_datasets
from common.data.validate import validate
from worm.importers.sources import compare_to_reference, has_reference, normalized_path
from worm.loader import parse_raw


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--organism", default="worm", choices=["worm"])
    p.add_argument("--dataset", action="append", help="dataset id; repeatable. Default: all.")
    p.add_argument(
        "--allow-unverified",
        action="store_true",
        help="write datasets that have no published totals to check against",
    )
    args = p.parse_args(argv)

    wanted = args.dataset or list_datasets()
    failures = 0

    for dataset_id in wanted:
        try:
            c = parse_raw(dataset_id)
        except FileNotFoundError as exc:
            print(f"{dataset_id:20s} SKIP  {exc}".replace("\n", " "), file=sys.stderr)
            failures += 1
            continue

        problems = validate(c)
        if problems:
            print(f"{dataset_id:20s} INVALID", file=sys.stderr)
            for v in problems[:10]:
                print(f"  {v}", file=sys.stderr)
            failures += 1
            continue

        totals = c.totals().to_dict()
        if has_reference(dataset_id):
            mismatches = compare_to_reference(dataset_id, totals)
            if mismatches:
                print(f"{dataset_id:20s} TOTALS DO NOT MATCH PUBLISHED FIGURES", file=sys.stderr)
                for m in mismatches:
                    print(f"  {m}", file=sys.stderr)
                failures += 1
                continue
            verdict = "verified"
        elif args.allow_unverified:
            verdict = "UNVERIFIED (no published totals)"
        else:
            print(
                f"{dataset_id:20s} has no entry in reference_totals.toml; "
                "pass --allow-unverified to write it anyway",
                file=sys.stderr,
            )
            failures += 1
            continue

        dest = normalized_path(dataset_id)
        save_connectome(c, dest)
        print(
            f"{dataset_id:20s} {verdict:32s} "
            f"cells={totals['cells']:3d} chem={totals['chemical_edges']:4d} "
            f"elec={totals['electrical_edges_undirected']:3d} -> {dest.name}/"
        )

    if failures:
        print(f"\n{failures} dataset(s) not written.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
