"""Download published source datasets listed in the organism's source manifest.

Source data is deliberately not committed to this repository. This script fetches
it into ``<organism>/data/raw/`` and verifies the bytes against the SHA-256
recorded in ``sources.toml``.

The hash check is not ceremony. Our committed derived data and its published-total
assertions are only meaningful if they were produced from exactly these bytes; if
upstream silently revises a spreadsheet, we want a loud failure and a deliberate
re-review, not a quiet change in the connectome.

Every download prints the licence and citation of the file it fetched, so that
nobody ends up using this data without knowing whose work it is.
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import urllib.error
import urllib.request
from pathlib import Path

from worm.importers.sources import RAW_DIR, source_manifest

_ORGANISMS = {"worm": (source_manifest, RAW_DIR)}

CHUNK = 1 << 16


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while chunk := f.read(CHUNK):
            h.update(chunk)
    return h.hexdigest()


def _status(path: Path, expected: str) -> str:
    if not path.exists():
        return "missing"
    if not expected:
        return "present (hash unpinned)"
    return "ok" if sha256_file(path) == expected else "HASH MISMATCH"


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--organism", default="worm", choices=sorted(_ORGANISMS))
    p.add_argument("--dataset", action="append", help="source id; repeatable")
    p.add_argument("--all", action="store_true", help="fetch every entry in the manifest")
    p.add_argument("--check", action="store_true", help="verify what is on disk; download nothing")
    p.add_argument(
        "--print-hashes",
        action="store_true",
        help="print the SHA-256 of each present file, for populating sources.toml",
    )
    p.add_argument("--force", action="store_true", help="re-download even if the file is valid")
    args = p.parse_args(argv)

    manifest_fn, raw_dir = _ORGANISMS[args.organism]
    manifest = manifest_fn()

    if args.dataset:
        unknown = [d for d in args.dataset if d not in manifest]
        if unknown:
            print(f"unknown source id(s): {unknown}\nknown: {sorted(manifest)}", file=sys.stderr)
            return 2
        wanted = list(args.dataset)
    elif args.all or args.check or args.print_hashes:
        wanted = sorted(manifest)
    else:
        p.error("give --dataset ID, or --all, or --check")

    raw_dir.mkdir(parents=True, exist_ok=True)
    failures = 0

    for source_id in wanted:
        entry = manifest[source_id]
        filename = str(entry["filename"])
        expected = str(entry.get("sha256") or "")
        dest = raw_dir / filename

        if args.print_hashes:
            print(f"{source_id:32s} {sha256_file(dest) if dest.exists() else '<missing>'}")
            continue

        if args.check:
            state = _status(dest, expected)
            print(f"{source_id:32s} {state:24s} {filename}")
            if state == "HASH MISMATCH" or state == "missing":
                failures += 1
            continue

        ok_states = ("ok", "present (hash unpinned)")
        if dest.exists() and not args.force and _status(dest, expected) in ok_states:
            print(f"{source_id:32s} already present, hash ok")
            continue

        url = str(entry["url"])
        print(f"{source_id:32s} downloading {url}")
        try:
            with urllib.request.urlopen(url, timeout=120) as r:  # noqa: S310 - https from manifest
                data = r.read()
        except urllib.error.URLError as exc:
            print(f"  FAILED: {exc}", file=sys.stderr)
            failures += 1
            continue

        got = hashlib.sha256(data).hexdigest()
        if expected and got != expected:
            print(
                f"  HASH MISMATCH: expected {expected}\n"
                f"                 got      {got}\n"
                "  The upstream file changed. Not writing it. Review the change, then update "
                "worm/data/sources.toml and rebuild the derived data deliberately.",
                file=sys.stderr,
            )
            failures += 1
            continue

        dest.write_bytes(data)
        print(f"  wrote {dest.relative_to(raw_dir.parent.parent.parent)} ({len(data):,} bytes)")
        print(f"  licence:  {entry['license']}")
        print(f"  citation: {entry['citation']}")
        if entry.get("mirror_note"):
            print(f"  note:     {entry['mirror_note']}")

    if failures:
        print(f"\n{failures} problem(s).", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
