"""Deterministic on-disk format for normalized connectomes.

Layout, one directory per dataset::

    worm/data/normalized/witvliet_2021_7/
        cells.csv        id,category,class_name,source_labels
        connections.csv  pre,post,synapse_type,weight
        meta.json        provenance, sources, schema_version, totals, csv hashes

Why CSV plus a JSON sidecar rather than one JSON blob or Parquet:

* These graphs are small (hundreds of nodes, thousands of edges), so a binary
  columnar format buys nothing while costing ``git diff``, ``git blame`` on an
  individual edge, and reviewability of a parser change in a pull request.
* Sorted CSV means a changed synapse count churns exactly one line. A JSON array
  of objects churns four to six re-indented lines for the same change.
* Provenance is nested and heterogeneous and is read once, so it belongs in JSON.

Only *unannotated* anatomy is stored. Role, neurotransmitter and polarity
annotations are applied at load time by :mod:`common.data.overlay`, so a
correction to an annotation table never dirties a connectome's golden files.

Determinism rules, all of which are asserted by the test suite:

* newline is always ``\\n`` (a CRLF checkout would break every recorded hash);
* rows are sorted by a fixed key;
* column order is fixed;
* JSON is written with sorted keys and a trailing newline.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
from collections.abc import Mapping, Sequence
from pathlib import Path

from common.data.schemas import (
    SCHEMA_VERSION,
    Cell,
    CellCategory,
    Confidence,
    Connection,
    Connectome,
    Provenance,
    Scope,
    Sign,
    SynapseType,
)

CELLS_FILE = "cells.csv"
CONNECTIONS_FILE = "connections.csv"
META_FILE = "meta.json"

# `category_source` is carried explicitly because a cell's category (neuron /
# muscle / glia / other) is NOT measured by a connectome reconstruction — it comes
# from a cell registry. Dropping the citation on the way to disk would silently
# turn an annotation into an apparent measurement.
CELL_COLUMNS = ("id", "category", "category_source", "source_labels")
CONNECTION_COLUMNS = ("pre", "post", "synapse_type", "weight")

MULTI_SEP = "|"


# ---------------------------------------------------------------------------
# Rendering (pure string production, so tests can compare bytes without disk I/O)
# ---------------------------------------------------------------------------


def _write_csv(columns: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    buf = io.StringIO(newline="")
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(columns)
    w.writerows(rows)
    return buf.getvalue()


def render_cells_csv(c: Connectome) -> str:
    rows = [
        (
            cell.id,
            str(cell.category),
            cell.field_sources.get("category", ""),
            MULTI_SEP.join(sorted(cell.source_labels)),
        )
        for cell in sorted(c.cells, key=lambda x: x.id)
    ]
    return _write_csv(CELL_COLUMNS, rows)


def render_connections_csv(c: Connectome) -> str:
    rows = [
        (e.pre, e.post, str(e.synapse_type), str(e.weight))
        for e in sorted(c.connections, key=lambda e: (str(e.synapse_type), e.pre, e.post))
    ]
    return _write_csv(CONNECTION_COLUMNS, rows)


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def render_meta_json(c: Connectome, cells_csv: str, connections_csv: str) -> str:
    meta = {
        "id": c.id,
        "schema_version": c.schema_version,
        "organism": c.organism,
        "sex": c.sex,
        "stage": c.stage,
        "scope": str(c.scope),
        "provenance": c.provenance.to_dict(),
        "sources": {k: v.to_dict() for k, v in sorted(c.sources.items())},
        "totals": c.totals().to_dict(),
        "files": {
            CELLS_FILE: {"sha256": sha256_text(cells_csv)},
            CONNECTIONS_FILE: {"sha256": sha256_text(connections_csv)},
        },
    }
    return json.dumps(meta, indent=2, sort_keys=True, ensure_ascii=False) + "\n"


def render(c: Connectome) -> dict[str, str]:
    """Full on-disk representation as ``{filename: text}``. No I/O."""
    cells_csv = render_cells_csv(c)
    connections_csv = render_connections_csv(c)
    return {
        CELLS_FILE: cells_csv,
        CONNECTIONS_FILE: connections_csv,
        META_FILE: render_meta_json(c, cells_csv, connections_csv),
    }


# ---------------------------------------------------------------------------
# Disk I/O
# ---------------------------------------------------------------------------


def save_connectome(c: Connectome, dest: Path) -> dict[str, str]:
    """Write ``c`` into directory ``dest``. Returns the rendered files."""
    dest.mkdir(parents=True, exist_ok=True)
    files = render(c)
    for name, text in files.items():
        (dest / name).write_text(text, encoding="utf-8", newline="")
    return files


def load_connectome(src: Path) -> Connectome:
    """Read a normalized dataset directory written by :func:`save_connectome`.

    Recorded CSV hashes are verified: silent corruption of committed derived data
    is exactly the failure mode this format exists to prevent.
    """
    meta_path = src / META_FILE
    if not meta_path.exists():
        raise FileNotFoundError(f"not a normalized dataset directory (no {META_FILE}): {src}")
    meta = json.loads(meta_path.read_text(encoding="utf-8"))

    version = str(meta.get("schema_version"))
    if version != SCHEMA_VERSION:
        raise ValueError(
            f"{src} was written by schema version {version!r}, this build expects "
            f"{SCHEMA_VERSION!r}; re-run tools/build_datasets.py"
        )

    cells_csv = (src / CELLS_FILE).read_text(encoding="utf-8")
    connections_csv = (src / CONNECTIONS_FILE).read_text(encoding="utf-8")
    for name, text in ((CELLS_FILE, cells_csv), (CONNECTIONS_FILE, connections_csv)):
        want = meta.get("files", {}).get(name, {}).get("sha256")
        got = sha256_text(text)
        if want and want != got:
            raise ValueError(
                f"{src / name} does not match the hash recorded in {META_FILE} "
                f"(expected {want[:12]}…, got {got[:12]}…). The file was edited by hand, "
                "or checked out with CRLF line endings — see .gitattributes."
            )

    sources = {
        k: Provenance.from_dict(v) for k, v in (meta.get("sources") or {}).items()
    }
    provenance = Provenance.from_dict(meta["provenance"])

    cells = tuple(
        Cell(
            id=row["id"],
            category=CellCategory(row["category"]),
            # class_name is deliberately not stored here: it is an annotation
            # (Wang et al. 2024), applied at load time by an overlay.
            class_name=None,
            source_labels=tuple(x for x in row["source_labels"].split(MULTI_SEP) if x),
            field_sources=(
                {"category": row["category_source"]} if row["category_source"] else {}
            ),
        )
        for row in _read_csv(cells_csv, CELL_COLUMNS, src / CELLS_FILE)
    )
    connections = tuple(
        Connection(
            pre=row["pre"],
            post=row["post"],
            synapse_type=SynapseType(row["synapse_type"]),
            weight=int(row["weight"]),
            sign=Sign.UNKNOWN,
            sign_confidence=Confidence.ASSUMED,
            field_sources={},
        )
        for row in _read_csv(connections_csv, CONNECTION_COLUMNS, src / CONNECTIONS_FILE)
    )

    return Connectome(
        id=str(meta["id"]),
        organism=str(meta["organism"]),
        sex=str(meta["sex"]),
        stage=str(meta["stage"]),
        scope=Scope(meta["scope"]),
        cells=cells,
        connections=connections,
        provenance=provenance,
        sources=sources,
        schema_version=version,
    )


def _read_csv(text: str, expected: Sequence[str], path: Path) -> list[Mapping[str, str]]:
    reader = csv.DictReader(io.StringIO(text, newline=""))
    if tuple(reader.fieldnames or ()) != tuple(expected):
        raise ValueError(
            f"{path}: expected columns {list(expected)}, found {reader.fieldnames}"
        )
    return list(reader)
