"""One entry point for getting an annotated *C. elegans* connectome.

    from worm.loader import load
    connectome, reports = load("witvliet_2021_7")

Resolution order: prefer the committed normalized data in
``worm/data/normalized/``, fall back to parsing the raw source file if it is
present. So a fresh clone works offline, and a developer who has fetched raw data
can check that the committed files still reproduce.

Annotations are applied here rather than being baked into the stored files. The
default set is anatomy-adjacent and well-published (neuron class, coarse role,
neurotransmitter). Synaptic polarity is deliberately not in the default set —
see :mod:`worm.annotations.overlays`.
"""

from __future__ import annotations

from common.data.io import load_connectome
from common.data.overlay import OverlayReport, apply_overlays
from common.data.registry import get_importer
from common.data.schemas import Connectome
from worm.annotations.overlays import DEFAULT_OVERLAYS, get_overlays
from worm.importers import cook_2019 as _cook  # noqa: F401  (registers importer)
from worm.importers import witvliet_2021 as _witvliet  # noqa: F401  (registers importers)
from worm.importers.sources import normalized_path, raw_path


def parse_raw(dataset_id: str) -> Connectome:
    """Parse the raw published file. Requires ``tools/fetch_datasets.py`` to have run."""
    importer = get_importer(dataset_id)
    files = {spec.key: raw_path(dataset_id) for spec in importer.source_files()}
    missing = [str(p) for p in files.values() if not p.exists()]
    if missing:
        raise FileNotFoundError(
            f"raw source file(s) not present for {dataset_id}: {missing}\n"
            f"Run: python tools/fetch_datasets.py --organism worm --dataset {dataset_id}"
        )
    return importer.parse(files)


def load_anatomy(dataset_id: str, *, prefer: str = "normalized") -> Connectome:
    """Unannotated connectome, from committed data or from the raw file."""
    normalized = normalized_path(dataset_id)
    if prefer == "normalized" and (normalized / "meta.json").exists():
        return load_connectome(normalized)
    if prefer == "raw":
        return parse_raw(dataset_id)
    if (normalized / "meta.json").exists():
        return load_connectome(normalized)
    return parse_raw(dataset_id)


def load(
    dataset_id: str,
    *,
    annotations: tuple[str, ...] = DEFAULT_OVERLAYS,
    prefer: str = "normalized",
) -> tuple[Connectome, list[OverlayReport]]:
    """Annotated connectome plus a report of what each overlay could and could not cover."""
    c = load_anatomy(dataset_id, prefer=prefer)
    if not annotations:
        return c, []
    return apply_overlays(c, get_overlays(annotations))
