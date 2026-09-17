#!/usr/bin/env python3
"""Strip interpreter-dependent metadata from executed Jupyter notebooks.

Executing a notebook stamps `metadata.language_info.version` with whichever
CPython happened to run it. That single field made an otherwise byte-identical
notebook differ between a system interpreter and a virtual environment, which
would defeat the determinism gate and the checksum manifest for a reason that
has nothing to do with the computation.

The field is optional in nbformat 4, so it is removed rather than pinned to a
fabricated value. Everything that records what was actually computed --- the
cell sources, the execution counts and the outputs --- is left untouched.

Run from the repository root:

    .venv/bin/python scripts/normalize_notebooks.py 2>&1 | tee logs/normalize_notebooks.log
"""

from __future__ import annotations

from pathlib import Path

import nbformat

REPOSITORY_ROOT = Path(__file__).resolve().parent.parent

NOTEBOOKS = [
    "notebooks/gpt5_6_cosmology.ipynb",
    "notebooks/fable_spinor_dark_energy.ipynb",
]

#: Metadata paths removed because they record the machine, not the result.
VOLATILE = [("language_info", "version")]


def normalize(path: Path) -> list[str]:
    # nbformat's own reader/writer is used rather than plain json so the file is
    # byte-identical to what nbconvert writes: same indent, same separators, and
    # real UTF-8 rather than \uXXXX escapes.
    notebook = nbformat.read(path, as_version=nbformat.NO_CONVERT)
    removed: list[str] = []
    for section, key in VOLATILE:
        block = notebook.get("metadata", {}).get(section)
        if isinstance(block, dict) and key in block:
            removed.append(f"{section}.{key}={block.pop(key)!r}")
    before = path.read_bytes()
    nbformat.write(notebook, path)
    if path.read_bytes() == before and not removed:
        return []
    return removed or ["reformatted to nbformat canonical form"]


def main() -> int:
    total = 0
    for name in NOTEBOOKS:
        path = REPOSITORY_ROOT / name
        if not path.is_file():
            print(f"  skip  {name} (absent)")
            continue
        removed = normalize(path)
        total += len(removed)
        detail = ", ".join(removed) if removed else "already normalized"
        print(f"  ok    {name}  [{detail}]")
    print()
    print(f"notebooks normalized : {len(NOTEBOOKS)}")
    print(f"volatile fields removed : {total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
