#!/usr/bin/env python3
"""Audit every notebook in the repository for unexecuted or failed cells.

Covers both Jupyter notebooks and the Mathematica notebook. Exits nonzero if any
code cell was never executed, any cell produced an error, or the Mathematica
notebook carries no evaluated output.

Run from the repository root:

    .venv/bin/python scripts/audit_notebooks.py 2>&1 | tee logs/audit_notebooks.log
"""

from __future__ import annotations

import json
import re
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parent.parent

JUPYTER = [
    "notebooks/gpt5_6_cosmology.ipynb",
    "notebooks/fable_spinor_dark_energy.ipynb",
]
MATHEMATICA = ["notebooks/gpt-5.6_bridge.nb"]


def audit_jupyter(name: str) -> tuple[bool, str]:
    notebook = json.loads((REPOSITORY_ROOT / name).read_text(encoding="utf-8"))
    cells = notebook["cells"]
    code = [c for c in cells if c["cell_type"] == "code"]
    interactive = [
        c for c in code if "interactive" in c.get("metadata", {}).get("tags", [])
    ]
    unrun = [
        c
        for c in code
        if c.get("execution_count") is None and c not in interactive
    ]
    errors = [
        o
        for c in cells
        for o in c.get("outputs", [])
        if o.get("output_type") == "error"
    ]
    figures = [
        o
        for c in cells
        for o in c.get("outputs", [])
        if "image/png" in o.get("data", {})
    ]
    produced = [c for c in code if c.get("outputs")]
    ok = not unrun and not errors
    detail = (
        f"{len(cells):>3} cells | {len(code):>2} code | "
        f"{len(code) - len(unrun):>2} executed | {len(unrun):>2} not executed | "
        f"{len(produced):>2} with output | {len(figures):>2} figures | "
        f"{len(errors):>2} errors"
    )
    if interactive:
        detail += f" | {len(interactive)} interactive (skipped by design)"
    return ok, detail


def audit_mathematica(name: str) -> tuple[bool, str]:
    text = (REPOSITORY_ROOT / name).read_text(encoding="utf-8", errors="replace")
    inputs = len(re.findall(r'"Input"', text))
    # An evaluated Wolfram cell yields Output (a returned value), Print (console
    # text) or Message. A cell body ending in ";" legitimately returns nothing,
    # so counting only "Output" would mark a correctly evaluated notebook unrun.
    returned = len(re.findall(r'"Output"', text))
    printed = len(re.findall(r'"Print"', text))
    messages = len(re.findall(r'"Message"', text))
    outputs = returned + printed
    cell_total = len(re.findall(r"\bCell\[", text))
    ok = inputs > 0 and outputs >= inputs and messages == 0
    detail = (
        f"{cell_total:>3} cells | {inputs:>2} Input | {outputs:>2} evaluated "
        f"output ({returned} Out, {printed} Print) | {messages:>2} messages | "
        f"{'evaluated' if ok else 'NOT EVALUATED'}"
    )
    return ok, detail


def main() -> int:
    failures: list[str] = []
    width = max(len(n) for n in JUPYTER + MATHEMATICA)

    print("Jupyter notebooks")
    for name in JUPYTER:
        ok, detail = audit_jupyter(name)
        print(f"  {'PASS' if ok else 'FAIL'}  {name:<{width}}  {detail}")
        if not ok:
            failures.append(name)

    print("\nMathematica notebooks")
    for name in MATHEMATICA:
        ok, detail = audit_mathematica(name)
        print(f"  {'PASS' if ok else 'FAIL'}  {name:<{width}}  {detail}")
        if not ok:
            failures.append(name)

    print()
    total = len(JUPYTER) + len(MATHEMATICA)
    print(f"notebooks audited : {total}")
    print(f"fully executed    : {total - len(failures)}")
    print(f"not executed      : {len(failures)}")
    if failures:
        print("\nFAILED: " + ", ".join(failures))
        return 1
    print("\nSUCCESS: every notebook is fully executed with zero errors.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
