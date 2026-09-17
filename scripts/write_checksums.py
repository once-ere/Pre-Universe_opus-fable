#!/usr/bin/env python3
"""Regenerate artifacts/SHA256SUMS over every tracked file.

The manifest records every file in the commit except itself and the submodule
gitlink. Only run this after scripts/determinism_gate.py passes: pinning an
artifact that is not byte-reproducible turns the integrity gate into a flaky
test that people learn to ignore.

Run from the repository root:

    .venv/bin/python scripts/write_checksums.py 2>&1 | tee logs/write_checksums.log
    sha256sum -c artifacts/SHA256SUMS
"""

from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parent.parent
MANIFEST = Path("artifacts/SHA256SUMS")
EXCLUDED = {MANIFEST.as_posix(), "vendor/sundials_rs"}


def tracked_files() -> list[str]:
    completed = subprocess.run(
        ["git", "ls-files"],
        cwd=REPOSITORY_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    names = []
    for name in completed.stdout.splitlines():
        if name in EXCLUDED:
            continue
        if not (REPOSITORY_ROOT / name).is_file():
            continue
        names.append(name)
    return sorted(names)


def main() -> int:
    names = tracked_files()
    lines = [
        f"{hashlib.sha256((REPOSITORY_ROOT / name).read_bytes()).hexdigest()}  {name}"
        for name in names
    ]
    (REPOSITORY_ROOT / MANIFEST).write_text(
        "\n".join(lines) + "\n", encoding="utf-8", newline="\n"
    )
    print(f"manifest : {MANIFEST}")
    print(f"entries  : {len(lines)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
