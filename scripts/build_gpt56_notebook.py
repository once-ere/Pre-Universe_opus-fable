#!/usr/bin/env python3
"""Build the self-contained gpt-5.6_bridge Mathematica notebook."""

from __future__ import annotations

from pathlib import Path
import re


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATH = REPOSITORY_ROOT / "wolfram" / "gpt56_bridge.wls"
NOTEBOOK_PATH = REPOSITORY_ROOT / "notebooks" / "gpt-5.6_bridge.nb"
SECTION_MARKER = re.compile(r"^\(\* ::NotebookSection:: (.+?) \*\)$")


def wolfram_string(value: str) -> str:
    replacements = {
        "\\": "\\\\",
        '"': '\\"',
        "\n": "\\n",
        "\t": "\\t",
    }
    return '"' + "".join(replacements.get(character, character) for character in value) + '"'


def source_sections(source: str) -> list[tuple[str, str]]:
    sections: list[tuple[str, str]] = []
    title: str | None = None
    lines: list[str] = []

    for line in source.splitlines():
        marker = SECTION_MARKER.match(line)
        if marker is None:
            lines.append(line)
            continue
        if title is not None:
            sections.append((title, "\n".join(lines).strip() + "\n"))
            lines = []
        title = marker.group(1)

    if title is None:
        raise ValueError(f"No notebook section markers found in {SOURCE_PATH}")
    sections.append((title, "\n".join(lines).strip() + "\n"))
    return sections


def notebook_cell(style: str, body: str) -> str:
    if style == "Input":
        return f'Cell[BoxData[{wolfram_string(body)}], "Input"]'
    return f"Cell[{wolfram_string(body)}, {wolfram_string(style)}]"


def build_notebook() -> tuple[int, int]:
    source = SOURCE_PATH.read_text(encoding="utf-8")
    sections = source_sections(source)
    cells = [
        notebook_cell("Title", "gpt-5.6_bridge"),
        notebook_cell(
            "Subtitle",
            "A curvature-torsion bridge in a local 4+4 pseudo-orthogonal frame",
        ),
        notebook_cell(
            "Text",
            "The visible frame name is gpt-5.6f_rame. Internal Wolfram symbols avoid punctuation: "
            "gpt56Frame, gpt56Bridge, and omegaGPT56. The frame is a local SO(4,4) gauge choice; "
            "the new geometric object is the metric-compatible connection family whose midpoint "
            "has both nonzero torsion and nonzero curvature.",
        ),
    ]
    for title, code in sections:
        cells.append(notebook_cell("Section", title))
        cells.append(notebook_cell("Input", code))
    cells.extend(
        [
            notebook_cell("Section", "Interpretation"),
            notebook_cell(
                "Text",
                "At t=0 the family is Levi-Civita. At t=1 it is the curvature-free "
                "Weitzenboeck connection used by FABLE-5.1. The distinguished new connection "
                "omegaGPT56 is t=1/2. Its value is not claimed to be experimentally preferred; "
                "its advantage is that one covariant, exactly tested family contains both endpoint "
                "geometries and a nontrivial mixed curvature-torsion interior.",
            ),
        ]
    )

    notebook = "\n".join(
        [
            "(* Content-type: application/vnd.wolfram.mathematica *)",
            "",
            "(*** Wolfram Notebook File ***)",
            "(* http://www.wolfram.com/nb *)",
            "",
            "(* Generated deterministically by scripts/build_gpt56_notebook.py. *)",
            "",
            "Notebook[{",
            ",\n".join(cells),
            "},",
            "WindowSize->{1400, 900},",
            f"WindowTitle->{wolfram_string('gpt-5.6_bridge')},",
            "WindowMargins->{{Automatic, 0}, {Automatic, 0}},",
            'StyleDefinitions->"Default.nb"',
            "]",
            "",
        ]
    )
    NOTEBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    NOTEBOOK_PATH.write_text(notebook, encoding="utf-8", newline="\n")
    return len(cells), len(sections)


def main() -> int:
    cell_count, section_count = build_notebook()
    print(f"source: {SOURCE_PATH.relative_to(REPOSITORY_ROOT)}")
    print(f"sections: {section_count}")
    print(f"cells: {cell_count}")
    print(f"notebook: {NOTEBOOK_PATH.relative_to(REPOSITORY_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
