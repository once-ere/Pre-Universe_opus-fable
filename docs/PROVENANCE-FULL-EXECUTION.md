# Provenance — executing every notebook in this repository

**Effort.** This page records a single, specific action: finding that one of the
three notebooks in this repository had never actually been evaluated, fixing the
tooling that allowed that to happen, executing all three end to end, and making
the result a permanent, self-checking gate. Everything needed to repeat the work
is on this page. No other file needs to be consulted, and no command below is
abbreviated.

---

## 1. Environment actually used

| item | value |
|---|---|
| OS | Ubuntu 26.04.1 LTS, x86-64 |
| kernel | Linux 7.0.0-31-generic |
| C library | GNU libc 2.43 |
| shell | bash |
| Wolfram Language | 15.0.1 for Linux x86 (64-bit), 2 July 2026, Professional licence |
| Python (virtual environment) | 3.12.13, created with `uv` 0.11.26 |
| numpy / scipy / matplotlib | 2.5.3 / 1.18.1 / 3.11.2 |
| nbformat / nbconvert | 5.11.1 / via `jupyter nbconvert` |
| rustc / cargo | 1.96.1 (31fca3adb 2026-06-26) |
| pdfTeX | 3.141592653-2.6-1.40.28 (TeX Live 2025/Debian) |
| cores | 24 |
| git | 2.53.0 |
| SUNDIALS submodule commit | `c1be0885d1009357993e81b65e484c3a55743f79` |

---

## 2. What was wrong

Three notebooks live in this repository:

| notebook | kind |
|---|---|
| `notebooks/gpt5_6_cosmology.ipynb` | Jupyter |
| `notebooks/fable_spinor_dark_energy.ipynb` | Jupyter |
| `notebooks/gpt-5.6_bridge.nb` | Mathematica |

An audit of all three reported:

```
Jupyter notebooks
  PASS  notebooks/gpt5_6_cosmology.ipynb           24 cells |  9 code |  9 executed |  0 not executed |  9 with output |  0 figures |  0 errors
  PASS  notebooks/fable_spinor_dark_energy.ipynb   21 cells | 10 code | 10 executed |  0 not executed | 10 with output |  2 figures |  0 errors

Mathematica notebooks
  FAIL  notebooks/gpt-5.6_bridge.nb                15 cells |  5 Input |  0 Output | NOT EVALUATED

notebooks audited : 3
fully executed    : 2
not executed      : 1

FAILED: notebooks/gpt-5.6_bridge.nb
```

**Five Input cells, zero Output cells.** The Mathematica notebook shipped
unevaluated.

### Why it happened

`scripts/run_gpt56_notebook.wls` reads the notebook, evaluates each Input cell,
and checks the resulting test report. That is a genuine check — but it
**discards the results**. Nothing was ever written back, so the `.nb` file on
disk stayed exactly as `scripts/build_gpt56_notebook.py` had generated it:
input-only. The verification suite passed because it inspected the *log*, never
the notebook.

Two further defects surfaced while fixing this:

1. Cell 1 raised `ClearAll::clloc: Cannot clear local variable i.` The notebook
   opens with `ClearAll["Global`*"]`, and the executor's loop iterator `i` was a
   `Block`-scoped `Global`i`. Clearing a locally scoped symbol is an error.
2. Cells 1 to 4 each end in `;`, so they return `Null` and print nothing. Even a
   perfect executor would have stored only one output cell, which is thin
   evidence that anything was computed.

---

## 3. What was changed

### 3.1 A new in-place executor

`scripts/execute_gpt56_notebook.wls` evaluates each Input cell, captures both
the returned value and anything printed, and **writes the notebook back** with
an output cell after every input cell. It gates on messages: any Wolfram message
raised by any cell fails the run.

The iterator fix is one line, and the comment in the file says why it matters:

```wolfram
(* The iterator must NOT live in the Global context: cell one evaluates
   ClearAll["Global`*"], which would try to clear a Block-scoped Global`i and
   raise ClearAll::clloc. *)
Do[
  Module[{cell = GPT56Executor`cells[[GPT56Executor`index]], ...
```

### 3.2 Real per-section summaries

`wolfram/gpt56_bridge.wls` gained one summary `Print` at the end of each of its
first four sections. They report only quantities the section has already
computed, so the cost is negligible, and they make each notebook cell produce
genuine evidence instead of silence:

| cell | what it now prints |
|---|---|
| 1 | dimension, tangent signature, and `Det g` |
| 2 | gamma dimensions, that every entry is an integer, and the Clifford residual |
| 3 | that the Hadamard generator squares to the identity, `Det Lambda`, the rapidity |
| 4 | the interior torsion and curvature witnesses |
| 5 | the full 36-assertion verification report (already present) |

### 3.3 Interpreter metadata stripped from Jupyter notebooks

Executing a Jupyter notebook stamps `metadata.language_info.version` with
whichever CPython ran it. That one field made an otherwise byte-identical
notebook differ between the system interpreter (3.14.4) and the project virtual
environment (3.12.13) — a reproducibility failure with nothing to do with the
computation.

`scripts/normalize_notebooks.py` removes it. The field is optional in nbformat
4, so it is **removed rather than pinned to a fabricated value**. The script uses
nbformat's own reader and writer rather than plain `json`, so the output is
byte-identical to what `nbconvert` writes: same indentation, same separators,
and real UTF-8 rather than `\uXXXX` escapes. An earlier attempt using
`json.dumps` re-encoded every non-ASCII character and is the reason this
paragraph exists.

### 3.4 The gate itself

`scripts/verify_all.sh` grew from 12 stages to 13:

- **Stage 7** now calls `scripts/execute_gpt56_notebook.wls` instead of
  `scripts/run_gpt56_notebook.wls`, so the notebook is evaluated *in place*
  every time the suite runs. It greps for five exact lines, including
  `Output cells stored      : 5` and `Cells with messages      : 0`.
- **Stages 3 and 12** call `scripts/normalize_notebooks.py` immediately after
  each Jupyter execution.
- **Stage 13** is new: `scripts/audit_notebooks.py` walks every notebook,
  Jupyter and Mathematica alike, and fails if any code cell was never executed,
  any cell produced an error, or the Mathematica notebook carries fewer output
  cells than input cells.

`scripts/run_gpt56_notebook.wls` is retained for anyone who wants to evaluate
the notebook *without* rewriting the file.

---

## 4. Reproducing the whole thing

### 4.1 One-time setup

Install the Rust toolchain from <https://rustup.rs> if `cargo` is missing.
Install the Wolfram Engine or Mathematica 14 or newer, activated, if
`wolframscript` is missing.

```bash
git clone --recurse-submodules https://github.com/once-ere/Pre-Universe_opus-fable.git
cd Pre-Universe_opus-fable
mkdir -p logs build
```

If the repository was cloned without `--recurse-submodules`:

```bash
cd Pre-Universe_opus-fable
git submodule update --init --recursive
```

Create the Python environment. On Debian and Ubuntu the system Python refuses
system-wide `pip` installs, so a virtual environment is required:

```bash
cd Pre-Universe_opus-fable
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt
```

If `python3 -m venv` reports that `ensurepip` is unavailable, either install the
system package or use `uv`, which needs no administrator rights:

```bash
cd Pre-Universe_opus-fable
uv venv .venv
uv pip install --python .venv/bin/python -r requirements.txt
```

Build the Rust reference integrator:

```bash
cd Pre-Universe_opus-fable/fable_cosmo_rs
cargo build --release 2>&1 | tee ../logs/fable_cosmo_build.log
echo "exit status: ${PIPESTATUS[0]}"
cd ..
```

### 4.2 Audit before touching anything

```bash
cd Pre-Universe_opus-fable
.venv/bin/python scripts/audit_notebooks.py 2>&1 | tee logs/audit_notebooks.log
echo "exit status: ${PIPESTATUS[0]}"
```

### 4.3 Execute the Mathematica notebook in place

```bash
cd Pre-Universe_opus-fable
.venv/bin/python scripts/build_gpt56_notebook.py 2>&1 | tee logs/build_gpt56_notebook.log
wolframscript -file scripts/execute_gpt56_notebook.wls 2>&1 | tee logs/execute_gpt56_notebook.log
echo "exit status: ${PIPESTATUS[0]}"
```

### 4.4 Execute both Jupyter notebooks

```bash
cd Pre-Universe_opus-fable
.venv/bin/python scripts/build_cosmology_notebook.py
.venv/bin/python -m jupyter nbconvert --to notebook --execute --inplace \
  notebooks/gpt5_6_cosmology.ipynb \
  --ExecutePreprocessor.timeout=180 \
  --ExecutePreprocessor.record_timing=False 2>&1 | tee logs/execute_cosmology_notebook.log

.venv/bin/python scripts/build_fable_notebook.py
.venv/bin/python -m jupyter nbconvert --to notebook --execute --inplace \
  notebooks/fable_spinor_dark_energy.ipynb \
  --ExecutePreprocessor.timeout=900 \
  --ExecutePreprocessor.record_timing=False 2>&1 | tee logs/execute_fable_notebook.log

.venv/bin/python scripts/normalize_notebooks.py 2>&1 | tee logs/normalize_notebooks.log
```

### 4.5 Audit again

```bash
cd Pre-Universe_opus-fable
.venv/bin/python scripts/audit_notebooks.py 2>&1 | tee logs/audit_notebooks.log
echo "exit status: ${PIPESTATUS[0]}"
```

### 4.6 Or do all of it with one command

```bash
cd Pre-Universe_opus-fable
bash scripts/verify_all.sh 2>&1 | tee logs/verify_all_full.log
echo "exit status: ${PIPESTATUS[0]}"
```

---

## 5. The result

### 5.1 The Mathematica notebook, cell by cell

Full output of `logs/execute_gpt56_notebook.log` from the recorded run:

```
CELL 1 seconds=0.610606 printed=123 messages=0
CELL 2 seconds=0.128638 printed=123 messages=0
CELL 3 seconds=1.174485 printed=116 messages=0
CELL 4 seconds=20.345539 printed=178 messages=0
CELL 5 seconds=29.064920 printed=252 messages=0

Notebook path            : notebooks/gpt-5.6_bridge.nb
Input cells evaluated    : 5
Output cells stored      : 5
Total evaluation seconds : 51.324188
Cells with messages      : 0
Notebook tests succeeded : 36
Notebook tests failed    : 0
SUCCESS: the notebook is evaluated and its outputs are stored.
```

### 5.2 What the notebook now contains

The five stored output cells, read back out of the `.nb` file itself:

```
[1] canonical 4+4 geometry: dimension = 8, tangent signature =
    {1, 1, 1, 1, -1, -1, -1, -1}, Det g = Sec[6*H*x0]^2

[2] real Clifford algebra: gamma dimensions = {8, 16, 16},
    all entries integer = True, max |Clifford residual| = 0

[3] gpt-5.6f_rame: Hadamard generator squares to the identity = True,
    Det Lambda = 1, rapidity = kappa*x0*x4

[4] gpt-5.6_bridge family: interior torsion witness T^1_(0 1) = -1/2*(H*Cot[6*H*x0]),
    interior curvature witness R^1_(0 0 1) = (H^2*(25 + Cos[12*H*x0])*Csc[6*H*x0]^2)/4

Tests run: 36
Canonical nonzero spin-connection components: 24
Interior torsion witness T^1_(0 1): -1/2*(H*Cot[6*H*x0])
Interior curvature witness R^1_(0 0 1): (H^2*(25 + Cos[12*H*x0])*Csc[6*H*x0]^2)/4
Tests succeeded: 36
Tests failed: 0
```

These are computed values, not transcribed ones: `Det g = Sec[6Hx^0]^2` is the
determinant of the canonical 4+4 metric, the Clifford residual is exactly zero
over integer matrices, the Hadamard generator squares to the identity (which is
what makes the closed-form $\Lambda=\cosh u+\sinh u\,J$ possible), and both
interior witnesses are nonzero, which is the whole point of the bridge's
midpoint.

### 5.3 The audit, after

```
Jupyter notebooks
  PASS  notebooks/gpt5_6_cosmology.ipynb           24 cells |  9 code |  9 executed |  0 not executed |  9 with output |  0 figures |  0 errors
  PASS  notebooks/fable_spinor_dark_energy.ipynb   21 cells | 10 code | 10 executed |  0 not executed | 10 with output |  2 figures |  0 errors

Mathematica notebooks
  PASS  notebooks/gpt-5.6_bridge.nb                20 cells |  5 Input |  5 evaluated output (0 Out, 5 Print) |  0 messages | evaluated

notebooks audited : 3
fully executed    : 3
not executed      : 0

SUCCESS: every notebook is fully executed with zero errors.
```

`0 Out, 5 Print` is correct and expected: every cell body ends in `;`, so each
returns `Null` and its evidence is the printed text. `gpt5_6_cosmology.ipynb`
reports `0 figures` because it writes its plots to `artifacts/figures/` as files
rather than embedding them inline; its HTML rendering embeds them, and stage 4
of the suite asserts the alt text of all four.

### 5.4 Totals across the whole execution

| quantity | value |
|---|---|
| notebooks executed | 3 of 3 |
| notebooks not executed | 0 |
| Wolfram assertions, `wolfram/gpt56_bridge.wls` | 36 passed, 0 failed |
| Wolfram assertions, `wolfram/fable_spinor.wls` | 30 passed, 0 failed |
| Wolfram assertions, `wolfram/gpt56_bridge_dynamic.wls` | 17 passed, 0 failed |
| **Wolfram assertions, total** | **83 passed, 0 failed** |
| Python tests | 37 passed |
| Wolfram messages raised, any cell | 0 |
| Jupyter cell errors | 0 |
| verification stages | 13 of 13 passed |
| wall-clock, full suite | 2 min 16 s on 24 cores |

### 5.5 Determinism

Executing a notebook rewrites it, so the notebooks are only safe to pin in an
integrity manifest if execution is reproducible. `scripts/determinism_gate.py`
runs the entire pipeline twice and compares every watched artifact byte for
byte:

```
watched artifacts : 27
byte-identical    : 27
drifted           : 0

SUCCESS: every watched artifact is byte-identical across two full executions.
```

That includes all three executed notebooks, all ten figures and all three PDFs.

Run it yourself — it takes about two minutes:

```bash
cd Pre-Universe_opus-fable
.venv/bin/python scripts/determinism_gate.py 2>&1 | tee logs/determinism_gate.log
echo "exit status: ${PIPESTATUS[0]}"
```

### 5.6 Integrity manifest

With determinism established, every tracked file is pinned:

```bash
cd Pre-Universe_opus-fable
sha256sum -c artifacts/SHA256SUMS
```

Expected: `73` lines ending in `: OK` and no warning. The manifest is
regenerated only by `scripts/write_checksums.py`, and only after the determinism
gate passes — pinning an artifact that is not byte-reproducible would turn the
integrity check into a flaky test that people learn to ignore.

### 5.7 The full suite

Final line of `logs/verify_all_full.log`:

```
All repository verification gates passed.
```

Stage headings from the same run:

```
[1/13] Python tests
[2/13] Numerical artifacts
[3/13] Jupyter notebook build and execution
[4/13] Structured cosmology checks
[5/13] Markdown and Python hygiene
[6/13] LaTeX/PDF report
[7/13] Wolfram source and generated notebook
[8/13] Required deliverables
[9/13] fableSpinor symbolic proofs (Wolfram)
[10/13] fable_cosmo_rs: pure-Rust SUNDIALS reference integration
[11/13] Three-way cross-check and figures
[12/13] fableSpinor notebook, report and deliverables
[13/13] Every notebook executed, Jupyter and Mathematica
```

---

## 6. Viewing the executed notebooks in a web browser

### 6.1 The Jupyter notebooks, live

```bash
cd Pre-Universe_opus-fable
.venv/bin/python -m jupyter lab notebooks/
```

JupyterLab prints an address of the form `http://localhost:8888/lab?token=<hex>`.
Open it, choose a notebook, and use **Run > Run All Cells**. Stop the server with
`Ctrl-C` twice.

### 6.2 The Jupyter notebooks, as static pages

```bash
cd Pre-Universe_opus-fable
mkdir -p build
.venv/bin/python -m jupyter nbconvert --to html --output-dir build \
  notebooks/gpt5_6_cosmology.ipynb
.venv/bin/python -m jupyter nbconvert --to html --output-dir build \
  notebooks/fable_spinor_dark_energy.ipynb
python3 -m http.server 8911 --directory build
```

Open <http://127.0.0.1:8911/gpt5_6_cosmology.html> and
<http://127.0.0.1:8911/fable_spinor_dark_energy.html>. Stop with `Ctrl-C`.
Port 8911 is used deliberately: the related rustSolveIt project occupies
8895–8907.

### 6.3 The Mathematica notebook

On a machine with the Mathematica graphical front end:

```bash
cd Pre-Universe_opus-fable
mathematica notebooks/gpt-5.6_bridge.nb
```

The stored output cells are visible immediately, with no re-evaluation. To
recompute them in the front end, choose **Evaluation > Evaluate Notebook**;
expect about 51 seconds and the same 36 passing assertions.

### 6.4 This page and the figures

```bash
cd Pre-Universe_opus-fable
python3 -m http.server 8911
```

Open <http://127.0.0.1:8911/docs/PROVENANCE-FULL-EXECUTION.md> for this page
and <http://127.0.0.1:8911/artifacts/fable/figures/> for the figures. Stop with
`Ctrl-C`.

---

## 7. Honest limits of this page

- This page certifies that every notebook **ran to completion without errors or
  messages**, and that the results are byte-reproducible. It does not certify
  that the physics is correct; that is the business of the symbolic proofs and
  the three-way numerical cross-check, each documented on its own page.
- `0 figures` for `gpt5_6_cosmology.ipynb` is a property of how that notebook
  emits plots (to files, not inline). It is not a defect, and the HTML gate in
  stage 4 checks the figures another way.
- The Mathematica notebook stores `Print` cells rather than `Out` cells because
  every cell body ends in `;`. That is ordinary Wolfram behaviour, not a
  workaround.
- Determinism was established on one machine across two runs. Byte-identical
  output across *different* machines is not claimed; the interpreter-version
  normalisation removes the one cause of drift that was actually observed.
- Timings are from a 24-core machine and will differ elsewhere. Only the cell
  count, message count and assertion totals are invariant.
