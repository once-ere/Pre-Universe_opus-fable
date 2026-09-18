# Provenance — the fableSpinor notebook (`notebooks/fable_spinor_dark_energy.ipynb`)

**Effort.** This page records how the executable notebook is authored, executed,
checked and displayed in a web browser, and exactly what it proves. Everything
needed to repeat the work is on this page. No other file needs to be consulted,
and no command below is abbreviated.

---

## 1. Environment actually used

| item | value |
|---|---|
| OS | Ubuntu 26.04.1 LTS, x86-64 |
| kernel | Linux 7.0.0-31-generic |
| C library | GNU libc 2.43 |
| shell | bash |
| Python | 3.14.4 |
| numpy / scipy / matplotlib | 2.5.3 / 1.18.1 / 3.11.2 |
| nbformat / IPython / ipykernel | 5.11.1 / 9.17.1 / 7.3.0 |
| rustc / cargo | 1.96.1 (31fca3adb 2026-06-26) |
| cores | 24 |
| git | 2.53.0 |
| repository commit | `bae14ee24b58fb35f8a029f0b5528342d7a99599` |
| SUNDIALS submodule commit | `c1be0885d1009357993e81b65e484c3a55743f79` |

---

## 2. One-time setup

Install the Rust toolchain from <https://rustup.rs> if `cargo` is missing. On
Debian or Ubuntu, `python3-venv` is needed because the system Python refuses
system-wide `pip` installs:

```bash
sudo apt install git python3 python3-venv build-essential
```

Then get the repository, create the environment, and build the reference
integrator the notebook drives:

```bash
git clone --recurse-submodules https://github.com/once-ere/Pre-Universe_opus-fable.git
cd Pre-Universe_opus-fable
mkdir -p logs
uv venv .venv && uv pip install --python .venv/bin/python -r requirements.txt
# or, if you prefer the standard library: python3 -m venv .venv && .venv/bin/python -m pip install -r requirements.txt
cd fable_cosmo_rs && cargo build --release 2>&1 | tee ../logs/fable_cosmo_build.log && cd ..
echo "build exit status: ${PIPESTATUS[0]}"
```

If the repository was cloned without `--recurse-submodules`, repair it first:

```bash
cd Pre-Universe_opus-fable
git submodule update --init --recursive
```

The Rust build is **required, not optional**. The notebook's helper raises
`FileNotFoundError` if the binary is absent rather than quietly substituting the
SciPy answer, so a stale or missing build can never masquerade as a verified
reference result.

---

## 3. Rebuilding and executing the notebook

The notebook is generated from a deterministic builder, so re-running the
builder always produces the identical file and the repository never churns:

```bash
cd Pre-Universe_opus-fable
.venv/bin/python scripts/build_fable_notebook.py 2>&1 | tee logs/build_fable_notebook.log
echo "exit status: ${PIPESTATUS[0]}"
```

Expected output:

```
cells written : 21 {'markdown': 11, 'code': 10}
notebook      : /path/to/Pre-Universe_opus-fable/notebooks/fable_spinor_dark_energy.ipynb
```

Execute it headlessly, embedding the real outputs in place:

```bash
cd Pre-Universe_opus-fable
.venv/bin/python -m jupyter nbconvert --to notebook --execute --inplace \
  notebooks/fable_spinor_dark_energy.ipynb \
  --ExecutePreprocessor.timeout=900 \
  --ExecutePreprocessor.record_timing=False \
  --KernelManager.transport_encryption=auto 2>&1 | tee logs/execute_fable_notebook.log
echo "exit status: ${PIPESTATUS[0]}"
```

Expected output, with no warnings:

```
[NbConvertApp] Converting notebook notebooks/fable_spinor_dark_energy.ipynb to notebook
[NbConvertApp] Writing 437618 bytes to notebooks/fable_spinor_dark_energy.ipynb
```

The `--KernelManager.transport_encryption=auto` option makes `jupyter_client`
provision CurveZMQ keys for the kernel it launches, so the five loopback ZeroMQ
sockets are encrypted. Without it, `ipykernel` 7 prints
`Kernel is running over TCP without encryption` at every launch; with it the
execution log is empty of warnings. The option changes only the transport and
has no effect on any cell output.

Then strip the one piece of interpreter-specific metadata that execution
stamps into the file, so the notebook is byte-identical whichever CPython ran it:

```bash
cd Pre-Universe_opus-fable
.venv/bin/python scripts/normalize_notebooks.py 2>&1 | tee logs/normalize_notebooks.log
```

After normalization the file is 437605 bytes; that is the size the integrity
manifest pins. The two figures are embedded at 100 dpi and each carries a
textual `alt` description, so the exported HTML contains no placeholder
"No description has been provided for this image" text.

Runtime is under a minute. The notebook contains `assert` statements at every
checkpoint, so execution **fails loudly** if any tolerance is exceeded; a
successful conversion is itself the pass certificate.

---

## 4. Checking the executed notebook

This command reports the cell counts, the embedded figures, and the number of
errors:

```bash
cd Pre-Universe_opus-fable
.venv/bin/python -c "
import json
from pathlib import Path
nb = json.loads(Path('notebooks/fable_spinor_dark_energy.ipynb').read_text())
code = [c for c in nb['cells'] if c['cell_type'] == 'code']
errors = [o for c in nb['cells'] for o in c.get('outputs', []) if o.get('output_type') == 'error']
images = [o for c in nb['cells'] for o in c.get('outputs', []) if 'image/png' in o.get('data', {})]
unrun = [c for c in code if c.get('execution_count') is None]
print('cells           :', len(nb['cells']))
print('code cells      :', len(code))
print('executed        :', len(code) - len(unrun))
print('embedded figures:', len(images))
print('errors          :', len(errors))
assert not errors and not unrun
print('OK')
"
```

Expected output:

```
cells           : 21
code cells      : 10
executed        : 10
embedded figures: 2
errors          : 0
OK
```

---

## 5. What the notebook establishes

| section | claim | how it is checked |
|---|---|---|
| 3 | the eight $\mathrm{Cl}(4,4)$ generators are real integer matrices obeying the Clifford relation, with signature $(4,4)$ | residual printed and asserted zero |
| 3 | $C$ symmetric, every $C\gamma^a$ antisymmetric | printed as `True` |
| 3 | scalar and pseudoscalar bilinears survive, the vector current vanishes identically | `max |vector current| : 0.0` |
| 4 | commutant of $\mathrm{Pin}(4,4)$ has dimension **1** — absolutely irreducible over $\mathbb{R}$ | asserted |
| 4 | commutant of $\mathrm{Spin}(4,4)$ has dimension **2** — reducible, $16=8\oplus8$ | asserted |
| 4 | the 256 Clifford words span the full real $16\times16$ algebra; the 128 even words span half | rank printed |
| 5 | $n-1$ equals the benchmark $-0.764$ with zero difference | printed |
| 6 | the pure-Rust SUNDIALS CVODE reference runs and passes its own gates | subprocess output embedded |
| 7 | SciPy Radau agrees with CVODE on all 19 physical columns, and both agree with the closed form | asserted below $10^{-10}$ and $10^{-9}$ |
| 8 | the sector runs from $w\approx0$ (dark matter) to $w\to n-1$ (dark energy), monotonically, never crossing $-1$ | printed and asserted |
| 9 | the bridge field moves $w$ at fixed potential index, and the effect vanishes exactly when $\xi=0$ | spreads printed |

### The numbers it printed on the recorded run

```
commutant dimension, Pin(4,4) : 1  -> irreducible
commutant dimension, Spin(4,4): 2  -> reducible, 8 + 8

rank of the 256 Clifford words: 256 (the full real 16x16 matrix algebra)
rank of the 128 even words    : 128
```

```
SciPy Radau vs pure-Rust CVODE, largest |a-b|/(1+|b|) over 19 columns:
  bilinear            : 2.858131e-11
  bridge              : 3.612555e-11
  bridge_h            : 5.434065e-11
  deceleration        : 4.167814e-12
  density_dust        : 2.836593e-11
  density_fable       : 2.779980e-11
  density_potential   : 5.071802e-12
  dilution            : 3.093799e-12
  e_folds             : 0.000000e+00
  hubble_over_h0      : 1.053052e-11
  omega_fable         : 2.376723e-12
  pressure_dust       : 1.676068e-11
  pressure_fable      : 1.532891e-11
  pressure_potential  : 4.542126e-12
  redshift            : 0.000000e+00
  scale_factor        : 0.000000e+00
  w_dust              : 3.949653e-12
  w_fable             : 3.541414e-12
  w_potential         : 5.396753e-13

both integrators vs the closed form : 3.233014e-11
covariant conservation residual     : 1.194042e-15

All three agree.
```

```
w of the whole sector at a = 6.144e-06 : -0.000000   (dust is 0)
w of the whole sector today               : -0.558316
w of the whole sector at a = 4.034e+02 : -0.764000   (target -0.764)

monotonic, and never crosses -1 (no phantom): True
```

```
spread of w_potential, xi = 0.00 : 0.0
spread of w_potential, xi = 0.25 : 0.012249074407234994
spread of w_dust,      xi = 0.00 : 0.0
spread of w_dust,      xi = 0.25 : 0.05190285765777514
```

The last block is the control for the second mechanism: with the torsion
coupling off, both equations of state are exactly constant; with it on, both
move, at a fixed potential index.

---

## 6. Running the underlying tests

`pytest.ini` puts `src/` on the import path, so no `PYTHONPATH` is needed:

```bash
cd Pre-Universe_opus-fable
.venv/bin/python -m pytest tests/test_fable_spinor.py -q 2>&1 | tee logs/pytest_fable_spinor.log
echo "exit status: ${PIPESTATUS[0]}"
```

Expected output when the Rust reference integrator has been built:

```
...................................                                      [100%]
35 passed in 0.91s
```

One of those tests runs the Rust binary and compares it with SciPy on all 19
physical columns. If the binary has not been built yet, that single test is
skipped rather than failed, and the output reads `34 passed, 1 skipped`; the
skip message names the build command. Nothing is silently substituted.

To run every test in the repository, including the earlier `gpt5_6` suite:

```bash
cd Pre-Universe_opus-fable
.venv/bin/python -m pytest -q 2>&1 | tee logs/pytest_all.log
echo "exit status: ${PIPESTATUS[0]}"
```

Expected: `45 passed`.

---

## 7. Viewing the notebook in a web browser

### Live, interactive

```bash
cd Pre-Universe_opus-fable
.venv/bin/python -m jupyter lab notebooks/fable_spinor_dark_energy.ipynb
```

JupyterLab prints an address of the form
`http://localhost:8888/lab?token=<hex>`. Open it in a browser, then choose
**Run > Run All Cells**. Stop the server with `Ctrl-C` twice in the terminal.

### Static HTML, no kernel needed

```bash
cd Pre-Universe_opus-fable
mkdir -p build
.venv/bin/python -m jupyter nbconvert --to html \
  --output-dir build \
  --output fable_spinor_dark_energy.html \
  notebooks/fable_spinor_dark_energy.ipynb 2>&1 | tee logs/nbconvert_html.log
python3 -m http.server 8911 --directory build
```

Then open <http://127.0.0.1:8911/fable_spinor_dark_energy.html>. The page is
self-contained: prose, code, printed output and both figures. Stop the server
with `Ctrl-C`.

Port 8911 is used deliberately, because the related rustSolveIt project occupies
ports 8895–8907. `build/` is ignored by version control, so this step never
dirties the repository.

### The standalone figures

```bash
cd Pre-Universe_opus-fable
python3 -m http.server 8911
```

Then open <http://127.0.0.1:8911/artifacts/fable/figures/> and click through
`fable_equation_of_state.png`, `fable_dark_sector.png`,
`fable_second_mechanism.png`, `fable_expansion.png` and
`fable_cross_check.png`.

---

## 8. Honest limits of this page

- The notebook computes a **background** only. No perturbations are evolved, so
  nothing in it constrains structure formation, the CMB, or any observable that
  depends on more than the expansion history.
- $w=-0.764$ appears as an **exact limit of the model**, not as a fit to data.
  No supernova, BAO or CMB likelihood is evaluated anywhere in this repository,
  and no observational dataset is shipped with it.
- The agreement figures compare three routes to the **same equations**. They
  establish numerical correctness, not physical correctness.
- The logistic law for the bridge field is a phenomenological closure, not a
  consequence of an action.
- No claim is made that a new fundamental particle has been detected.
