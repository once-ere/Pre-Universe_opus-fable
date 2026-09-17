# Provenance — the reference integrator (`fable_cosmo_rs`)

**Effort.** This page records how the fableSpinor background is integrated by a
pure-Rust program driving the vendored pure-Rust SUNDIALS 7.8.0 CVODE solver,
how its results are cross-checked by two completely independent routes, and how
to reproduce every number it prints. Everything needed to repeat the work is on
this page. No other file needs to be consulted, and no command below is
abbreviated.

---

## 1. Environment actually used

| item | value |
|---|---|
| OS | Ubuntu 26.04.1 LTS, x86-64 |
| kernel | Linux 7.0.0-31-generic |
| C library | GNU libc 2.43 |
| shell | bash |
| rustc / cargo | 1.96.1 (31fca3adb 2026-06-26) |
| Python | 3.14.4 |
| numpy / scipy / matplotlib | 2.5.3 / 1.18.1 / 3.11.2 |
| cores | 24 |
| git | 2.53.0 |
| repository commit | `bae14ee24b58fb35f8a029f0b5528342d7a99599` |
| SUNDIALS submodule commit | `c1be0885d1009357993e81b65e484c3a55743f79` |

---

## 2. One-time setup

Install the Rust toolchain from <https://rustup.rs> if `cargo` is missing (one
command, no administrator rights). Then:

```bash
git clone --recurse-submodules https://github.com/once-ere/Pre-Universe-GPT5_6_Sol.git
cd Pre-Universe-GPT5_6_Sol
mkdir -p logs
cargo --version
python3 --version
```

If the repository was cloned **without** `--recurse-submodules`, the engine is
missing and nothing will build. Repair it with:

```bash
cd Pre-Universe-GPT5_6_Sol
git submodule update --init --recursive
```

Confirm the engine is present and pinned to the exact reference commit:

```bash
cd Pre-Universe-GPT5_6_Sol
git -C vendor/sundials_rs rev-parse HEAD
ls vendor/sundials_rs/crates | sort
```

Expected output:

```
c1be0885d1009357993e81b65e484c3a55743f79
arkode_rs
cvode_rs
cvodes_rs
ida_rs
idas_rs
kinsol_rs
sundials_core
```

Install the Python packages used by the cross-check and the figures:

```bash
cd Pre-Universe-GPT5_6_Sol
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

---

## 3. Build and run

```bash
cd Pre-Universe-GPT5_6_Sol/fable_cosmo_rs
cargo build --release 2>&1 | tee ../logs/fable_cosmo_build.log
echo "exit status: ${PIPESTATUS[0]}"
```

Zero errors and **zero warnings** are expected; the crate declares
`#![forbid(unsafe_code)]` and `#![deny(warnings)]`, so a warning is a build
failure.

```bash
cd Pre-Universe-GPT5_6_Sol
./fable_cosmo_rs/target/release/fable_cosmo_rs --out artifacts/fable 2>&1 | tee logs/fable_cosmo_run.log
echo "exit status: ${PIPESTATUS[0]}"
```

Exit status `0` means every gated invariant held. The program returns a failure
status and prints the offending residual otherwise; it never prints a result it
has not checked.

Options, all optional:

```bash
./fable_cosmo_rs/target/release/fable_cosmo_rs --help
./fable_cosmo_rs/target/release/fable_cosmo_rs --xi 0.0 --out artifacts/fable-xi0
./fable_cosmo_rs/target/release/fable_cosmo_rs --n-start -12 --n-end 6 --samples 4001 --out artifacts/fable-wide
```

---

## 4. What it integrates

Independent variable is the number of e-folds $N=\ln a$; the state is
$[\ln S,\;s]$, where $S$ is the fableSpinor scalar bilinear and $s$ is the
gpt-5.6_bridge field:

$$\frac{\mathrm{d}\ln S}{\mathrm{d}N}=-\nu(s),
\qquad
\frac{\mathrm{d}s}{\mathrm{d}N}=\gamma\,s(1-s),$$

$$\nu(s)=3-\xi\bigl(1-h(s)\bigr),\qquad h(s)=3s^{2}-2s^{3}.$$

Derived algebraically at every output point:

| quantity | expression |
|---|---|
| energy density | $\rho=M S+\lambda S^{n}$ |
| pressure | $p=\bigl(\tfrac{\nu}{3}-1\bigr)MS+\bigl(\tfrac{n\nu}{3}-1\bigr)\lambda S^{n}$ |
| equation of state | $w=p/\rho$ |
| dust-like component | $w_{\text{dust}}=\nu/3-1$ |
| dark-energy-like component | $w_{\text{potential}}=n\nu/3-1$ |
| expansion rate | $H/H_0=\sqrt{\Omega_r a^{-4}+\Omega_b a^{-3}+\rho}$ |

### Solver configuration

| setting | value |
|---|---|
| method | `CV_BDF` with Newton iteration |
| linear solver | `SUNDenseMatrix` + `SUNLinSol_Dense` |
| relative tolerance | `1e-12` |
| absolute tolerances | `1e-14` per component |
| step ceiling | `500000` |
| user data | `Option<Box<dyn Any>>`, downcast in the callback |
| right-hand side | a plain `fn` pointer matching `CVRhsFn` |

The `N_Vector` borrow guard rule is obeyed throughout: a live `RefMut` from
`N_VGetArrayPointer` is never held across a `CVode` call. Every returned flag is
checked and converted into a named error; no flag is discarded.

---

## 5. The three-way check

The program is only one of three routes to the same numbers.

1. **Closed form.** The system integrates exactly:

   $$s(N)=\frac{s_0e^{\gamma N}}{1-s_0+s_0e^{\gamma N}},\qquad
   \ln S(N)=\ln S_0-3N+\frac{\xi}{\gamma}\Bigl[\ln s+s-s^{2}\Bigr]_{s(0)}^{s(N)} .$$

2. **Pure-Rust SUNDIALS CVODE** — the reference, this program.
3. **SciPy Radau**, `rtol=1e-12`, `atol=1e-14` — an independent cross-check.

Run all three and compare:

```bash
cd Pre-Universe-GPT5_6_Sol
.venv/bin/python scripts/run_fable_cosmology.py 2>&1 | tee logs/run_fable_cosmology.log
echo "exit status: ${PIPESTATUS[0]}"
```

### Gated invariants

| gate | tolerance | meaning |
|---|---|---|
| continuity residual | $10^{-13}$ | $\mathrm{d}\rho/\mathrm{d}N+3(\rho+p)=0$ |
| closed-form residual | $10^{-9}$ | integrator versus the exact solution |
| bridge residual | $10^{-9}$ | integrated $s$ versus the logistic solution |
| cross-check, $\lvert a-b\rvert/(1+\lvert b\rvert)$ | $10^{-10}$ | SciPy versus CVODE |
| benchmark error | $10^{-12}$ | $w_{\text{potential}}=-0.764$ when $\xi=0$ |
| dust error | $10^{-12}$ | $w_{\text{dust}}=0$ when $\xi=0$ |

The mixed measure $\lvert a-b\rvert/(1+\lvert b\rvert)$ is used rather than a
bare absolute difference because $S$ spans nine decades over the integration
range; for the bounded quantities the denominator is close to one and the
measure reduces to the absolute difference.

---

## 6. The result

Full output of `logs/fable_cosmo_run.log` from the recorded run:

```
fable_cosmo_rs — fableSpinor unified dark sector
engine                        : sundials_rs 7.8.0 (pure Rust, Linux x86-64)
method                        : CVODE BDF + Newton + dense
e-fold range                  : -7.003 .. 1.098612289  (1601 samples)
CVODE steps / rhs evaluations  : 261 / 291

today (N = 0):
  bridge field s             : 0.5
  dilution exponent nu       : 2.925
  w of the whole sector      : -0.569361147743
  w of the dust-like part    : -0.025
  w of the potential part    : -0.7699
  Omega_fable                : 0.95091
  deceleration parameter     : -0.3120718135

benchmark limit (xi = 0):
  w of the potential part    : -0.764
  target                     : -0.764
  |difference|               : 0.000000e+00
  w of the dust-like part    : 0.000000e+00

gated invariants:
  max |continuity residual|  : 1.194042e-15   (limit 1.000000e-13)
  max |closed-form residual| : 3.233014e-11   (limit 1.000000e-09)
  max |bridge residual|      : 5.396339e-11   (limit 1.000000e-09)
  |benchmark - (-0.764)|     : 0.000000e+00   (limit 1.000000e-12)

SUCCESS: every gated invariant holds.
```

Cross-check section of `logs/run_fable_cosmology.log` from the recorded run:

```
SciPy Radau vs pure-Rust CVODE, largest |a-b|/(1+|b|):
  bilinear          : 2.858131e-11
  bridge            : 3.612555e-11
  e_folds           : 0.000000e+00
  hubble_over_h0    : 1.053052e-11
  w_fable           : 3.541414e-12
  w_potential       : 5.396753e-13

max |closed-form residual| : 3.233014e-11  (gate 1.0e-09)
max |continuity residual|  : 1.194042e-15  (gate 1.0e-13)
benchmark |w_pot - (-0.764)|: 0.000000e+00  (gate 1.0e-12)
benchmark |w_dust|          : 0.000000e+00  (gate 1.0e-12)

SUCCESS: every gate holds.
```

Reading the two headline rows: switching the torsion coupling off makes
$w_{\text{potential}}$ exactly $-0.764$ and $w_{\text{dust}}$ exactly $0$ — a
flat, non-evolving dark energy alongside pressureless dark matter, from one
field. Switching it on ($\xi=0.15$) moves them to $-0.7699$ and $-0.025$: the
second mechanism at work.

---

## 7. Files written

| path | content | tracked? |
|---|---|---|
| `artifacts/fable/fable_background.csv` | 1601 rows, 22 columns, 16-digit reals | no — regenerable bulk data |
| `artifacts/fable/fable_summary.json` | parameters, today's values, invariants, CVODE statistics | yes |
| `artifacts/fable/fable_published_summary.json` | the above plus all cross-check gates | yes |
| `artifacts/fable/figures/*.png`, `*.pdf` | five figures | yes |

Inspect the CSV directly:

```bash
cd Pre-Universe-GPT5_6_Sol
head -1 artifacts/fable/fable_background.csv | tr ',' '\n' | nl
wc -l artifacts/fable/fable_background.csv
```

Read the summary:

```bash
cd Pre-Universe-GPT5_6_Sol
cat artifacts/fable/fable_published_summary.json
```

---

## 8. Checking the build posture yourself

No `unsafe`, and no external crates:

```bash
cd Pre-Universe-GPT5_6_Sol
grep -rn "unsafe" fable_cosmo_rs/src/
grep -c '^\[\[package\]\]' fable_cosmo_rs/Cargo.lock
grep -A1 '^\[\[package\]\]' fable_cosmo_rs/Cargo.lock | grep '^name'
```

Expected: the only `unsafe` occurrences are the `#![forbid(unsafe_code)]`
attributes, and `Cargo.lock` lists exactly three packages — `cvode_rs`,
`fable_cosmo_rs` and `sundials_core` — all local, none from a registry.

Rebuild from scratch to confirm a clean tree builds warning-free:

```bash
cd Pre-Universe-GPT5_6_Sol/fable_cosmo_rs
cargo clean
cargo build --release 2>&1 | tee ../logs/fable_cosmo_rebuild.log
grep -c warning ../logs/fable_cosmo_rebuild.log
```

Expected: `0`.

---

## 9. Viewing the result in a web browser

Serve the repository over HTTP from its root:

```bash
cd Pre-Universe-GPT5_6_Sol
python3 -m http.server 8911
```

Then open:

- <http://127.0.0.1:8911/docs/PROVENANCE-FABLE_COSMO_RS.md> — this page
- <http://127.0.0.1:8911/artifacts/fable/figures/fable_equation_of_state.png>
- <http://127.0.0.1:8911/artifacts/fable/figures/fable_dark_sector.png>
- <http://127.0.0.1:8911/artifacts/fable/figures/fable_second_mechanism.png>
- <http://127.0.0.1:8911/artifacts/fable/figures/fable_expansion.png>
- <http://127.0.0.1:8911/artifacts/fable/figures/fable_cross_check.png>
- <http://127.0.0.1:8911/artifacts/fable/fable_published_summary.json>

Stop the server with `Ctrl-C`. Port 8911 is used deliberately: the related
rustSolveIt project occupies 8895–8907.

For the same physics as a live, executed notebook with embedded plots:

```bash
cd Pre-Universe-GPT5_6_Sol
.venv/bin/python -m jupyter lab notebooks/fable_spinor_dark_energy.ipynb
```

JupyterLab prints a `http://localhost:8888/lab?token=...` address; open it and
choose **Run > Run All Cells**.

---

## 10. Honest limits of this page

- The integration is of a **background** only: one homogeneous field in a
  spatially flat universe. No perturbations, so nothing here bears on structure
  formation or the CMB.
- $w=-0.764$ is an **exact limit of the model**, not a fit. No supernova, BAO or
  CMB likelihood is evaluated anywhere in this repository, and no dataset is
  shipped with it.
- The logistic law for the bridge field is a **phenomenological closure**. It is
  smooth, stays in $(0,1)$ and integrates in closed form, which is why it was
  chosen; it is not derived from an action.
- The agreement figures in section 6 measure **numerical** agreement between
  three routes to the same equations. They say nothing about whether those
  equations describe the real universe.
