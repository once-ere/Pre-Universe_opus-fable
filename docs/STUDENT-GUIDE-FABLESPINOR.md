# fableSpinor — complete instructions for a student who has never done any of this

**What this is.** The full, self-contained guide to setting up and solving the
`fableSpinor` numerical cosmology from nothing: from a bare Linux machine to
rendered graphs in your web browser. Everything you need is on this page. You
will not be sent to another document for any command.

**What you will have when you finish.** A computer that has, on its own:

1. built an eight-dimensional real Clifford algebra and *proved by computation*
   that a certain 16-component field cannot be split into smaller pieces;
2. integrated a cosmology **three independent ways** — by hand in closed form,
   by a pure-Rust SUNDIALS solver, and by SciPy — and confirmed all three agree;
3. drawn five figures showing a single field behaving like **dark matter** early
   in cosmic history and like **dark energy** now;
4. reproduced the number $w=-0.764$ exactly, as a limit of the model.

---

## Words you need, each explained once

| word | meaning |
|---|---|
| **spinor** | a kind of field that needs a $720^\circ$ turn to come back to itself, unlike a vector which needs $360^\circ$. Electrons are described by one. |
| **real** | built only from ordinary numbers, no $\sqrt{-1}$ anywhere. Here every matrix entry is an *integer*. |
| **Clifford algebra** | the set of matrices $\gamma^a$ obeying $\gamma^a\gamma^b+\gamma^b\gamma^a = 2\eta^{ab}$. They are the machinery that makes spinors possible. |
| **signature $(4,4)$** | four directions that measure like space and four that measure like time. Written $\eta=\mathrm{diag}(1,1,1,1,-1,-1,-1,-1)$. |
| **irreducible** | cannot be broken into independent smaller pieces. If a field is irreducible, you are not secretly studying several smaller fields glued together. |
| **commutant** | all the matrices that commute with everything in a given family. If the only ones are multiples of the identity, the family is irreducible. That is *Schur's lemma*, and it is the whole proof. |
| **bilinear** | a number built from two copies of the field, here $S=\Psi^{\mathsf T}C\Psi$. |
| **scale factor $a$** | how stretched the universe is. $a=1$ is today; $a=0.5$ is when everything was half as far apart. |
| **redshift $z$** | another way to say the same thing: $z = 1/a - 1$. Today $z=0$; the early universe is large $z$. |
| **e-fold $N$** | $N=\ln a$. The natural clock of an expanding universe, and the variable we integrate in. |
| **equation of state $w$** | pressure divided by energy density. The one number that says how a substance gravitates: $w=0$ is dust (dark matter), $w=-1$ is a cosmological constant, $w\approx-0.76$ is dark-energy-like. |
| **ODE** | ordinary differential equation: a rule saying how fast something changes. |
| **solver / integrator** | a program that follows such a rule forward in time. |
| **BDF** | backward differentiation formulas, the solver method used here; built for problems that are numerically stiff. |
| **torsion** | a twisting of geometry that ordinary Einstein gravity sets to zero but which is allowed in general. |
| **quintessence** | a slowly-changing field proposed as dark energy, as an alternative to a fixed cosmological constant. |

---

## Part A — One-time setup

You need a Linux machine. This guide was verified on Ubuntu 26.04.1 LTS,
x86-64, GNU libc 2.43, with 24 cores; any Linux with the same tools works.

### A1. Check what you already have

Open a terminal and run each line:

```bash
git --version
python3 --version
cargo --version
cc --version
```

You want `git` 2.x, Python 3.12 or newer, `cargo` 1.80 or newer, and any C
compiler (`cargo` uses it only as the linker).

### A2. Install anything that is missing

On Ubuntu or Debian:

```bash
sudo apt update
sudo apt install git python3 python3-venv build-essential
```

If `cargo` was missing, install Rust — one command, no administrator rights
needed:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Accept the default answer, then either open a new terminal or run:

```bash
source "$HOME/.cargo/env"
cargo --version
```

### A3. Get the project

The `--recurse-submodules` flag is **not optional**: it fetches the pure-Rust
SUNDIALS solver that does the integration.

```bash
git clone --recurse-submodules https://github.com/once-ere/Pre-Universe-GPT5_6_Sol.git
cd Pre-Universe-GPT5_6_Sol
mkdir -p logs
```

If you already cloned it the wrong way, repair it in place:

```bash
cd Pre-Universe-GPT5_6_Sol
git submodule update --init --recursive
```

Check the solver arrived and is pinned to the exact tested version:

```bash
git -C vendor/sundials_rs rev-parse HEAD
ls vendor/sundials_rs/crates | sort
```

You should see:

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

### A4. Make a private Python environment

Modern Debian and Ubuntu refuse to install Python packages system-wide. A
virtual environment is a private folder of packages that cannot break your
system:

```bash
cd Pre-Universe-GPT5_6_Sol
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt
```

Check it worked:

```bash
cd Pre-Universe-GPT5_6_Sol
.venv/bin/python -c "import numpy, scipy, matplotlib; print(numpy.__version__, scipy.__version__, matplotlib.__version__)"
```

Something like `2.5.3 1.18.1 3.11.2` means you are ready.

### A5. Build the solver program

This compiles the SUNDIALS engine and the cosmology program. The first build
takes a few minutes; later builds take seconds.

```bash
cd Pre-Universe-GPT5_6_Sol/fable_cosmo_rs
cargo build --release 2>&1 | tee ../logs/fable_cosmo_build.log
echo "exit status: ${PIPESTATUS[0]}"
cd ..
```

You want exit status `0` and **no warnings at all**. The program is compiled
with `#![deny(warnings)]`, so a warning would stop the build.

---

## Part B — The physics, in plain language

You may skip to Part C and come back; nothing below is needed to run anything.

### B1. What the field is

`fableSpinor`, written $\Psi$, has **16 real components**. It lives in a world
with four space-like and four time-like directions.

The important and easily-mistaken point: those 16 components are **not** four
copies of the familiar 4-component electron spinor. If they were, you would be
studying four smaller fields wearing a trench coat. They are genuinely one
object, and the notebook proves it by computing a commutant.

The proof has two halves, and both matter:

- Under the **rotations** of this 8-dimensional world — the connected group
  $\mathrm{Spin}(4,4)$ — the commutant has dimension **2**. That means the 16
  *does* split, into $8+8$.
- Under the **full** symmetry group including reflections — the universal cover
  of the real $O(4,4)$, called $\mathrm{Pin}(4,4)$ — the commutant has dimension
  **1**. That means it does **not** split.

A reflection swaps the two halves of 8, welding them into one irreducible 16.
So the answer depends entirely on whether you include reflections, and the
question asks for the full real $O(4,4)$, so we do.

### B2. Why only one number survives

Build a number from two copies of the field, $S=\Psi^{\mathsf T}C\Psi$, where
$C$ is a fixed matrix. Because $C$ is **symmetric**, this number is generally
not zero. Because every $C\gamma^a$ is **antisymmetric**, the corresponding
"current" $\Psi^{\mathsf T}C\gamma^a\Psi$ is exactly zero for every field
configuration — you can check this yourself in Part C.

So the entire cosmology is controlled by a single number $S$.

### B3. The one rule $S$ obeys

Put the field into a uniformly expanding universe and the field equation
collapses to one statement:

$$\frac{\mathrm{d}S}{\mathrm{d}N} = -\nu\,S ,
\qquad
\nu = 3-\xi\bigl(1-h\bigr).$$

Read it as: *as the universe expands by one e-fold, $S$ shrinks by a factor
$e^{-\nu}$.* If $\nu = 3$, that is exactly how ordinary matter dilutes — spread
the same stuff over a volume that grew by $a^3$.

Here $h$ runs from 0 to 1 and says which geometry you are using: $h=0$ is
ordinary Einstein geometry, $h=1$ is the twisted "teleparallel" geometry where
all the content sits in torsion. $\xi$ says how strongly torsion couples. When
$\xi=0$, or when $h=1$, you get $\nu=3$ exactly.

### B4. Where dark matter and dark energy come from

The energy density has two pieces:

$$\rho = \underbrace{M\,S}_{\text{dilutes fast}} \;+\; \underbrace{\lambda\,S^{n}}_{\text{dilutes slowly}} ,
\qquad 0<n<1 .$$

Because $n$ is less than 1, raising the small number $S$ to the power $n$ makes
it shrink much more slowly. So:

- **Early**, when $S$ is huge, the first piece dominates. Its equation of state
  is $w=\nu/3-1$, which is exactly $0$ when $\nu=3$. That is **dust — dark
  matter**.
- **Late**, when $S$ is tiny, the second piece dominates. Its equation of state
  is $w=n\nu/3-1$, which is exactly $n-1$ when $\nu=3$. That is **dark energy**.

One field. Both behaviours. Nothing added by hand.

And the headline: choose $n=0.236$ and the late-time equation of state is
exactly $n-1=-0.764$ — a flat, **non-evolving** dark energy.

### B5. The second mechanism

Ordinary quintessence gives one way to make $w$ change: change the potential,
that is, change $n$. This model has a second, independent way.

The bridge field $s$ enters through $\nu$. Turning on $\xi$ shifts **both**
equations of state while leaving $n$ untouched, and it gives the dark-matter-like
piece a small **negative pressure** — something no pure quintessence can do.
Turn $\xi$ off and the effect vanishes exactly. That is the control experiment,
and you will run it in Part C.

---

## Part C — Doing it

### C1. Run the reference solver

```bash
cd Pre-Universe-GPT5_6_Sol
./fable_cosmo_rs/target/release/fable_cosmo_rs --out artifacts/fable 2>&1 | tee logs/fable_cosmo_run.log
echo "exit status: ${PIPESTATUS[0]}"
```

Read the output. The lines that matter:

```
benchmark limit (xi = 0):
  w of the potential part    : -0.764
  target                     : -0.764
  |difference|               : 0.000000e+00
  w of the dust-like part    : 0.000000e+00
```

Zero difference, not "small" difference. And:

```
SUCCESS: every gated invariant holds.
```

The program checks itself and refuses to exit successfully otherwise.

### C2. Cross-check it with a completely different solver

```bash
cd Pre-Universe-GPT5_6_Sol
.venv/bin/python scripts/run_fable_cosmology.py 2>&1 | tee logs/run_fable_cosmology.log
echo "exit status: ${PIPESTATUS[0]}"
```

This runs the Rust program again, then integrates the *same* equations with
SciPy's Radau method, then compares both against a formula solved by hand. All
three must agree. Expect differences around $10^{-11}$ and the final line:

```
SUCCESS: every gate holds.
```

Why bother? Because a solver that is wrong in an interesting way will usually be
wrong differently from a second, unrelated solver — and neither will match a
closed-form answer. Agreement of all three is strong evidence.

### C3. Prove the irreducibility yourself

```bash
cd Pre-Universe-GPT5_6_Sol
PYTHONPATH=src .venv/bin/python -c "
import numpy as np, fable_spinor as fs
spin = np.array([fs.LORENTZ[a, b] for a in range(8) for b in range(a + 1, 8)])
print('commutant of Pin(4,4)  :', fs.commutant_dimension(fs.GAMMA), '-> irreducible')
print('commutant of Spin(4,4) :', fs.commutant_dimension(spin), '-> reducible, 8+8')
print('rank of 256 words      :', fs.word_span_rank())
print('rank of 128 even words :', fs.word_span_rank(even_only=True))
"
```

Expected:

```
commutant of Pin(4,4)  : 1 -> irreducible
commutant of Spin(4,4) : 2 -> reducible, 8+8
rank of 256 words      : 256
rank of 128 even words : 128
```

### C4. Watch the current vanish

```bash
cd Pre-Universe-GPT5_6_Sol
PYTHONPATH=src .venv/bin/python -c "
import numpy as np, fable_spinor as fs
rng = np.random.default_rng(20260917)
for trial in range(3):
    psi = rng.standard_normal(16)
    print(f'trial {trial}:  S = {fs.fable_bilinear(psi): .6f}   max|current| = {abs(fs.fable_current(psi)).max():.1e}')
"
```

Expected — $S$ is different each time, the current is exactly zero every time:

```
trial 0:  S =  5.756479   max|current| = 0.0e+00
trial 1:  S = -2.412645   max|current| = 0.0e+00
trial 2:  S = -7.052059   max|current| = 0.0e+00
```

### C5. Run the control experiment for the second mechanism

```bash
cd Pre-Universe-GPT5_6_Sol
PYTHONPATH=src .venv/bin/python -c "
import numpy as np, fable_spinor as fs
base = vars(fs.FableParameters())
for xi in (0.0, 0.10, 0.25):
    sol = fs.solve_background(fs.FableParameters(**{**base, 'xi': xi}))
    print(f'xi = {xi:4.2f}   spread of w_potential = {np.ptp(sol.w_potential):.6f}   spread of w_dust = {np.ptp(sol.w_dust):.6f}')
"
```

Expected — exactly zero spread with the coupling off, nonzero with it on:

```
xi = 0.00   spread of w_potential = 0.000000   spread of w_dust = 0.000000
xi = 0.10   spread of w_potential = 0.004900   spread of w_dust = 0.020761
xi = 0.25   spread of w_potential = 0.012249   spread of w_dust = 0.051903
```

### C6. Run the tests

```bash
cd Pre-Universe-GPT5_6_Sol
PYTHONPATH=src .venv/bin/python -m pytest -q 2>&1 | tee logs/pytest_all.log
echo "exit status: ${PIPESTATUS[0]}"
```

Expected: `37 passed`.

### C7. Optional — the symbolic proofs

Only if you have Mathematica or the Wolfram Engine installed and activated:

```bash
cd Pre-Universe-GPT5_6_Sol
wolframscript -code '{$Version, $LicenseType}'
wolframscript -file wolfram/fable_spinor.wls 2>&1 | tee logs/fable_spinor.log
wolframscript -file wolfram/gpt56_bridge_dynamic.wls 2>&1 | tee logs/gpt56_bridge_dynamic.log
wolframscript -file wolfram/gpt56_bridge.wls 2>&1 | tee logs/gpt56_bridge.log
```

Expect 30, 17 and 36 assertions respectively, all passing, all closed
symbolically rather than numerically. If `wolframscript` is missing or
unlicensed, skip this step: nothing in Parts A to C depends on it.

---

## Part D — Seeing it in a web browser

### D1. The notebook, live

```bash
cd Pre-Universe-GPT5_6_Sol
.venv/bin/python -m jupyter lab notebooks/fable_spinor_dark_energy.ipynb
```

JupyterLab prints an address like `http://localhost:8888/lab?token=abc123...`.
Open it in your browser. Then choose **Run > Run All Cells** from the menu and
watch every number above be recomputed in front of you, with plots.

Stop the server by pressing `Ctrl-C` twice in the terminal.

### D2. The notebook, as a static page

If you just want to read it without running a kernel:

```bash
cd Pre-Universe-GPT5_6_Sol
mkdir -p build
.venv/bin/python -m jupyter nbconvert --to html \
  --output-dir build \
  --output fable_spinor_dark_energy.html \
  notebooks/fable_spinor_dark_energy.ipynb 2>&1 | tee logs/nbconvert_html.log
python3 -m http.server 8911 --directory build
```

Open <http://127.0.0.1:8911/fable_spinor_dark_energy.html>. Stop with `Ctrl-C`.

### D3. The figures

```bash
cd Pre-Universe-GPT5_6_Sol
python3 -m http.server 8911
```

Open <http://127.0.0.1:8911/artifacts/fable/figures/> and click through:

| file | what to look at |
|---|---|
| `fable_equation_of_state.png` | the solid curve sliding from about $0$ at high redshift down toward $-0.764$ |
| `fable_dark_sector.png` | two straight lines on a log-log plot crossing near today: the handover from dark matter to dark energy |
| `fable_second_mechanism.png` | four curves that collapse onto one flat line when $\xi=0$ |
| `fable_expansion.png` | the deceleration parameter crossing zero — the moment expansion started speeding up |
| `fable_cross_check.png` | all residuals far below the dashed gate line |

Stop with `Ctrl-C`. Port 8911 is chosen deliberately, since the related
rustSolveIt project uses 8895–8907.

---

## Part E — If something goes wrong

| symptom | cause | fix |
|---|---|---|
| `cargo: command not found` | Rust not installed or not on the path | `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \| sh` then `source "$HOME/.cargo/env"` |
| `failed to read .../vendor/sundials_rs/crates/cvode_rs/Cargo.toml` | cloned without submodules | `git submodule update --init --recursive` |
| `FileNotFoundError: The reference binary ... does not exist` | the Rust program was never built | `cd fable_cosmo_rs && cargo build --release && cd ..` |
| `error: externally-managed-environment` from pip | Debian/Ubuntu blocking system-wide installs | use the virtual environment from step A4, and call `.venv/bin/python` not `python3` |
| `ModuleNotFoundError: No module named 'fable_spinor'` | `PYTHONPATH` not set | prefix the command with `PYTHONPATH=src` |
| `Address already in use` on port 8911 | a previous server is still running | `python3 -m http.server 8912` and use 8912 in the address |
| `wolframscript` prints a licensing error | Wolfram not activated | skip Part C7; nothing else depends on it |
| JupyterLab does not open a browser | headless machine or no default browser | copy the printed `http://localhost:8888/lab?token=...` address by hand |

---

## Part F — What you may and may not conclude

**You may say:** this is a self-consistent background model in which one real,
irreducible 16-component field reproduces both a pressureless early component
and a late-time accelerating component; in a clean limit its dark-energy piece
has a flat, non-evolving $w=-0.764$ exactly; and it contains a second,
geometric mechanism for a time-varying $w$ that is independent of the potential.

**You may not say:** that this has been fitted to data, that it is preferred by
observation, or that a new particle has been found. None of those is true.
Specifically:

- Only the **background** is computed. No perturbations, so nothing here bears
  on galaxy formation or the cosmic microwave background.
- No observational dataset is shipped or used anywhere in this repository. The
  number $-0.764$ is treated as a **benchmark to be reproduced as a limit**, not
  as a measurement being fitted.
- The rule $\mathrm{d}s/\mathrm{d}N=\gamma s(1-s)$ governing the bridge field is
  a **phenomenological closure**, chosen because it is smooth, stays between 0
  and 1, and can be integrated by hand. It is not derived from an action.
- The agreement between the three solvers proves the **arithmetic** is right. It
  says nothing about whether the equations describe the real universe.
