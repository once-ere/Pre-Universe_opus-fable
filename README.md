# fableSpinor, the gpt-5.6 bridge, and the gpt5_6 spinor cosmology

This repository contains three tested constructions:

- `fableSpinor`: a **real** 16-component field carrying an **irreducible**
  representation of the universal cover of the **real** $O(4,4)$ — proved by
  computing a commutant, not asserted — whose background cosmology is
  dark-matter-like early and dark-energy-like today, and whose massless
  torsion-free limit has a flat, non-evolving $w=n-1$, exactly $-0.764$ at
  $n=0.236$.
- `gpt-5.6_bridge` and `gpt-5.6f_rame`: an exact 4+4 canonical-to-Weitzenbock
  connection and frame construction in Wolfram Language, whose interpolation
  parameter is promoted to a scalar **field**, making it the only connection in
  the comparison that a spinor cannot redefine away.
- `gpt5_6`: an earlier four-flavor, 16-component nonlinear Dirac condensate whose
  homogeneous FLRW background reconstructs a chosen time-dependent dark-energy
  equation of state. It is retained for provenance; note that its
  $\mathbb{1}_4\otimes\gamma^\mu$ construction is a direct sum of four $3{+}1$
  Dirac spinors and is therefore **not** $O(4,4)$-irreducible, which is exactly
  what `fableSpinor` corrects.

Start with [the theory section below](#the-theory-behind-this-repository) or the
[complete student guide](docs/STUDENT-GUIDE-FABLESPINOR.md).

The cosmology is an effective background model, not an observational fit or a
claim that a new fundamental particle has been discovered.

## Start here

- [Complete student guide](docs/gpt5_6_cosmology.md)
- [Dark-sector relationships and conclusions](docs/gpt5_6_dark_sector_relationships.md)
- [Compiled scientific report](docs/gpt5_6_cosmology.pdf)
- [Compiled dark-sector technical note](docs/gpt5_6_dark_sector_relationships.pdf)
- [Executable cosmology notebook](notebooks/gpt5_6_cosmology.ipynb)
- [Cosmology source](src/gpt5_6_cosmology.py)
- [Numerical summary](artifacts/gpt5_6_summary.json)
- [Build and source provenance](PROVENANCE.md)

## Quick verification

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m pytest -q
.venv/bin/python scripts/run_cosmology.py
.venv/bin/python scripts/build_cosmology_notebook.py
.venv/bin/python -m jupyter nbconvert --to notebook --execute --inplace \
  notebooks/gpt5_6_cosmology.ipynb --ExecutePreprocessor.timeout=180 \
  --ExecutePreprocessor.record_timing=False
bash scripts/build_documentation.sh
wolframscript -file wolfram/gpt56_bridge.wls
wolframscript -file scripts/run_gpt56_notebook.wls
```

The expected cosmology result is 10 passing tests and three numerical invariant
errors below `2e-9`. The Wolfram source and generated notebook each report 36
passing tests, zero failures, and no messages. The documentation build publishes
the 11-page scientific report and the 7-page dark-sector technical note without
TeX diagnostics.

---

# The theory behind this repository

This section is the conceptual map. Every claim below is either proved
symbolically in Wolfram Language, verified numerically by two independent
integrators, or explicitly flagged as an assumption. Each is linked to the
self-contained provenance page that reproduces it.

## 1. The setting: a 4+4 spacetime

The geometry has eight coordinates $x^0,\dots,x^7$ and a tangent metric with
four space-like and four time-like directions,

$$\eta_{ab}=\operatorname{diag}(1,1,1,1,-1,-1,-1,-1).$$

A **frame field** (an eight-dimensional vielbein) $e_\mu{}^a$ attaches a flat
$4{+}4$ Minkowski frame to every point of the curved manifold, and links the two
metrics by

$$g_{\mu\nu}=e_\mu{}^{a}\,\eta_{ab}\,e_\nu{}^{b}.$$

The canonical frame is diagonal,

$$e_\mu{}^a=\operatorname{diag}\bigl(\tan 6Hx^0,\;q,q,q,\;1,\;p,p,p\bigr),\qquad
q=\frac{e^{-a_4(Hx^4)}}{\sin^{1/6}6Hx^0},\qquad
p=\frac{e^{+a_4(Hx^4)}}{\sin^{1/6}6Hx^0},$$

with $a_4$ a free function. The two four-blocks carry independent scale factors
constrained only by $pq=\sin^{-1/3}6Hx^0$ — the geometric root of the second
mechanism in section 5.

Solving the **vielbein postulate**

$$\partial_\mu e_\nu{}^a-\Gamma^\rho_{\mu\nu}e_\rho{}^a+\omega_\mu{}^a{}_b\,e_\nu{}^b=0$$

with the Levi-Civita $\Gamma$ yields the canonical spin connection: 37 nonzero
Christoffel symbols of 512, **24** nonzero spin-connection components of 512,
antisymmetric in the flat pair, metric compatible and torsion free. Acting on a
spinor,

$$D_\mu\psi=\partial_\mu\psi+\tfrac18\,\omega_{\mu ab}\,[\gamma^a,\gamma^b]\,\psi .$$

Because each $\gamma^a$ is block off-diagonal in the parity grading, every
commutator $[\gamma^a,\gamma^b]$ is block diagonal, so this derivative preserves
the eight-dimensional **type-one split-octonion block**.

*Reproduced by* [docs/PROVENANCE-GPT-5.6_BRIDGE-DYNAMIC.md](docs/PROVENANCE-GPT-5.6_BRIDGE-DYNAMIC.md).

## 2. gpt-5.6f_rame: a frame that mixes the blocks

Let $R$ be the normalised order-four Hadamard matrix and

$$J=\begin{pmatrix}0&R\\R^{\mathsf T}&0\end{pmatrix},\qquad
J^{\mathsf T}\eta+\eta J=0,\qquad J^{2}=\mathbb{1}_8 .$$

The first identity places $J$ in $\mathfrak{so}(4,4)$. The second is the useful
accident: because $J$ squares to the identity, the group element needs no
series,

$$\Lambda(u)=\cosh u\,\mathbb{1}+\sinh u\,J,\qquad \Lambda^{\mathsf T}\eta\,\Lambda=\eta,$$

a simultaneous equal-rapidity boost in all four orthogonal $(+,-)$ planes.
`gpt-5.6f_rame` is $e\cdot\Lambda(u)$. It reconstructs the *same* curved metric
for every $u$, reduces to the canonical frame at $u=0$, and — unlike the
diagonal canonical frame — is dense: it genuinely mixes the $+4$ and $-4$
blocks.

## 3. gpt-5.6_bridge: the parameter becomes a field

Let $K=\Gamma_{\mathrm W}-\Gamma_{\mathrm{LC}}$ be the difference between the
Weitzenböck and Levi-Civita connections. A difference of connections is a
**tensor**, so for any scalar $h$

$$\Gamma(h)=\Gamma_{\mathrm{LC}}+h\,K$$

is again a connection; with $h(s)=3s^2-2s^3$ the endpoints are exact — $h=0$ is
Levi-Civita, $h=1$ is the Weitzenböck geometry. One family therefore contains
**both** prior geometries exactly.

Nothing in that argument required $h$ to be constant. Promoting it to a scalar
**field** $h(s(x))$ keeps the product tensorial, so the object remains a
metric-compatible $\mathfrak{so}(4,4)$ connection — and it is not free. Two
exact transgression identities follow:

$$T=h\,T_{\mathrm{Weitzenb\ddot{o}ck}},\qquad
R=R_{\mathrm{LC}}+h\,D_{\mathrm{LC}}K+h^{2}K\wedge K+\underbrace{\mathrm{d}h\wedge K}_{\text{new}} .$$

The final term is identically zero for every constant parameter, so no member of
the constant-parameter family reproduces it.

### Why this is superior to a Weitzenböck connection — and the criterion that failed

The obvious claim would be that the bridge has nonzero **axial** torsion where
Weitzenböck's vanishes. **That claim is false, and the test suite proves it
false**: the bridge torsion is a scalar multiple of the Weitzenböck torsion, so
its totally antisymmetric part is also identically zero. The negative result is
kept as two named assertions so it cannot be quietly forgotten.

The criterion that *does* hold is **irremovability**. With torsion potential
$\mathcal T=\log\sin 6Hx^0$, the entire family contributes one scalar one-form
to the Dirac operator,

$$\gamma^\mu\Gamma^{\mathrm{spin}}_\mu=-\tfrac12(1-h)\,\partial_\mu\mathcal T\,\gamma^\mu ,$$

which a spinor rescaling $\psi\to e^{-F}\psi$ absorbs **if and only if** that
one-form is exact — that is, if and only if the obstruction two-form vanishes:

$$\Omega=\mathrm{d}\bigl[(1-h)\,\mathrm{d}\mathcal T\bigr]=\tfrac12\,\mathrm{d}h\wedge\mathrm{d}\mathcal T .$$

| connection | $\Omega$ | consequence |
|---|---|---|
| canonical Levi-Civita, $h=0$ | $0$ | removable by a rescaling |
| Weitzenböck, $h=1$ | $0$ | the connection vanishes outright |
| any **constant** bridge parameter | $0$ | removable by a rescaling |
| **bridge field**, $\partial h\nparallel\partial\mathcal T$ | $\neq 0$ | **irremovable** |

For $s=\tfrac12+\tfrac14\tanh x^4$ the exact witness is

$$\Omega_{x^0x^4}=\frac{9H\cot(6Hx^0)\,\operatorname{sech}^2(x^4)\bigl(\tanh^2(x^4)-4\bigr)}{32}\neq0 .$$

That, and nothing wider, is the measurable sense in which `gpt-5.6_bridge` is
superior: it is the only entry in the table that a spinor cannot be redefined to
ignore.

*Reproduced by* [docs/PROVENANCE-GPT-5.6_BRIDGE-DYNAMIC.md](docs/PROVENANCE-GPT-5.6_BRIDGE-DYNAMIC.md).

## 4. fableSpinor: a real 16-component field that does not split

`fableSpinor` is a **real** 16-component field. The generators of
$\mathrm{Cl}(4,4)$ are raising and lowering operators on the exterior algebra of
$\mathbb{R}^4$: four square to $+1$ and are symmetric, four square to $-1$ and
are antisymmetric, and **every entry is an integer**.

The representation-theoretic point is the one most easily got wrong, so it is
*computed* rather than asserted. By real Schur's lemma a representation is
absolutely irreducible over $\mathbb{R}$ exactly when its commutant is
one-dimensional:

| group | commutant dimension | verdict |
|---|---|---|
| $\mathrm{Pin}(4,4)$, the universal cover of the **real** $O(4,4)$ | **1** | absolutely irreducible |
| $\mathrm{Spin}(4,4)$, the connected subgroup | **2** | reducible, $16=8\oplus8$ |

The second row is why the *full* $O(4,4)$ is required: the chirality operator
commutes with every $S^{ab}$ but anticommutes with every $\gamma^a$, so a
reflection exchanges the two halves of eight and welds them into one irreducible
sixteen. A construction such as $\mathbb{1}_4\otimes\gamma^\mu$ — a direct sum of
four $3{+}1$ Dirac spinors — does **not** answer the question.

With the spinor metric $C=\gamma^1\gamma^2\gamma^3\gamma^4$ one has
$C^{\mathsf T}=C$ and $(C\gamma^a)^{\mathsf T}=-C\gamma^a$. For a commuting real
field that is precisely the right pairing: the scalar
$S=\Psi^{\mathsf T}C\Psi$ survives, the vector current
$\Psi^{\mathsf T}C\gamma^a\Psi$ vanishes **identically**, and the kinetic term is
a genuine first-order symplectic form rather than a total derivative.

$$\mathcal L=e\Bigl[\tfrac12\bigl(\Psi^{\mathsf T}C\gamma^\mu D_\mu\Psi-(D_\mu\Psi)^{\mathsf T}C\gamma^\mu\Psi\bigr)-MS-V(S)\Bigr],
\qquad \gamma^\mu D_\mu\Psi=\bigl(M+V'(S)\bigr)\Psi .$$

*Reproduced by* [docs/PROVENANCE-FABLESPINOR-ALGEBRA.md](docs/PROVENANCE-FABLESPINOR-ALGEBRA.md).

## 5. The cosmology: one field, two behaviours

Contracting the homogeneous field equation with $\Psi^{\mathsf T}C$ — using only
the symmetry of $C$, the antisymmetry of $C\gamma^0$, and a **generic**
sixteen-component $\Psi$ — collapses the whole system to one exact law:

$$\frac{\mathrm{d}S}{\mathrm{d}N}=-\nu(s)\,S,\qquad
\nu(s)=3-\xi\bigl(1-h(s)\bigr),\qquad N=\ln a .$$

The energy density is $\rho=MS+V(S)$, and covariant conservation then fixes the
pressure uniquely. For $V(S)=\lambda S^{n}$ with $0<n<1$:

| quantity | expression |
|---|---|
| kinetic energy density | $(M+V'(S))\,S$ |
| potential energy density | $\lambda S^{n}$ |
| energy density | $\rho=MS+\lambda S^{n}$ |
| pressure | $p=\bigl(\tfrac{\nu}{3}-1\bigr)MS+\bigl(\tfrac{n\nu}{3}-1\bigr)\lambda S^{n}$ |
| equation of state | $w=p/\rho$ |
| dust-like component | $w_{\mathrm{dust}}=\nu/3-1$ |
| dark-energy-like component | $w_{\mathrm{potential}}=n\nu/3-1$ |

### Dark matter and dark energy, from the same field

Early, $S$ is large and $MS$ dominates; since $w_{\mathrm{dust}}$ vanishes
exactly at $\nu=3$, the sector is pressureless and gravitates as **dark
matter**. Late, $S$ is small and — because $n<1$ — the term $\lambda S^{n}$ falls
far more slowly and takes over, gravitating as **dark energy**. The handover
needs no second field. The numerical solution confirms the transition is
monotonic and never crosses $w=-1$, so the model cannot counterfeit phantom
behaviour.

### Two independent mechanisms for a time-varying $w$

1. **Potential** (a generalisation of quintessence): the ratio
   $\lambda S^{\,n-1}/M$ grows as $a^{3(1-n)}$, sliding the total $w$ from $0$
   toward $n-1$.
2. **Geometric** (supplied by `gpt-5.6_bridge`): the bridge field enters $\nu$
   and shifts *both* components' $w$ at fixed $n$, giving the dark-matter-like
   component a small **negative pressure** — something no pure quintessence can
   do. Setting $\xi=0$ switches the effect off exactly, which is the control.

### The benchmark

Switching the torsion coupling off makes $w_{\mathrm{potential}}$ exactly $n-1$
— a flat, **non-evolving** equation of state. Hence $n=0.236$ reproduces
$w=-0.764$ with **zero** numerical error, while $w_{\mathrm{dust}}$ is exactly
zero. The Unite-only constant-$w$ benchmark is a clean *limit* of the model, not
a fitted output; no likelihood is evaluated anywhere in this repository.

*Reproduced by* [docs/PROVENANCE-FABLE_COSMO_RS.md](docs/PROVENANCE-FABLE_COSMO_RS.md)
and [docs/PROVENANCE-FABLESPINOR-NOTEBOOK.md](docs/PROVENANCE-FABLESPINOR-NOTEBOOK.md).

## 6. How the numbers are trusted

The background is integrated **three independent ways** — a closed-form solution
derived by hand, the pure-Rust SUNDIALS 7.8.0 CVODE reference, and SciPy Radau —
and every comparison is gated, with nonzero exit on violation.

| quantity | measured | gate |
|---|---|---|
| covariant conservation residual | $1.194042\times10^{-15}$ | $10^{-13}$ |
| integrators versus the closed form | $3.233014\times10^{-11}$ | $10^{-9}$ |
| SciPy versus CVODE, $\lvert a-b\rvert/(1+\lvert b\rvert)$ | $3.612555\times10^{-11}$ | $10^{-10}$ |
| $\lvert w_{\mathrm{potential}}-(-0.764)\rvert$ at $\xi=0$ | $0$ | $10^{-12}$ |
| $\lvert w_{\mathrm{dust}}\rvert$ at $\xi=0$ | $0$ | $10^{-12}$ |

The last two rows are **exact**, not merely small.

---

# Documentation index

Every provenance page is self-contained: it repeats the full environment, every
command, and the expected output, and never redirects you to another document.

| document | what it reproduces |
|---|---|
| [docs/STUDENT-GUIDE-FABLESPINOR.md](docs/STUDENT-GUIDE-FABLESPINOR.md) | bare Linux machine to rendered graphs, for a reader new to all of it |
| [docs/PROVENANCE-FABLESPINOR-ALGEBRA.md](docs/PROVENANCE-FABLESPINOR-ALGEBRA.md) | the real Cl(4,4) algebra, the irreducibility proof, Lagrangian, stress tensor, $\rho$, $p$, $w$ |
| [docs/PROVENANCE-GPT-5.6_BRIDGE-DYNAMIC.md](docs/PROVENANCE-GPT-5.6_BRIDGE-DYNAMIC.md) | `gpt-5.6f_rame`, `gpt-5.6_bridge`, the transgression identities, the failed axial criterion, the obstruction |
| [docs/PROVENANCE-FABLE_COSMO_RS.md](docs/PROVENANCE-FABLE_COSMO_RS.md) | the pure-Rust SUNDIALS CVODE reference integration and its gates |
| [docs/PROVENANCE-FABLESPINOR-NOTEBOOK.md](docs/PROVENANCE-FABLESPINOR-NOTEBOOK.md) | authoring, executing, checking and displaying the notebook |
| [docs/PROVENANCE-FULL-EXECUTION.md](docs/PROVENANCE-FULL-EXECUTION.md) | one end-to-end execution of every notebook and every gate in the repository |
| [docs/fable_spinor.md](docs/fable_spinor.md) / [.tex](docs/fable_spinor.tex) / [.pdf](docs/fable_spinor.pdf) | the compiled scientific report |
| [PROVENANCE.md](PROVENANCE.md) | build and source provenance of the earlier `gpt5_6` release |

## Reproducing everything with one command

```bash
git clone --recurse-submodules https://github.com/once-ere/Pre-Universe-GPT5_6_Sol.git
cd Pre-Universe-GPT5_6_Sol
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
mkdir -p logs build
bash scripts/verify_all.sh 2>&1 | tee logs/verify_all.log
```

That runs all twelve gate stages: Python tests, numerical artifacts, both
Jupyter notebooks built and executed, structured invariant checks, hygiene,
three LaTeX reports, the Wolfram source and generated Mathematica notebook, the
required deliverables, the fableSpinor symbolic proofs, the pure-Rust SUNDIALS
reference integration, the three-way cross-check, and the fableSpinor notebook,
report and deliverables. The final line is
`All repository verification gates passed.`

## What is not claimed

- Only the **background** is computed. No perturbations, so nothing here
  constrains structure formation or the microwave background.
- **No observational dataset is shipped or used**, and no likelihood is
  evaluated. $-0.764$ is a benchmark reproduced as a limit, never data being
  fitted.
- The flow law $\mathrm{d}s/\mathrm{d}N=\gamma s(1-s)$ for the bridge field is a
  **phenomenological closure**, chosen for smoothness, boundedness and analytic
  integrability. It is not derived from an action.
- `gpt-5.6f_rame` is a pseudo-orthogonal frame gauge transformation; by itself it
  changes no observable.
- The axial-torsion superiority criterion **failed**; the claim rests on
  containment, the gradient transgression term, and irremovability.
- Agreement among three solvers establishes numerical correctness, not physical
  correctness.
- No claim is made that a new fundamental particle has been detected.
