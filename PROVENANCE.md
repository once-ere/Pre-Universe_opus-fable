# Provenance and verification record

## Scope

This record covers the `gpt5_6` nonlinear-spinor cosmology and the existing
`gpt-5.6_bridge` / `gpt-5.6f_rame` 4+4 geometry. It distinguishes external
scientific inputs, source code, generated artifacts, and validation evidence.

Validation date: 2026-09-17.

## External scientific input

### Supplied PDF

The user supplied this file outside the repository:

`Gmail - w equation of state parameter =w -0.764 forcing the Unite supernova data by itself to fit a flat, non-evolving dark energy model.pdf`

- Type: six-page Gmail print, not the full journal or arXiv paper.
- Producer: Skia/PDF m153.
- Size: 213843 bytes.
- SHA-256:
  `8cd08da437c57262d77b276cfb6c20ce104e77bbf3c229f6b04df03b857a34e7`.
- Extraction: `pdftotext -layout`, producing 293 lines for local review.
- Distribution: not copied into this repository.

The email's descriptions were treated as secondary commentary. Parameter
values, datasets, model labels, and statistical caveats were checked against:

- R. Camilleri et al., *Supernovae Unite: Combining Pantheon+ and DES-SN5YR*,
  arXiv:2609.05053v2,
  <https://arxiv.org/abs/2609.05053> and
  <https://arxiv.org/html/2609.05053v2>.
- Sections 7.0.3 and 7.0.4 are identified in the code records.

The encoded benchmarks are deliberately separate:

- Unite-only flat-$w$CDM:
  $(\Omega_m,w)=(0.197^{+0.056}_{-0.054},-0.764^{+0.078}_{-0.096})$.
- Unite+BAO+CMB flat-$w_0w_a$CDM:
  $(\Omega_m,w_0,w_a)=(0.305\pm0.004,
  -0.861^{+0.044}_{-0.042},-0.60^{+0.17}_{-0.19})$.

No supernova, BAO, or CMB catalog is downloaded or fitted by this repository.

### Theory references

- Cai and Wang, arXiv:0806.3890, supplies prior literature for nonlinear
  spinor dark energy and phantom-divide crossing.
- Chevallier and Polarski (2001) and Linder, astro-ph/0208512, supply the CPL
  parameterization.

### Optional solver assessment

The requested Linux repository
<https://github.com/once-ere/rustSolveIt_linux_SUNDIALS_7_8_0> was inspected at
the then-current commit
`6f58e02e53717a51375bd4bc5918edc57088d922`. It contains a vendored pure-Rust
SUNDIALS 7.8.0 implementation, its own mechanics-oriented notebook language,
and Jupyter wrappers. It was not vendored or executed here because it does not
provide a drop-in Python API for this custom FLRW state, while the smooth
three-state system is already checked against exact invariants at about
$10^{-11}$. This decision is documented in the notebook and both reports.

## Source and artifact lineage

The authoritative numerical implementation is
`src/gpt5_6_cosmology.py`. It contains:

- the 16-by-16 Clifford representation;
- the four-flavor field representative;
- published benchmark metadata;
- the reconstructed potential and its analytic derivative;
- homogeneous thermodynamics and perfect-fluid stress tensor;
- the FLRW ODE and diagnostics;
- figure generation.

The generation chain is:

1. `scripts/run_cosmology.py` imports the source module and writes
   `artifacts/gpt5_6_summary.json` plus four PNG/PDF figure pairs.
2. `scripts/build_cosmology_notebook.py` deterministically creates
   `notebooks/gpt5_6_cosmology.ipynb`.
3. `nbconvert --execute --inplace` runs all nine Python cells and stores their
   outputs in the notebook with execution timing metadata disabled.
4. `nbconvert --to html` writes the local HTML rendering under `build/`.
5. `scripts/build_documentation.sh` performs three strict `pdflatex` passes on
   `docs/gpt5_6_cosmology.tex` and
   `docs/gpt5_6_dark_sector_relationships.tex`, publishing the 11-page main
   report and 7-page dark-sector technical note beside their sources.
6. `scripts/build_gpt56_notebook.py` generates the Mathematica notebook from
   `wolfram/gpt56_bridge.wls`.

The Markdown publications and LaTeX sources repeat equations for pedagogy, but
all displayed numerical values were checked against the generated JSON and
executed notebook. The `build/`, `logs/`, and `backups/` directories are
intentionally ignored. Source notebooks, figures, JSON, Markdown, TeX, and both
compiled reports remain repository deliverables.

Generated PDFs use the fixed UTC release epoch `1789603200` (2026-09-17), and
executed notebooks omit wall-clock timing metadata. Rebuilding with the pinned
environment therefore produces byte-identical JSON, figures, notebook, and
report artifacts.

## Numerical choices

- Independent variable: $N=\ln a$.
- Interval: $N\in[-4,1]$.
- Output points: 1201.
- State: $H_0(t-t_0)$, $S$, and $\rho_\Psi/\rho_{c0}$.
- Integrator: SciPy `solve_ivp`, method `DOP853`.
- Relative tolerance: $10^{-11}$.
- Absolute tolerance: $10^{-13}$.
- Present normalization: $S_0=1$.
- Demonstration radiation density: $\Omega_{r0}=9\times10^{-5}$.
- Demonstration Hubble constant: $70\,\mathrm{km\,s^{-1}\,Mpc^{-1}}$.

The numerical root locations use linear interpolation between adjacent output
points. The phantom crossing also has an exact CPL expression, providing an
independent check.

## Verified numerical results

- $w(0.1)=-1.401$, $w(1)=-0.861$, and $w(2)=-0.261$.
- Phantom crossing: $a=0.7683333333$, $z=0.3015184382$.
- Independent NEC root: $a=0.7683333357$, $z=0.3015184342$.
- Acceleration begins: $a=0.5719453847$, $z=0.7484186896$.
- Formal future acceleration end: $a=1.8012033794$,
  $z=-0.4448156097$.
- Maximum relative bilinear error: $4.5232\times10^{-11}$.
- Maximum relative density error: $2.2269\times10^{-11}$.
- Maximum relative Friedmann error: $1.7988\times10^{-12}$.
- The Unite-only constant-$w$ comparison has no $w=-1$ crossing.

The future acceleration endpoint is explicitly an extrapolation of CPL beyond
the fitted domain.

## Validation environment

- OS: Linux 7.0.0-31-generic, x86-64 GNU/Linux.
- Python: 3.14.4.
- NumPy: 2.5.3.
- SciPy: 1.18.1.
- Matplotlib: 3.11.2.
- nbformat: 5.11.1.
- pytest: 9.1.1.
- WolframScript: 1.14.0 for Linux x86-64.
- pdfTeX: 3.141592653-2.6-1.40.28, TeX Live 2025/Debian.
- Git: 2.53.0.

Package versions are pinned in `requirements.txt` where they are direct project
dependencies. The local virtual environment is excluded from version control.

## Verification command and acceptance criteria

Run:

```bash
bash scripts/verify_all.sh
```

The command performs these gates:

1. Full Python suite: 10 tests pass.
2. Numerical JSON and all eight figure files regenerate.
3. The 24-cell notebook executes with nine code-cell outputs and zero error
   outputs.
4. All three numerical invariants remain below $2\times10^{-9}$.
5. The HTML export contains descriptive alt text for all four figures.
6. Python compilation and `git diff --check` pass.
7. The 11-page LaTeX report and 7-page dark-sector technical note build with no
   warning, error, undefined-reference, overfull-box, or underfull-box
   diagnostic in either final log.
8. Wolfram source and generated notebook each report 36 successes and zero
   failures; the notebook runner reports zero message-producing cells.
9. Every required source and generated deliverable is present and nonempty.

The Jupyter kernel may print a local transport warning about an unencrypted TCP
connection during command-line execution. No remote kernel is used, and this
does not appear in the notebook output or HTML artifact.

## Integrity manifest

`artifacts/SHA256SUMS` records every file in the release commit except the
manifest itself. Verify the delivered bytes before regenerating artifacts with:

```bash
sha256sum -c artifacts/SHA256SUMS
```

Two consecutive complete artifact rebuilds were also compared and produced
identical SHA-256 hashes for the JSON summary, all PNG/PDF figures, the executed
notebook, and the compiled report.

## Scientific exclusions

The following statements are not established by this repository:

- that `gpt5_6` is a discovered particle or fundamental field;
- that the reconstructed potential is unique;
- that the phantom regime is perturbatively or quantum mechanically stable;
- that the model fits raw supernova, BAO, CMB, or structure-growth data;
- that the future CPL extrapolation is a prediction;
- that the same fitted component is simultaneously all dark matter and dark
  energy.
