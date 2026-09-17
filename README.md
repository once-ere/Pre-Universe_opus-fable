# gpt5_6 spinor cosmology and 4+4 bridge

This repository contains two tested constructions:

- `gpt5_6`: a four-flavor, 16-component nonlinear Dirac condensate whose
  homogeneous FLRW background reconstructs a chosen time-dependent dark-energy
  equation of state.
- `gpt-5.6_bridge` and `gpt-5.6f_rame`: an exact 4+4
  canonical-to-Weitzenbock connection and frame construction in Wolfram
  Language.

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
