# MRTOE: Finite-window thermal-activity criterion verification

This repository contains the numerical verification materials for the manuscript
*Finite-Window Effects and Nonoscillatory Positive Decisions in an Agarwood
Thermal-Activity Classifier*, prepared for the journal **Measurement**.

The experiment recorded agarwood temperature responses. The published mathematical
analysis separately evaluates whether a score threshold can be exceeded by a
strictly decreasing constructed cooling curve. The supplied self-contained
programs reproduce the constructed controls, boundary/interior comparisons,
constant-temperature case, baseline-offset identity, and quadratic curvature
condition.

## Verification

Python 3.11 or later:

```bash
python -m pip install -r requirements.txt
python verify_controls.py --output /tmp/mrtoe_controls.json
python verify_additional_checks.py > /tmp/mrtoe_additional.json
python verify_expected.py
```

The final command compares the computed output against the checked reference JSON
files with a numerical tolerance of `1e-10` and fails on disagreement.

## Principal mathematical result

A strictly decreasing quadratic constructed control gives an interior-only score
of approximately **0.155994 °C**, exceeding the tested threshold
**0.120062 °C**, without a directional reversal. Its distinct
window-normalized full-record score is approximately **−0.105551 °C**.
These quantities refer to different integration domains.

## Scope of the repository

The scripts operate on **constructed mathematical controls**. Original FLIR
records, laboratory reports, and the archived MARC implementation are retained
in the research owner's private records; they are not distributed in this
repository. The manuscript's access statement governs requests for those
materials. The two types of evidence must not be conflated.

## Repository structure

- `verify_controls.py` — mathematical controls, threshold, and boundary/interior checks
- `verify_additional_checks.py` — offset identity and quadratic curvature condition
- `verify_expected.py` — reproducibility comparison against expected JSON outputs
- `checks/` — reference numerical outputs supplied with the research supplement
- `.github/workflows/verify.yml` — automatic reproducibility run on the main branch

## Research ownership and rights

Copyright © 2026 **Meshaal Obaid Mouawad**. All rights reserved.
Publication of these files does not grant a license to the original research
software, FLIR records, or other restricted materials.

### Research reference

Original public deposit: https://doi.org/10.5281/zenodo.20683004

Journal manuscript: *Measurement*, Elsevier (submission preparation).
