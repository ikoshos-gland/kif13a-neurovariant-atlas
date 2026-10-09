# KIF13A K1122E local monomer protocol

This frozen exploratory protocol tests local helix and contact changes in the
AlphaFold Q9H1H9 v6 tail. It does not represent the dimer, an interacting partner,
the full protein or a clinically validated pathogenicity assay.

## Systems and controls

`protocol.json` v3 uses **ff19SB + OPC** as the main model and defines 30 tasks:
WT, K1122E, K1122R, K1104E and E1128K, each
with three independently seeded replicas at two boundaries. The primary fragment
is 1092-1146 (100 ns production); the padded sensitivity fragment is 1084-1154
(50 ns). Synthetic controls have unknown clinical effects.

ACE/NME caps inherit neighboring AF backbone coordinates. Bond lengths and Amber
template assignment must pass before production. OpenMM 8.3.1 requires its
standard hydrogen loader to run before the custom NME CH3-parent definition.
`data/interim/k1122e_md_v3/cap_validation_local_20261004.json` records successful
ff19SB parameterization of all ten variant/boundary combinations using
PDBFixer 1.12.0 and OpenMM 8.3.1. The earlier cap record is historical v2 evidence.

K1122 pLDDT is 80.31; the shared 1100-1138 core has mean 83.47 and minimum 74.94.
The padded fragment reaches 49.41 at its edge. Local confidence does not validate
oligomerization, side-chain contacts or global tail arrangement.

## Commands

Install `environment.yml` or use the TRUBA runtime assembled under
`scripts/truba/k1122e/`. With that interpreter, run from the repository root:

```text
python scripts/k1122e_md/common.py --output RUN/tasks.json
python scripts/k1122e_md/validate_caps.py --output RUN/cap_validation.json
python scripts/k1122e_md/build.py --task-index 0 --output-root RUN --platform CPU
python scripts/k1122e_md/run.py --task-dir RUN/primary_WT_r01 --platform CUDA
python scripts/k1122e_md/analyze.py --task-dir RUN/primary_WT_r01
python scripts/k1122e_md/analyze.py --run-root RUN --output-root RESULTS
```

Use separate output roots for `run.py --smoke` and `--duration-ns` pilots.
Smoke runs have 30 equilibration steps and 100 production steps with ten frames.
`analyze.py --task-dir ... --allow-smoke` produces diagnostics only. Production
uses 310 K, 1 bar, 0.15 M NaCl, ff19SB/OPC, 2 fs steps and 10 ps reporting.
OpenMM 8.3.1 already bundles the selected Amber19 files. `tip4pew` supplies
initial four-site box coordinates; `amber19/opc.xml` supplies OPC water and
compatible ion parameters. Constrained minimization corrects the initial water
geometry. Every water's charges, massless virtual site and weights are checked;
build records contain parameter-file hashes and ff19SB CMAP torsion counts.
The cubic box edge is the largest initial atom extent plus 2.4 nm, allowing free
tumbling without the short dimensions of an elongated rectangular box.
Resume with `--resume` on the same platform and exact runtime; mismatched DCD
frames, protocol hashes or durations are rejected.

## Frozen endpoints and interpretation

Primary endpoints are the core alpha-torsion fraction and geometric i,i+4
hydrogen-bond fraction. Alpha torsions use phi [-100,-30] and psi [-80,0] degrees;
this is not DSSP. Hydrogen bonds require O-N distance <=0.35 nm and N-H-O angle
>=150 degrees. Secondary endpoints describe site helix, heavy-atom contacts,
charged-group contacts and core CA RMSD. Exact definitions are in the config.

Discard the first 10 ns of production. Both boundaries use **10-50 ns** for
replica contrasts, avoiding a duration confound in the boundary comparison.
Full-duration means, contact tables and traces remain descriptive diagnostics.
Autocorrelation, complete 10 ns block
means and half-trajectory drift describe sampling within each run. Statistical
comparisons resample the three independent replica means, never frames or blocks.
Intervals are exploratory with this small replica count. Any solute self-image
contact blocks contrasts; primary endpoint direction must agree across both
boundaries before a locally consistent exploratory signal is reported.

The historical ff14SB/GBn2 static proxy is disabled in v3. ThermoMPNN remains
a separate stability predictor. Absolute potentials across mutants with
different chemistry must not be subtracted as folding ddG or binding energy.
All derived artifacts preserve input hashes in `.prov.json` sidecars.

## TRUBA deployment

The deployed code lives in `scripts/k1122e_md_v3`, with output roots ending in
`_v3`; v2 inputs and pilots remain separate. GPU pilots cover primary WT,
primary K1122E and padded WT, each with the full 1.2 ns equilibration and 0.1 ns
production. The WT pilot also validates all ten cap combinations.
`verify_pilot.py` requires finite states, 290-330 K production temperatures,
positive volumes, OPC particles and ff19SB CMAPs. The production array is
submitted with `afterok` on all three pilots; final analysis/archive depends
on all 30 production tasks. Pilots never enter scientific contrasts.
See `data/interim/k1122e_md/TRUBA_jobs_20261004.json` for verified job IDs/status.
Historical v2 settings are in `protocol_ff14sb_tip3p_v2_20261004.json`.

This isolated tail fragment includes water and ions. It contains no microtubule,
motor domain, cargo partner or dimer; results address local structure only.

## Verification

```text
python -m unittest discover -s tests -p test_k1122e_md.py -v
```

Engine tests run when OpenMM is installed; dependency-free accounting tests still
run without it. `validate_caps.py` checks all exact AF-derived mutant systems.
