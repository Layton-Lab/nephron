# Notes on the rodent diabetes / sex-difference / collecting-duct changes

Scope: rat and mouse, male and female, **non-diabetic (`Non`) and Moderate diabetes**, superficial nephron, PT → IMCD.
Branch `local/numpy-asmatrix-patch` (changes relative to upstream `main` at `c5906b8`).
**Not re-validated:** Severe diabetes, SGLT2 inhibition and empagliflozin in diabetes, juxtamedullary/multiple-nephron runs,
pregnancy, unilateral nephrectomy. The human case keeps its legacy parameters (see "Human" below).

## 1. What changed and why

### 1.1 Collecting-duct robustness (CCD / OMCD / IMCD)
| Change | File | Why |
|---|---|---|
| First-cell positivity safeguard, per-cell solver statistics (`Newton.STATS`) | `Newton.py` | CCD first cell stalled/diverged; need per-cell diagnostics |
| OMCD first cell re-solved with the safeguard only if the plain solve has negative formate | `Newton.py` | valid OMCD solutions have tiny negative H2CO2, so the safeguard cannot be always on there |
| Undamped Newton near the solution (IMCD, built into `_damping`) | `Newton.py` | IMCD solve speed |
| Convective flux kept when `|dC| <= 1e-6` (`NEPH_CONVFIX`) | `electrochemical.py` | log-mean → common concentration, so the flux tends to `(1-sig)*C*Jv`, not 0; omitting it made the residual discontinuous |
| IMCD basolateral Na–HPO4 cotransport ×0.05 (rat and mouse, M and F) | `datafiles/IMCDparams_*` | non-physical IMCD roots |
| IMCD LIS concentrations start at bath values (`NEPH_IMCD_LISBATH`, default 1) | `driver.py` | IMCD params file has LIS = cell column |
| IMCD geometry for `N != 200` | `Newton.py`/`equations.py` | N-independent geometry |
| `np.mat` → `np.asmatrix` | `glucose.py` | NumPy 2 |

### 1.2 Diabetes: SGLT, glucose, GFR, torque, diameter (all in `set_params.py` / `boundaryBath.py`)
1. **Duplicate SGLT scaling removed.** The diabetic factors used to be applied in `set_params.py` *and* inside `glucose.sglt1/sglt2`.
   The `glucose.py` copies were removed. Human values are unchanged (the product of the two former factors).
2. **SGLT1 unchanged (×1.0) in rodents.** Literature is mixed (Akita mouse ↓, ob/ob mouse ↑, human T2D biopsies unchanged).
3. **Plasma glucose** – `boundaryBath.diab_plasma_glucose(cell)` is the single source of truth and is used for **both**
   the PT inlet (filtrate) and the bath/interstitial profile (which must agree): rat/mouse Moderate **25 mM** (Layton, Vallon & Edwards 2016),
   human Moderate 8.6 mM, Severe 25 mM (rodent) / 20 mM (human), healthy 5 mM. The bath used to be 5× healthy (25 mM) for every species.
4. **Diabetic mouse SNGFR** = healthy × `MOUSE_DIAB_SNGFR_FACTOR` (1.25; Akita GFR +25 %, micropuncture SNGFR +11 %, males). The
   mouse block used to contain the *rat* values (×4.4–5.5 the healthy mouse).
5. **Torque reference (bug fix).** PT/S3 apical fluxes are scaled by torque relative to a reference flow. For diabetic rodents the
   reference was the cell's own diabetic inlet flow, cancelling the hyperfiltration effect the papers describe. `torq_ref_vol()` now
   returns the **healthy** inlet flow (`cell.vol_torq0`). Switch: `NEPH_TORQ_HEALTHYREF=0` restores the old behaviour.
6. **Diabetic PT/S3 diameter ×1.2** (25 → 30 µm, Layton 2016) instead of ×1.1. Switch: `NEPH_PTDIAM`.
7. **Extra diabetic SGLT2 factor** on top of the papers' +38 %, per (species, sex): **rat M 4.2, rat F 3.4, mouse M 4.9, mouse F 3.55**
   (`SGLT2_DIAB_GLUCOSE_FACTOR_BY`). Tuned (with items 3, 5, 6) so that ≈20 % of the filtered glucose leaves S3
   (paper: fractional glucose excretion 17.9 % at 25 mM, SNGFR 30 → 45 nl/min).
   *This factor is large (total SGLT2 ≈ ×5–7 baseline in diabetes) because our PT glucose capacity is lower than the 2015/2016 PT
   model's; see "Open issues".*

### 1.3 Sex differences in baseline PT SGLT2 (`datafiles/PTparams_F_*.dat`)
Female SGLT2 density relative to male (male = 0.52e-9): **rat ×1.6 (0.832e-9)** (was ×2.5), **mouse ×1.9 (0.988e-9)** (was ×1.0).
Chosen so that fractional PT glucose reabsorption is similar in the sexes (≈10–14 % remaining half-way along the female rat PT) and PT Na/water
at the PT outlet stay close to the previous values.

## 2. Environment switches (all default to the behaviour described above)
| Variable | Default | Effect |
|---|---|---|
| `NEPH_TORQ_HEALTHYREF` | 1 | torque reference = healthy inlet flow (0 = legacy) |
| `NEPH_PTDIAM` | 1.2 | diabetic Moderate PT/S3 diameter factor |
| `NEPH_SGLT2_DIAB` | 1.0 | extra SGLT2 factor for (species, sex) pairs not in the table (use for scans after emptying the table) |
| `NEPH_GLU_ND` | 0 (off) | test hook: plasma glucose (mM) of the *non-diabetic* case, e.g. 25 |
| `NEPH_CONVFIX` | 1 | convective-flux continuity fix (0 = old) |
| `NEPH_POSSAFE` | 1 | positivity safeguard of the Newton step (0 = off); the segments it applies to are the constants `POS_FIRST_CELL` (CCD, S3), `POS_EVERY_CELL` (CCD), `POS_IF_NEG_FORMATE` (OMCD) in `Newton.py` |
| `NEPH_IMCD_LISBATH` | 1 | IMCD LIS starts at bath |
| `NEPH_STRICT` | 0 | raise instead of continuing when the Newton iteration cap is reached |

## 3. Results at the committed settings (superficial nephron, 200 cells/segment, PT → IMCD)
No errors and no non-converged cells in any of the eight runs below.

Urine (IMCD outlet), Moderate vs non-diabetic:

| | Rat M | Rat F | Mouse M | Mouse F |
|---|---|---|---|---|
| Urine volume, D/ND | 2.1 | 1.9 | 2.0 | 1.7 |
| Urine volume, % of PT inflow, ND → D | 0.92 → 1.30 | 1.23 → 1.54 | 0.76 → 1.24 | 1.05 → 1.46 |
| Glucose excreted, % of filtered | 19.9 | 20.6 | 19.7 | 19.4 |
| Urine glucose, mM | 385 | 334 | 399 | 332 |
| Urine Na, mM, ND → D | 151 → 65 | 140 → 54 | 145 → 69 | 148 → 77 |
| Na flow, D/ND | 0.9 | 0.7 | 1.0 | 0.9 |
| K flow, D/ND | 1.5 | 1.8 | 1.4 | 1.4 |
| Na excreted, % of filtered, ND → D | 0.96 → 0.58 | 1.18 → 0.57 | 0.76 → 0.59 | 1.07 → 0.78 |

Non-diabetic urine (mM): rat M Na 151 / K 95 / Cl 99 / urea 236; rat F Na 140 / K 83 / Cl 45 / urea 255;
mouse M Na 145 / K 121 / Cl 55 / urea 201; mouse F Na 148 / K 123 / Cl 42 / urea 196.

## 4. Open issues / caveats
* **Extra SGLT2 factor is large.** Even with the papers' SGLT2 (1.3e-9, ×2.5 our male baseline) the model at 25 mM / SNGFR 30 leaves ~66 % of the
  glucose after S3, versus ~48 % in Layton et al. 2015 (Fig 4A). ~1.7× PT glucose capacity is unexplained (candidates: apical/basolateral
  area, S3 area, torque scaling constants). The literature supports only +40–80 % SGLT2 in diabetic rodents.
* **Water:** at the tuned glucose handling the diabetic PT+S3 take up *more* water (as a fraction of inflow) than the non-diabetic ones
  (rat M: 0.165 vs 0.237 left after S3), although SNGFR is ×1.5, so absolute S3 outflow is ~10 % higher.
* **Urine volume** is ≈2× non-diabetic. It was ×6–10 before items 5–7 above. A published value of ≈4× (rat, Layton 2016) was recalled but
  could not be confirmed from the text.
* **Na excretion falls** in the Moderate runs (Na flow D/ND 0.7–1.0); osmotic diuresis would normally raise it. Unverified against data.
* **Urine HCO3 flow** is 2–5× higher in diabetes (not investigated).
* **Not re-validated** with the new rodent settings: Severe diabetes (uses the same SGLT2 table, torque reference and 25 mM), SGLT2
  inhibition and empagliflozin in diabetes (they multiply by the same factors), juxtamedullary nephrons, unilateral nephrectomy.
* Female mouse CCD limit cycle near cell 49 was seen with older settings; not seen in the committed runs.

## 5. Human
Human diabetes keeps its legacy SGLT factors (SGLT2 ×1.38, SGLT1 (1−0.33)·0.67·3.5, no re-tuning, no torque/diameter change).
**One change affects it:** the Moderate bath glucose is now the PT-inlet value 8.6 mM (it was 5 × 5 = 25 mM). Human results should be re-checked
before being used.

## 6. How to reproduce
```
cd nephron
python3 parallel_simulate.py --sex Male --species rat --type superficial --diabetes Moderate      # normal entry point
```
For a single-segment chain (what was used for the numbers above), call `compute_sup_jux_segment.compute_segment('sup', 'Male',
'rat', 'superficial', 'Moderate', None, 'N', 'non', 'N', 'results')` from a working directory that contains `datafiles/`, `outlets/` and
`results/`; the segments run in the order PT, S3, SDL, mTAL, cTAL, DCT, CNT and the outlet of each is written to `outlets/`.
The collecting duct (CCD, OMCD, IMCD) is then run with `driver.compute(N, 'datafiles/<SEG>params_<M|F>_<rat|mou>.dat', 'Newton', ...)`.
Outlet files hold 15 solute concentrations (rows 0–14, glucose = 14), then volume (row 15; ×6000 gives nl/min), electrical potential (16) and pressure (17).

## CNT ENaC flow dependence in diabetes (mouse)
`ENaC.py` (CNT): `facFv = max(floor, 1+3*(flow/NaMPq0 - 1))`, `NaMPq0 = CNT inlet volume - reference flow` (mouse 1.6, rat 2.0 nl/min).
Healthy mouse male: inlet 1.86 nl/min -> NaMPq0 = 0.26 nl/min, factor ~19.5 at the inlet and ~5 at the outlet.  The Moderate mouse CNT inlet is only 1.12 (male) /
1.30 (female) nl/min, below the 1.6 reference, so with the run's own inlet `NaMPq0 < 0`: the factor was pinned at the 0.01 floor everywhere, ENaC was effectively off
and the CNT secreted Na paracellularly (20.0 -> 38.6 pmol/min, male).
Fix (diabetic mouse only; `set_params.enac_cnt_healthy_inlet_vol`, `enac_cnt_floor`): NaMPq0 uses the HEALTHY CNT inlet (1.860 male / 2.023 female nl/min), a fixed
positive denominator, so the factor depends on the absolute CNT flow as in the healthy nephron; the floor is raised from 0.01 to 1.0 (healthy factor never goes below ~5), so
ENaC cannot switch off at low flow.  Healthy, rat (CNT inlet unchanged in diabetes) and human: legacy.  `NEPH_ENAC_CNT_LEGACY=1` restores the old formula, `NEPH_ENAC_CNT_FLOOR` sets the floor.
Cap: the diabetic-mouse factor is also limited to its healthy value at the healthy inlet, 1+3*1.6/(inlet-1.6) = 19.5 male / 12.3 female (`NEPH_ENAC_CNT_CAP=0` disables).  Without it the
high CNT flows of the SGLT2-inhibitor runs (3.4-4.6 nl/min) gave factors of 35-50, drove CNT lumen Na down to 1.8 mM (empagliflozin run) and the CCD Newton solve then diverged at cell 1.
The cap never binds for CNT inlet flows below ~2 nl/min (baseline unchanged).
An earlier attempt (reference flow x fixed 0.60) was dropped: it only works if the diabetic inlet is ~0.6 x healthy; A4/A6 (inlet 0.12/0.69 nl/min) secreted Na again.
CNT-only results (mouse male, Moderate, Na pmol/min in -> out): baseline 20.0 -> 2.8 (was 38.6; healthy mouse 47.1 -> 5.6); A6 (Na-K-ATPase x2) 14.7 -> 4.0 (was 32.5);
A4 (SGLT2 x2, inlet 0.12 nl/min, outlier) 7.1 -> 3.8 (was 18.3); B1 (SGLT2 inhibited) 96.6 -> 20.5 (was 53.6).  Full-chain urine with this version not yet rerun
(with the earlier 0.60 version: 0.154 male / 0.145 female nl/min vs 0.158 / 0.149 before).  CCD uses the same formula and was not changed.

## Newton positivity safeguard: S3 first cell and every CCD cell (mouse, Moderate)
Two solver failures in the parameter-sensitivity runs, both caused by a full Newton step driving a very small concentration (~3e-5 mM) negative:
* S3, halved PT Na-K-ATPase (`naka=0.5`): OverflowError in `ATPase.hatpase` in the FIRST S3 cell (the PT outlet composition is far from the S3 file's non-lumen start values).
* CCD, SGLT2 inhibition + NHE3 x0.35 (CNT inlet flow 4.6 nl/min): divergence at cell 1 (|f| 252 -> 7e4 at iteration 3, minimum concentration negative).
Fix (`Newton.py`): the fraction-to-boundary step limit (NEPH_POSSAFE=1) is now also applied in the first cell of the S3 (`POS_FIRST_CELL` = CCD, S3) and in every CCD cell
(`POS_EVERY_CELL` = CCD).  Where a case already solved the results are identical (checked: baseline and healthy mouse S3, baseline CCD: same flow, Na, iteration counts).
Results now: A5 (naka 0.5) urine 0.301 nl/min, B3 0.734 nl/min, both 0 non-converged.  Still failing: A4 (sglt2 x2; CCD diverges in its first cell, CNT inflow 0.12 nl/min).
B4 (empagliflozin) completes with CCD up to 290 iterations and 2 non-converged OMCD cells (CNT outlet lumen Na only 2.1 mM).
NOTE: a CNT-only diagnostic must pass the same `inhibition` as the full run: empagliflozin changes ROMK (DCT/CNT), paracellular permeabilities (CNT), AQP2, NHE3, NKCC2 ... in
several segments, which is why a CNT-only test without it predicted 19.4 instead of 5.1 pmol/min Na for B4.

## Newton.py clean-up (no change in results)
The rat/mouse solver was rewritten for readability: the per-segment damping if-chain (with ~100 lines of dead or commented-out branches: the unreachable IMCD `male` branch, LDL `iter>100`, duplicate
mTAL branches, a commented-out 70-line CNT schedule) is now one function `_damping(segment, k, iter, residual)`; the positivity safeguard is configured by three segment sets instead of five environment
variables; unused experiments (`NEPH_CHORD`, `NEPH_RECX`/`NEPH_TRACEK` debug recording, `NEPH_POSVOL`, `NEPH_FASTAMP`/`NEPH_FASTSEGS`, `NEPH_POSSAFE=2`) and unused variables were removed.  The damping VALUES are
unchanged (hand-tuned per segment).  Checks: (1) `_damping` equals the old if-chain on 4.3 million combinations of segment, sex, cell type, cell index, iteration (0-319) and residual norm (every threshold +-1e-9);
(2) full PT->IMCD re-runs of rat male Moderate, mouse male `Non`, B3 (SGLT2 inhibition + NHE3 x0.35) and A5 (Na-K-ATPase x0.5) give bit-identical summaries (all outlet values, flows, iteration counts) to the previous code.
`newton_human` and `Newton_preg.py` were not touched (no human or pregnancy validation was run).

## Newton.py: adaptive damping replaces the hand-tuned tables (rat / mouse)

`NEPH_SOLVER=adaptive` is now the default for `newton_rat`; `NEPH_SOLVER=tables` restores the previous per-segment, per-species
damping (`_damping`) for comparison.  The damped step is `x -= amp * J^-1 F` as before; only `amp` changes:

* the step grows x1.5 when ||f|| fell since the previous iteration, halves when it rose (bounds 0.05 ... 1); no extra function evaluations;
* retry cascade instead of segment tables: the first attempt starts with a full step; if a solve diverges, gives NaN, does not converge
  within 80 iterations, or ends on an unphysical root (a concentration below -0.05), the cell is re-solved from the same start with a
  smaller first step (first cell 1.0, 0.3, 0.12, 0.04; later cells 1.0, 0.3, 0.1).  The last attempt gets the full budget.
  Over the 8 reference cases and the 12 stress cases the cascade was needed for only a handful of cells (mTAL/cTAL/DCT first cell, a few CNT cells).
* A monotone Armijo line search on ||f|| was tried and does NOT work here: in the CCD it stalls at local minima of ||f||.  The unphysical-root
  check matters: without it the rat male non-diabetic OMCD converges (to 1e-9) to a root with glucose down to -55 mM.
* The positivity safeguards (POS_* sets) are unchanged.

Validation (urine = IMCD outlet, against the validated table-damping runs): the 8 reference cases (mouse/rat x M/F x Non/Moderate) agree to <= 1.1e-4
in urine flow and solute concentrations, 0 non-converged cells; mouse Moderate sensitivity runs A1-A8 and B1-B3 reproduce the previous outlet flows
(4 digits), B4 (empagliflozin) is unchanged (the same 2 non-converged OMCD cells), A4 (the known outlier) now runs through.  Iterations drop strongly
(IMCD max ~165 -> 8-34, OMCD 33 -> 9-12) and the OMCD takes less than half the time.  Results are not bit-identical to the table solver.
Not validated with it: Severe diabetes, human (`newton_human` is unchanged), `Newton_preg.py`.

## Calcium and magnesium (rat only)

Ca (solute 15) and Mg (solute 16), both valence +2, are ported from the Hakimi et al. calcium/magnesium code (2024, renal Ca and Mg handling in pregnancy): TRPV5, NCX1, PMCA (Ca-ATPase), TRPM6, NaMgX (Na:Mg 2:1), paracellular Ca/Mg, CaSR effects (TAL paracellular Ca/Mg permeability, ROMK, NKCC2, NCC, IMCD water permeability; CNT ENaC), interstitial Ca/Mg gradients.

* Switches: `NEPH_CALCIUM=1/0`, `NEPH_MAGNESIUM=1/0`; default on iff `rat` is on the command line. Mouse/human run with both off, as before. `NEPH_CASR_OMCD_HATPASE=1` enables the OMCD H-ATPase CaSR branch (dead code in the reference, default off).
* Segments without transcellular transport of a divalent keep the unknown in the vector but pinned (`x - x0 = 0`), which equals deleting it. Mg: whole DCT solves cell Mg (TRPM6); other cells pinned.
* Pregnancy: `set_params_preg.py` applies the Hakimi 2024 Table 1 style ratios (PMCA/NCX1 1.4/1.3, TRPV5, TRPM6, NaMgX superficial and juxtamedullary, PT/S3 apical Ca permeability 1.4/1.2). `NEPH_PREG_TABLE1=1` replaces the code ratios by the published Table 1 values (default off). Pregnancy runs use the adaptive Newton solver (`NEPH_SOLVER=tables` restores `Newton_preg`).
* Known reference quirks kept: TRPM6 uses cell pH[0] which is never updated (constant); our S3params_F_rat.dat has a duplicate NHE3_female line that the Mg repo removed.

Validation
* Virgin female rat Ca+Mg vs the original Mg code: within 0.7% through cTAL, 2-4% in DCT/CNT.
* Ca/Mg off reproduces previous results (mouse regression).
* Mid pregnancy vs the original Mg code (change MP/virgin in segment net reabsorption): PT Ca +25% (ref +28%, paper +33%), PT Mg +18% (+20%, +25%), mTAL Ca -13% (-18%, paper -11%), cTAL Ca -19% (-22%, -19%), cTAL Mg +15% (+15%, +18%), DCT Ca +22% (+34%), DCT Mg +96% (+102%). DCT/CNT concentrations differ from the reference run by 12-46%; parameters were verified identical and the difference starts at 1-3% in cTAL (PT outflow +1.1%), so it is likely solver/convergence path (the reference has unconverged CNT cells), not unresolved physics.
* Not validated: late pregnancy, CCD/OMCD/IMCD Ca/Mg in pregnancy, diabetes with Mg, inhibition options, unx, juxtamedullary nephrons, flux-decomposition outputs.

Pregnancy vs Stadt & Layton 2022: PT inlet loads match; distal and urinary values differ (e.g. NP urine Na 41 vs 32, MP 54 vs 34; PT Na reabsorption 54% vs 58%). The pristine upstream code shows the same differences (it has since gained glucose/SGLT etc.), and applying Table 1 values does not remove them. PT Na transport agrees with Hakimi 2024 (virgin 44.7, MP 53.8 umol/min per kidney).

Known limitations
* CCD cells may fail to converge with Ca (also with the old tuned damping).
* Diabetic rat: outlet Ca can come out slightly negative (over-reabsorption); needs investigation.
