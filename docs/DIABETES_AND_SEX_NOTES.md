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
| Undamped Newton near the solution (`NEPH_FASTAMP`, IMCD only) | `Newton.py` | IMCD solve speed |
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
| `NEPH_POSSAFE` / `NEPH_POSSEGS` / `NEPH_POSFRAC` / `NEPH_POSVOL` | 1 / CCD / 0.1 / 0 | first-cell positivity safeguard; `NEPH_POSVOL=1` also keeps compartment volumes positive (was needed once for the rat female Moderate CCD with the older settings) |
| `NEPH_IMCD_LISBATH` | 1 | IMCD LIS starts at bath |
| `NEPH_FASTAMP` / `NEPH_FASTSEGS` | 1 / IMCD | undamped Newton near the solution |
| `NEPH_STRICT`, `NEPH_TRACEK`, `NEPH_RECX` | 0 | diagnostics |

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
`ENaC.py` (CNT): `facFv = max(0.01, 1+3*(flow/NaMPq0 - 1))`, `NaMPq0 = CNT inlet volume - reference flow` (mouse 1.6, rat 2.0 nl/min).
The reference flow was tuned on the healthy CNT inlet (mouse male 1.86 nl/min -> factor 19.5 at the inlet, ~5 at the outlet).  The Moderate mouse CNT
inlet is only 1.12 (male) / 1.30 (female) nl/min, below 1.6, so `NaMPq0 < 0`: the factor was pinned at the 0.01 floor everywhere, ENaC was effectively
off and the CNT secreted Na by the paracellular route (20.0 -> 38.6 pmol/min, male).  Fix: `set_params.enac_cnt_flowref_scale` multiplies the reference
flow by (diabetic/healthy CNT inlet flow) = 0.60 (male), 0.64 (female); rat 1.0 (inlet flow unchanged in diabetes); human 1.0.  Factor profile becomes
19.3 -> 9.9 (healthy 19.5 -> 5.2).  Result (male, CNT only): Na 20.0 -> 1.9 pmol/min (healthy mouse: 47.1 -> 5.6; 88% reabsorbed vs 91%); scale 0.55 gives 2.6, so it is
not knife-edge, though the factor is sensitive to the denominator when scale -> 0.7.  Full chain, Moderate mouse: urine 0.154 (was 0.158) male, 0.145 (0.149) female;
IMCD Na 62.1 (69.1) / 65.6 (77.2) mM.  `NEPH_ENAC_CNT_DIAB_SCALE=1` restores the legacy behaviour.  The CCD uses the same form (ref 0.08 nl/min-type constant) and was not changed.
