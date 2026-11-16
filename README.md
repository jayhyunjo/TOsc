# TOsc — Wire-Cell 3+1 Sterile Neutrino Oscillation Fitter

C++/ROOT package implementing the MicroBooNE Wire-Cell 3+1 sterile-neutrino
oscillation analysis (BNB + NuMI, 7 channels × 2 beams), behind the published
Nature result. Original author: **Xiangpan Ji** (winxji).

This repository (`jayhyunjo/TOsc`) is a fork of
[winxji/ToAll_basic_nospeedup](https://github.com/winxji/ToAll_basic_nospeedup).

## Branch map

| Branch | Content |
|---|---|
| `xiangpan-final` (**default**) | Xiangpan's final version (via [jpmendez79/TOsc](https://github.com/jpmendez79/TOsc) commit `ec53810`, "Initial commit of new version from xiangpan", 2025-10). Loads **all 5 NuMI periods** (run1 FHC/RHC, run2 FHC/RHC, run3 RHC) + BNB runs 1–3. Use this. |
| `frozen-2022-local` | The original 2022 GitHub snapshot + local build fixes. Loads NuMI **run1 FHC only** — inconsistent with the full dataset (NuMI appearance template only ~19% covered). Kept for reference; do not use for physics. |
| `master` | Xiangpan's untouched `init` commit (Oct 2022). |

## What TOsc does

`TOsc` is a single class providing:
- **χ² calculation** (`FCN`) — full systematic covariance (flux, geant4, Xs,
  detector, dirt, MC-stat) with CNP statistical errors,
- **Minuit2 minimization** (`Minimization_OscPars_FullCov`),
- **pseudo-experiment generation** (`Set_toy_variations`, `Set_toy2fitdata`),
- goodness-of-fit utilities (`Exe_Goodness_of_fit`).

Oscillation effects are applied **event-by-event** from event lists (true
energy + baseline per event), not by histogram scaling: the CV contributions of
all oscillatable components are subtracted (`Set_oscillation_base`) and added
back weighted by `Prob_oscillaion(Etrue, baseline, mode)` at the requested
oscillation parameters (`Apply_oscillation`).

### Parameter convention (IMPORTANT)

`Set_oscillation_pars(dm2_41, p2, p3, p4)` in the **active** convention
("dm2 vs. t14 vs. t24", see `TOsc::Prob_oscillaion` in `src/TOsc.cxx`):

- `dm2_41` — Δm²₄₁ [eV²]
- `p2` — **sin²θ₁₄** (despite the variable name `sin2_2theta_14`!).
  Internally: effective sin²2θ₁₄ = 4·p2·(1−p2).
- `p3` — **sin²θ₂₄**
- `p4` — sin²θ₃₄ (unused, θ₃₄ = 0)

Probabilities implemented (θ₃₄ = 0):
- ν_e → ν_e : 1 − sin²2θ₁₄ · sin²Δ
- ν_μ → ν_μ : 1 − 4·cos²θ₁₄·sin²θ₂₄·(1 − cos²θ₁₄·sin²θ₂₄) · sin²Δ
- ν_μ → ν_e : sin²2θ₁₄ · sin²θ₂₄ · sin²Δ   (appearance amplitude sin²2θ_μe = sin²2θ₁₄·sin²θ₂₄)
- NC modes: survival of active flavors with the corresponding ν_s mixing.

Two alternative input conventions — (dm2, sin²2θ_μe, sin²θ₂₄) and
(dm2, sin²2θ_ee, sin²θ₂₄) — are present as commented blocks at the top of
`Prob_oscillaion`; swap the active block to use them for appearance /
disappearance-style scans.

## Requirements

- ROOT 6 **with Minuit2** (check: `root-config --features | grep minuit2`, or
  verify `libMinuit2.so` exists in `$ROOTSYS/lib`)
- g++ with C++17 (must be ABI/glibc-compatible with the ROOT build — see
  machine notes below)

## Build

```bash
export ROOTSYS=/path/to/root            # Makefile reads ROOTSYS from the environment
make clean
make                                     # produces ./read_oscillation_v01
```

### Machine-specific notes (jjo workstation)

Two ROOT installs exist; the conda one (first on PATH) cannot link the system
ROOT. Always build with:

```bash
cd /home/jjo/Work/MicroBooNE/ToAll_basic_nospeedup
PATH=/usr/local/bin:/usr/bin:$PATH ROOTSYS=/home/jjo/root-build make
```

The binary rpaths `/home/jjo/root-build/lib`, so it runs without sourcing
anything.

## Input files

`inc/WCPLEEANA/Configure_Osc.h` defines all input paths, by default under
`./data_inputs/dataset_input/`:

| Config entry | File | Content |
|---|---|---|
| `default_cv_file` | `merge.root` | Measured data (`hdata_obsch_1-14`, 7 BNB + 7 NuMI channels × 26 bins) + CV predictions (`histo_1-42`: per beam — 7 from regular MC, 7 from EXT, 7 from the ν_μ→ν_e full-oscillation sample) + `matrix_transform` (1092→364 collapse) |
| `default_dirtadd_file` | `merge.root` | `cov_mat_add` (1092², dirt systematics) |
| `default_mcstat_file` | `0.log` | MC-stat uncertainty per bin (from a no-oscillation framework run) |
| `default_fluxXs_dir` | `hist_rootfiles/XsFlux_edit/` | `cov_1..13` flux, `cov_14..16` geant4, `cov_17` Xs — fractional covariances `frac_cov_xf_mat_N`, **1092×1092** (42 histos × 26 bins), NuMI-knob hack applied |
| `default_detector_dir` | `hist_rootfiles/DetVar_edit/` | 9 detector variations `frac_cov_det_mat_N`, 1092², appearance-template blocks filled from overlay (see its `read_obj.cc`) |
| `default_eventlist_dir` | top level | 24 event lists `roofile_obj_{BNB,NuMI}_*.root` (trees of Etrue/Ereco/weight/baseline per component) + POT normalization files `checkout_prodgenie_*.root`, `run*_data_*.root` |

### Option A — use the prepared dataset (this workstation)

The canonical dataset is `/home/jjo/Work/MicroBooNE/LEEana/BNBNuMI_TOsc_input_combined/`
(see its README.md for provenance — assembled from two handover copies, the
`numifluxpatched` CV taking precedence; verified 2026-11-16 to match the
published Nature spectra exactly, bin by bin, in all four numuCC channels —
the other snapshot's NuMI CV is a pre-flux-patch version, +30–35% off).
Wire it up with a single symlink:

```bash
mkdir -p data_inputs
ln -s /home/jjo/Work/MicroBooNE/LEEana/BNBNuMI_TOsc_input_combined data_inputs/dataset_input
```

### Option B — regenerate inputs from the LEEana framework

Summary of Xiangpan's recipe (docs: "TOsc Instruction.pdf",
"TOsc_preparation.rtf", data structure in docdb-38700 p12–18; needs the
oscillation-modified wcp-uboone-bdt framework with full-oscillation samples):

1. Configure with fullosc samples, **breakdown mode** in `file_ch.txt`;
   `./convert_hist` (default mode) → `checkout_*.root`; `merge_hist` → `merge.root`.
2. Reconfigure **without** fullosc, non-breakdown, MC-stat mode
   (`cov_input.txt` + `file_ch.txt`); `./convert_hist`; `./run_mcstat.pl` → `0.log`.
3. Back to fullosc + breakdown; switch `convert_checkout_hist.cxx` to roofile
   mode (overlays only); produce each BNB roofile one-by-one via `cv_input.txt`
   (BNB mode), then each NuMI roofile (NuMI mode) → `roofile_obj_*.root`.
4. Run `xf_sys` (flux/Xs covariances, must come out 1092² = 42 histos × 26 bins),
   then apply `XsFlux_edit/read_obj.cc`: zeroes BNB×NuMI cross blocks and keeps
   only knob 3 for the NuMI block ("one knob represents all NuMI knobs") → `XsFlux_edit/`.
5. Run `det_sys`, then `DetVar_edit/read_obj.cc`: copies overlay covariance
   blocks into the appearance-template blocks → `DetVar_edit/`.

## Run

```bash
./read_oscillation_v01            # initialization + whatever example blocks are enabled
```

CLI options: `-pbnb F` / `-pnumi F` (POT scale factors), `-d 1` (interactive
graphics), `-f N`, `-it14/-idm2/-it24 N` (job indices for grid scans — parsed
but unused by the stock driver; intended for batch scan drivers).

The driver `read_oscillation_v01.cxx` is a template: initialization is always
run, then `if( 0 )`-gated example blocks show the API —

```cpp
osc_test->Set_oscillation_pars(dm2, s2t14, s2t24, 0);
osc_test->Apply_oscillation();          // event-by-event reweight
osc_test->Set_apply_POT();              // POT scale + collapse 1092 -> 364
osc_test->Set_meas2fitdata();           // fit real data,  OR
osc_test->Set_asimov2fitdata();         // fit the Asimov set built from current pars
osc_test->Set_toy_variations(N);        // or generate N pseudo-experiments
osc_test->Set_toy2fitdata(i);           //   and fit the i-th one
double chi2 = osc_test->FCN(pars4);     // chi2 at pars4 (internally re-applies oscillation)
osc_test->Minimization_OscPars_FullCov(dm2, t14, t24, 0, "");  // Minuit2 fit ("" = all free;
                                        //  "dm2"/"t14"/"t24" substrings fix parameters)
```

Flip a block to `if( 1 )`, `make`, run.

## Validation (2026-11-16, combined dataset)

- Loads with **zero ROOT errors** (17 systematic matrices, all 1092²).
- No-oscillation closure (self-check block near `src/TOsc.cxx:1657`,
  uncomment to use; it prints once per `Apply_oscillation` — including those
  inside `FCN` — so analyze only the first 546-line block): regular bins close
  to machine zero in both beams; NuMI appearance residue after event-list
  subtraction < 0.003 events.
- Asimov χ² at generation point = 0 exactly.
- Reference values: Asimov@(1.19, 0.8, 0.015) vs (2, 0.3, 0.1) = **2618.71**;
  real data vs (1.19, 0.8, 0.015) = **420.451**.

## Performance (i9-10900X, 10 cores / 20 threads, 31 GB)

- Initialization (one-time per process): ~1.4 s
- One `FCN` evaluation (full event-by-event reweight + collapse + χ²): **0.16 s**
- Peak RSS per process: **0.46 GB** → ~20 parallel single-threaded jobs fit
  comfortably in RAM.

Rule of thumb: a 60×60 Wilks Δχ² grid ≈ 10 CPU-min (≈1 min on 10–20 workers);
a frequentist CLs scan with ~500 toys/point on a ~45×45 grid ≈ overnight on
this machine (matches the "~12 h" experience of the downstream fork).

## Contour production (nue appearance)

`scan_appearance_cls.cxx` computes the χ² grid for the appearance-only
("simple case") analysis: `TOsc::flag_osc_appearance_only` switches
`Prob_oscillaion` to pure ν_μ→ν_e (second oscillation parameter becomes
sin²2θ_μe directly, all survival probabilities 1). Grid: 60×60 log10 bin
centers, Δm²₄₁ ∈ [0.01, 100] eV², sin²2θ_μe ∈ [1e-4, 1].

```bash
make                                              # builds scan_appearance_cls too
for j in $(seq 0 9); do ./scan_appearance_cls -j $j -n 10 & done; wait   # ~2.5 min on 10 cores
root -l -b -q plot_scan_appearance.C              # Gaussian-CLs 95% contours + comparison
```

Per grid point the scan stores χ²(4ν) on the 3ν-Asimov sample (= Δχ²_Asimov)
and on data, plus χ²(3ν|data). `plot_scan_appearance.C` converts these to 95%
CLs contours with the Gaussian CLs approximation (arXiv:1407.5052; sensitivity
contour ⇔ Δχ²_Asimov = 3.84) and overlays the published Nature curves from
`natureReproducer/data/`. First result (2026-11-16): sensitivity and data
curves track the published ones closely on the low-Δm² arm (0.03–2 eV²);
deviations above ~3 eV² are expected since the published appearance plot
profiles the full 3+1 over sin²θ₂₄ and uses toy-based CLs. Outputs:
`contour_appearance_cls.{png,pdf,root}`.

## Related

- `jpmendez79/TOsc` (Jesse Mendez, LSU) — active downstream fork: distributed
  best-fit / χ²-grid tools, CLs machinery, invisible-decay branch. Registered
  here as git remote `jesse`.
- Local docs: `/home/jjo/Work/MicroBooNE/docs/` (TOsc Instruction.pdf,
  TOsc_preparation.rtf, WireCell_OSC_Collab.pdf, talks).
- Published-result curves for comparison: `/home/jjo/Work/MicroBooNE/natureReproducer/data/`.
