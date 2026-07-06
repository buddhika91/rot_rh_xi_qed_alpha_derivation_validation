# ROT RH / Xi-QED Alpha Derivation Validation

This repository contains a reproducibility script for a proposed Xi-QED derivation chain for the low-energy fine-structure constant inverse, (\alpha(0)^{-1}).

The main script is:

```text
rot_rh_xi_qed_alpha_derivation_validation.py
```

It validates the full numerical derivation path:

```text
Xi curvature invariant
→ base spectral count 2π/K₂
→ exact spacelike one-loop QED correction
→ muon threshold residual
→ Xi heat-curvature S-fraction damping
→ Fibonacci/binary lower-scale refinement
→ final predicted alpha inverse
```

The current result is a stable high-precision near-identity for (\alpha(0)^{-1}), reaching approximately (10^{-27}) relative error against the reference value used in the script.

This repository is intended for independent numerical validation, criticism, and theoretical follow-up.

---

## Important Status

This is **not claimed as a completed unconditional analytic proof**.

The script validates a highly structured numerical derivation and produces data files that researchers can inspect. The remaining analytic proof obligations are explicitly tracked in the generated validation tables.

Current status:

```text
Stable numerical near-identity: yes
Closed-form QED kernel: yes
Xi invariants frozen: yes
Precision stability: yes
Rounding stability: yes
Independent QED recomputation: yes

Unconditional analytic proof: not yet
Xi → QED boundary measure theorem: open
Heat-curvature S-fraction derivation: open
```

---

## Requirements

Python 3.10+ is recommended.

Install dependencies:

```bash
pip install mpmath
```

The script uses only standard-library modules plus `mpmath`.

---

## Quick Start

Run the main validation:

```bash
python rot_rh_xi_qed_alpha_derivation_validation.py \
  --dps 100 \
  --dps-list "60,80,100,120,160" \
  --rounding-digits-list "20,30,40,50,60" \
  --random-control-trials 5000 \
  --random-denom-min 1 \
  --random-denom-max 20 \
  --random-sample-rows 500 \
  --seed 314159 \
  --out-prefix xi_qed_alpha_derivation_validation
```

On Windows PowerShell:

```powershell
python rot_rh_xi_qed_alpha_derivation_validation.py `
  --dps 100 `
  --dps-list "60,80,100,120,160" `
  --rounding-digits-list "20,30,40,50,60" `
  --random-control-trials 5000 `
  --random-denom-min 1 `
  --random-denom-max 20 `
  --random-sample-rows 500 `
  --seed 314159 `
  --out-prefix xi_qed_alpha_derivation_validation
```

Optional deeper Xi derivative check:

```powershell
python rot_rh_xi_qed_alpha_derivation_validation.py `
  --dps 100 `
  --recompute-xi-low-orders `
  --out-prefix xi_qed_alpha_derivation_validation
```

---

## What the Script Produces

The script generates only validation data files, not manuscript or email files.

Typical outputs include:

```text
xi_qed_alpha_derivation_validation_axioms.csv
xi_qed_alpha_derivation_validation_xi_invariants.csv
xi_qed_alpha_derivation_validation_qed_kernel_validation.csv
xi_qed_alpha_derivation_validation_qed_channel_components.csv
xi_qed_alpha_derivation_validation_correction_factors.csv
xi_qed_alpha_derivation_validation_coefficient_refinement.csv
xi_qed_alpha_derivation_validation_derivation_stages.csv
xi_qed_alpha_derivation_validation_precision_stability.csv
xi_qed_alpha_derivation_validation_rounding_stability.csv
xi_qed_alpha_derivation_validation_random_control_sample.csv
xi_qed_alpha_derivation_validation_proof_status.csv
```

Depending on options, JSON summary/control files may also be written for machine-readable validation.

---

## Core Definitions

The central Xi curvature invariant is:

```math
K_2
=
-\left.\frac{d^2}{d\gamma^2}
\log \Xi\!\left(\frac12+i\gamma\right)
\right|_{\gamma=0}.
```

The base Xi spectral count is:

```math
A_0=\frac{2\pi}{K_2}.
```

The script uses the frozen numerical value:

```text
K2 ≈ 0.0462099862308379415778676208606780280067635
```

which gives:

```text
A0 ≈ 135.970291698259424392746169290913234...
```

The reference value used for validation is:

```text
alpha(0)^(-1) = 137.035999177
```

So the required gap is:

```math
\Delta_\alpha
=
\alpha(0)^{-1}
-
\frac{2\pi}{K_2}.
```

---

## Exact QED Kernel

The QED correction uses the exact spacelike one-loop vacuum-polarization kernel.

For a lepton mass (m) and spacelike scale (Q), define:

```math
a=\left(\frac{Q}{m}\right)^2.
```

The kernel is:

```math
I(a)
=
-\frac{5}{18}
+
\frac{2}{3a}
+
\frac{a-2}{3a}
\frac{\operatorname{atanh}\sqrt{a/(a+4)}}{\sqrt{a/(a+4)}}.
```

The two frozen scales are:

```math
Q_1=\frac{4}{\pi}m_\pi,
```

and

```math
Q_2=\sqrt{2\pi}\frac{m_\mu m_\pi}{m_\mu+m_\pi}.
```

The second scale uses the electron and muon one-loop contributions.

---

## Xi Heat-Curvature Invariants

The normalized Xi heat-curvature invariants are:

```math
R_{2n}
=
\frac{L_{2n}}{K_2^n},
```

where (L_{2n}) are even derivatives of the local Xi log-profile.

The script freezes:

```text
R4  ≈ -0.208897141829072608475473247120663...
R6  ≈ -0.350663421125273269083846578271025...
R8  ≈ -1.465721490477072574803784790016129...
R10 ≈ -11.0691512339896275471965503980481...
```

The heat-curvature damping factors are:

```math
\frac{1}{1+K_2},
```

```math
\frac{1}{1+R_4/(6\pi^2)},
```

```math
\frac{1}{1+R_6/(360\pi^3)},
```

```math
\frac{1}{1+R_8/(2520\pi^5)},
```

and

```math
\frac{1}{1+R_{10}/(113400\pi^8)}.
```

The phase powers:

```text
2, 3, 5, 8
```

follow the Fibonacci pattern.

The denominator branch is mostly consistent with:

```math
D_{2k}=\frac{(2k)!}{2^k},
```

except for the isolated (R_6) anomaly:

```math
D_6=4\cdot \frac{6!}{2^3}=360.
```

---

## Lower-Scale S-Fraction Refinement

The current strongest near-closure uses an effective heat scale:

```math
K_2^{\mathrm{eff}}=K_2(1+s).
```

The correction is:

```math
s
=
-\frac{b_4^2}{2+K_2+c\,b_4},
```

where:

```math
b_4=\frac{R_4}{6\pi^2},
```

and:

```math
b_8=\frac{R_8}{2520\pi^5}.
```

The current best frozen coefficient is:

```math
c
=
\frac85
-
\frac{K_2}{16}
+
\frac{b_8}{7+4b_4+K_2^2/22}.
```

This coefficient is not fitted inside the script. It is a frozen symbolic candidate from the previous audit chain.

The script validates whether this frozen expression reproduces the alpha inverse value stably.

---

## Final Prediction Formula

The prediction has the form:

```math
\alpha_{\mathrm{pred}}^{-1}
=
\frac{2\pi}{K_2}
+
\Delta_{\mathrm{QED}}^{e+\mu}(Q_2)
-
\frac{\operatorname{Res}_\mu(Q_2)K_2^{\mathrm{eff}}}{4\pi^2}
\lambda_{\mathrm{heat}}.
```

Here:

```math
\lambda_{\mathrm{heat}}
=
\frac{1}{1+K_2^{\mathrm{eff}}}
\frac{1}{1+R_4/(6\pi^2)}
\frac{1}{1+R_6/(360\pi^3)}
\frac{1}{1+R_8/(2520\pi^5)}
\frac{1}{1+R_{10}/(113400\pi^8)}.
```

The predicted value is compared against the reference (\alpha(0)^{-1}) inside the generated validation tables.

---

## Validation Philosophy

This repository separates three different levels of evidence:

### 1. Direct computation

The script computes:

```text
Xi curvature count
QED correction
muon residual
heat-curvature damping
S-fraction refinement
final alpha prediction
```

### 2. Numerical stability

The script checks:

```text
multi-precision stability
rounding stability of Xi invariants
independent QED kernel recomputation
random denominator/control behavior
```

### 3. Proof obligations

The script explicitly marks which steps are:

```text
proved analytically
verified numerically
conditionally supported
still open
```

This is important because the current result is a high-precision near-identity, not a completed proof.

---

## Main Open Problems

The validation script tracks the two main remaining theorem obligations.

### Open Problem 1 — Xi Boundary Measure Selection

The QED Feynman parameter measure is equivalent to a Jacobian-weighted circular projection:

```math
x=\frac{1+\cos\theta}{2},
```

so:

```math
x(1-x)=\frac{\sin^2\theta}{4},
```

and:

```math
dx=\frac{\sin\theta}{2}d\theta.
```

The open theorem is to prove that the Xi heat/phase boundary measure selects (dx), or equivalently the Jacobian phase measure, from Xi principles alone.

### Open Problem 2 — Heat-Curvature S-Fraction Law

The current finite correction structure is:

```math
K_2,\quad R_4,\quad R_6,\quad R_8,\quad R_{10}.
```

The open theorem is to derive the full S-fraction law that produces:

```text
Fibonacci phase powers
factorial denominators
the R6 factor-4 anomaly
the lower-scale coefficient 8/5 - K2/16 + b8/(7+4b4+K2^2/22)
```

without relying on reverse-engineered numerical matching.

---

## Interpreting the Result

A successful run means:

```text
The frozen Xi-QED formula is numerically stable and reproducible.
```

It does not mean:

```text
The fine-structure constant has been unconditionally proved from Xi.
```

The strongest accurate claim is:

```text
This script validates a stable Xi-QED high-precision near-identity for alpha inverse,
with a transparent chain of numerical derivation and explicit analytic proof obligations.
```

---

## Suggested Repository Layout

```text
.
├── README.md
├── rot_rh_xi_qed_alpha_derivation_validation.py
├── data/
│   └── optional generated validation outputs
└── LICENSE
```

Recommended usage:

```text
Commit the script and README.
Generate validation data locally.
Commit selected data files only if needed for reproducibility snapshots.
```

---

## Reproducibility Notes

Use a fixed seed for random controls:

```text
--seed 314159
```

Use at least 100 decimal digits:

```text
--dps 100
```

Use multi-precision validation:

```text
--dps-list "60,80,100,120,160"
```

Use rounding validation:

```text
--rounding-digits-list "20,30,40,50,60"
```

A robust run should show that the final near-identity is stable under increased precision and under rounding of the frozen Xi invariants beyond approximately 50 digits.

---

## License

Choose a license appropriate for your intended use.

For open mathematical/scientific reproducibility, MIT, Apache-2.0, or BSD-3-Clause are common choices.

---

## Contact / Purpose

This repository is intended for peer review, criticism, replication, and theoretical development of the Xi-QED alpha derivation program.

The most valuable contributions are:

```text
1. independent reproduction of the numerical validation,
2. attacks on the QED/Xi boundary-measure bridge,
3. derivation or falsification of the heat-curvature S-fraction law,
4. identification of hidden fitting or non-canonical assumptions,
5. conversion of the strongest numerical identities into rigorous theorems.
```
