#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ROT RH / Xi-QED — Alpha Derivation Validation Data Generator
============================================================

This script is intended for GitHub/researcher validation of the current
Xi-QED alpha-derivation chain.  It does NOT generate a manuscript, email,
README, or narrative report.  It only writes machine-readable validation data
files: CSV and JSON.

What this script validates
--------------------------
The script reconstructs the full current numerical chain:

  1. Xi origin:
       K2 = - d^2/dγ^2 log Xi(1/2+iγ)|_{γ=0}
       A0 = 2π/K2

  2. QED boundary kernel:
       I(a) = ∫_0^1 x(1-x) log(1+a x(1-x)) dx

     with the exact spacelike one-loop closed form used in the audits.

  3. Muon threshold residual:
       epsilon0 = - Res_mu(Q2) K2/(4π^2)

  4. Heat-curvature/S-fraction damping:
       b4  = R4 /(6π^2)
       b6  = R6 /(360π^3)
       b8  = R8 /(2520π^5)
       b10 = R10/(113400π^8)

       lambda_heat = Π 1/(1+b_j)

  5. Final lower-scale renormalization:
       c  = 8/5 - K2/16 + b8/(7 + 4 b4 + K2^2/22)
       s  = -b4^2/(2 + K2 + c b4)
       K2_eff = K2(1+s)

  6. Final alpha prediction:
       alpha_pred^{-1} = A0 + Delta2 + epsilon0 * lambda_final

       where lambda_final = lambda_heat/(1+K2_eff).

Important honesty note
----------------------
This script validates a highly stable, structured numerical near-identity.  It
is NOT a formal proof of the fine-structure constant and it does not prove the
open theorem obligations.  The remaining theorem obligations are exported as a
CSV data table for reviewers:

  - Xi boundary-measure selection theorem: open.
  - Heat-curvature S-fraction theorem: open.
  - Last-mile coefficient residual derivation: open.

The goal is reproducibility and falsifiability: every numerical component is
written to data files so that a researcher can independently inspect the chain.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import random
from datetime import datetime, timezone
from typing import Any, Dict, Iterable, List, Tuple

import mpmath as mp

LINE = "=" * 120
DASH = "-" * 120

# -------------------------------------------------------------------------------------------------
# Frozen constants
# -------------------------------------------------------------------------------------------------
# These constants are intentionally frozen.  The purpose of this validation script is to evaluate
# the stated theory, not to search or refit parameters.  Values are taken from the closed-form
# kernel and final-formula stability audits.

FROZEN_STR: Dict[str, str] = {
    # Xi heat-curvature invariants.  K2 is the curvature scale and R_{2n} are normalized log-Xi
    # derivatives R_{2n}=L_{2n}/K2^n, where L_{2n}=d^{2n}/dγ^{2n} log Xi(1/2+iγ)|_{γ=0}.
    "K2": "0.0462099862308379415778676208606780280067635207948441802463650011215275",
    "R4": "-0.208897141829072608475473247120663288512610605",
    "R6": "-0.350663421125273269083846578271025228260426921",
    "R8": "-1.46572149047707257480378479001612943977596606",
    "R10": "-11.0691512339896275471965503980481068492534248",
    "R12": "-130.177921348738620887201271845299268907653484",

    # Reference alpha inverse used in the complete audit chain.
    "alpha_obs_inv": "137.035999177",

    # Physical masses in MeV.  These are not varied in this script.
    "m_e_MeV": "0.51099895000",
    "m_mu_MeV": "105.6583755",
    "m_pi_charged_MeV": "139.57039",
}

# Denominators and phase powers from the current heat-curvature law.
HEAT_DENOMS = {"b4": 6, "b6": 360, "b8": 2520, "b10": 113400}
HEAT_PI_POWERS = {"b4": 2, "b6": 3, "b8": 5, "b10": 8}


# -------------------------------------------------------------------------------------------------
# Generic output helpers
# -------------------------------------------------------------------------------------------------

def mpstr(x: Any, n: int = 50) -> str:
    """Return a stable decimal string for mpmath values."""
    if isinstance(x, mp.mpf) or isinstance(x, mp.mpc):
        return mp.nstr(x, n)
    return str(x)


def parse_int_list(text: str) -> List[int]:
    return [int(t.strip()) for t in str(text).split(",") if t.strip()]


def protocol_hash(payload: Dict[str, Any]) -> str:
    """Hash the validation protocol so outputs are tied to a reproducible configuration."""
    blob = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def write_csv(path: str, rows: List[Dict[str, Any]], fields: List[str] | None = None) -> None:
    """Write rows to CSV.  Values are converted to strings before writing."""
    if fields is None:
        fields = sorted({k for row in rows for k in row.keys()}) if rows else []
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: mpstr(row.get(k, ""), 90) for k in fields})


def write_json(path: str, obj: Any) -> None:
    """Write JSON data.  mpmath values are converted through mpstr."""
    def default(o: Any) -> str:
        if isinstance(o, (mp.mpf, mp.mpc)):
            return mpstr(o, 90)
        return str(o)

    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, sort_keys=True, default=default)


def round_sig_mpf(x: mp.mpf, digits: int) -> mp.mpf:
    """Round an mpf to a requested number of significant decimal digits."""
    if x == 0:
        return mp.mpf("0")
    return mp.mpf(mp.nstr(x, int(digits), strip_zeros=False))


# -------------------------------------------------------------------------------------------------
# Xi function utilities
# -------------------------------------------------------------------------------------------------
# These are included so researchers can optionally recompute the lowest Xi derivatives directly.
# High derivatives beyond R8 are expensive and therefore kept frozen by default.

def xi_completed(s: mp.mpc) -> mp.mpc:
    """Completed Riemann Xi function in the standard normalization."""
    return mp.mpf("0.5") * s * (s - 1) * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def recompute_low_xi_invariants(dps: int) -> List[Dict[str, Any]]:
    """Optionally recompute K2, R4, R6, R8 from numerical differentiation of log Xi."""
    mp.mp.dps = int(dps)
    xi0 = xi_completed(mp.mpf("0.5"))

    def L(gamma: mp.mpf) -> mp.mpf:
        z = xi_completed(mp.mpf("0.5") + 1j * gamma) / xi0
        return mp.re(mp.log(z))

    K2 = -mp.diff(L, mp.mpf("0"), 2)
    rows: List[Dict[str, Any]] = []
    rows.append({
        "quantity": "K2_recomputed",
        "definition": "-d^2 log Xi(1/2+i gamma)/d gamma^2 at gamma=0",
        "value": K2,
        "frozen_value": FROZEN_STR["K2"],
        "relative_error_vs_frozen": abs(K2 - mp.mpf(FROZEN_STR["K2"])) / abs(mp.mpf(FROZEN_STR["K2"])),
    })
    for order, key in [(4, "R4"), (6, "R6"), (8, "R8")]:
        derivative = mp.diff(L, mp.mpf("0"), order)
        R = derivative / (K2 ** (order // 2))
        rows.append({
            "quantity": key + "_recomputed",
            "definition": f"d^{order} log Xi / K2^{order//2}",
            "value": R,
            "frozen_value": FROZEN_STR[key],
            "relative_error_vs_frozen": abs(R - mp.mpf(FROZEN_STR[key])) / abs(mp.mpf(FROZEN_STR[key])),
        })
    return rows


# -------------------------------------------------------------------------------------------------
# QED one-loop kernel utilities
# -------------------------------------------------------------------------------------------------
# The exact spacelike one-loop kernel is a central part of the chain.  The script validates the
# closed form against numerical quadrature and against the Jacobian phase substitution.

def qed_kernel_integral_closed(a: mp.mpf) -> mp.mpf:
    """Closed form for I(a)=∫_0^1 x(1-x) log(1+a x(1-x)) dx."""
    a = mp.mpf(a)
    if a == 0:
        return mp.mpf("0")
    if abs(a) < mp.mpf("1e-8"):
        # Stable power-series fallback for extremely small a.
        total = mp.mpf("0")
        for n in range(1, 400):
            term = ((-1) ** (n + 1)) * (a ** n) * mp.beta(n + 2, n + 2) / mp.mpf(n)
            total += term
            if abs(term) < mp.eps * max(1, abs(total)):
                break
        return total
    s = mp.sqrt(a / (a + 4))
    return -mp.mpf(5) / 18 + mp.mpf(2) / (3 * a) + ((a - 2) / (3 * a)) * (mp.atanh(s) / s)


def qed_kernel_integral_x_quad(a: mp.mpf) -> mp.mpf:
    """Direct x-space integral for I(a)."""
    f = lambda x: x * (1 - x) * mp.log(1 + a * x * (1 - x))
    return mp.quad(f, [0, 1])


def qed_kernel_integral_theta_jacobian_quad(a: mp.mpf) -> mp.mpf:
    """Phase-variable integral using x=(1+cos θ)/2 and dx=(sin θ/2)dθ."""
    f = lambda th: (mp.sin(th) ** 2 / 4) * mp.log(1 + a * (mp.sin(th) ** 2 / 4)) * (mp.sin(th) / 2)
    return mp.quad(f, [0, mp.pi])


def qed_one_loop_spacelike(Q: mp.mpf, m: mp.mpf) -> mp.mpf:
    """Dimensionless one-loop spacelike vacuum-polarization contribution used by the chain."""
    return (mp.mpf(2) / mp.pi) * qed_kernel_integral_closed((Q / m) ** 2)


def qed_one_loop_asymptotic(Q: mp.mpf, m: mp.mpf) -> mp.mpf:
    """Asymptotic high-Q approximation used only to define the finite muon residual."""
    return (mp.mpf(2) / (3 * mp.pi)) * (mp.log(Q / m) - mp.mpf(5) / 6)


def reduced_mass(a: mp.mpf, b: mp.mpf) -> mp.mpf:
    return a * b / (a + b)


# -------------------------------------------------------------------------------------------------
# Core derivation math
# -------------------------------------------------------------------------------------------------

def constants() -> Dict[str, mp.mpf]:
    return {k: mp.mpf(v) for k, v in FROZEN_STR.items()}


def compute_qed_channel(c: Dict[str, mp.mpf]) -> Dict[str, mp.mpf]:
    """Compute QED correction, muon residual, epsilon0, and diagnostic lambda_exact."""
    K2 = c["K2"]
    alpha_obs_inv = c["alpha_obs_inv"]
    A0 = 2 * mp.pi / K2
    delta_needed = alpha_obs_inv - A0

    # The Q2 scale is the frozen two-channel reduced-mass scale used throughout the audits.
    Q2 = mp.sqrt(2 * mp.pi) * reduced_mass(c["m_mu_MeV"], c["m_pi_charged_MeV"])

    q2_e = qed_one_loop_spacelike(Q2, c["m_e_MeV"])
    q2_mu = qed_one_loop_spacelike(Q2, c["m_mu_MeV"])
    q2_mu_asymp = qed_one_loop_asymptotic(Q2, c["m_mu_MeV"])
    q2_mu_resid = q2_mu - q2_mu_asymp

    Delta2 = q2_e + q2_mu
    epsilon_target = delta_needed - Delta2
    epsilon0 = -q2_mu_resid * K2 / (4 * mp.pi ** 2)
    lambda_exact = epsilon_target / epsilon0

    return {
        "A0_xi_curvature_count": A0,
        "alpha_obs_inv": alpha_obs_inv,
        "delta_needed": delta_needed,
        "Q2_MeV": Q2,
        "q2_electron": q2_e,
        "q2_muon_exact": q2_mu,
        "q2_muon_asymptotic": q2_mu_asymp,
        "q2_muon_residual": q2_mu_resid,
        "Delta2_QED_e_plus_mu": Delta2,
        "epsilon_target": epsilon_target,
        "epsilon0": epsilon0,
        "lambda_exact_diagnostic": lambda_exact,
    }


def heat_components(c: Dict[str, mp.mpf]) -> Dict[str, mp.mpf]:
    """Compute heat-curvature b-terms and S-fraction/final coefficient corrections."""
    K2, R4, R6, R8, R10 = c["K2"], c["R4"], c["R6"], c["R8"], c["R10"]

    b4 = R4 / (6 * mp.pi ** 2)
    b6 = R6 / (360 * mp.pi ** 3)
    b8 = R8 / (2520 * mp.pi ** 5)
    b10 = R10 / (mp.mpf(113400) * mp.pi ** 8)
    heat_product = 1 / ((1 + b4) * (1 + b6) * (1 + b8) * (1 + b10))

    # The Fibonacci-binary refinement found in the last-mile residual audit.
    c0 = mp.mpf(8) / 5 - K2 / 16
    dc = b8 / (7 + 4 * b4 + K2 ** 2 / 22)
    c_final = c0 + dc
    denominator = 2 + K2 + c_final * b4
    s_final = -b4 ** 2 / denominator
    K2_eff = K2 * (1 + s_final)

    lambda_K2 = 1 / (1 + K2)
    lambda_R4 = lambda_K2 / (1 + b4)
    lambda_R6 = lambda_R4 / (1 + b6)
    lambda_R8 = lambda_R6 / (1 + b8)
    lambda_R10 = lambda_R8 / (1 + b10)
    lambda_c0 = heat_product / (1 + K2 * (1 - b4 ** 2 / (2 + K2 + c0 * b4)))
    lambda_final = heat_product / (1 + K2_eff)

    return {
        "b4": b4,
        "b6": b6,
        "b8": b8,
        "b10": b10,
        "heat_product_R4_to_R10": heat_product,
        "lambda_K2": lambda_K2,
        "lambda_R4": lambda_R4,
        "lambda_R6": lambda_R6,
        "lambda_R8": lambda_R8,
        "lambda_R10": lambda_R10,
        "c0": c0,
        "dc": dc,
        "c_final": c_final,
        "s_final": s_final,
        "K2_eff": K2_eff,
        "lambda_c0": lambda_c0,
        "lambda_final": lambda_final,
    }


def alpha_from_lambda(q: Dict[str, mp.mpf], lam: mp.mpf) -> mp.mpf:
    """Alpha inverse prediction at a particular damping lambda."""
    return q["A0_xi_curvature_count"] + q["Delta2_QED_e_plus_mu"] + q["epsilon0"] * lam


def stage_rows(c: Dict[str, mp.mpf], q: Dict[str, mp.mpf], h: Dict[str, mp.mpf]) -> List[Dict[str, Any]]:
    """Create the full derivation-stage table from Xi base count to final prediction."""
    target = c["alpha_obs_inv"]
    rows: List[Dict[str, Any]] = []

    def add(stage: int, name: str, method: str, lam: mp.mpf | None, counted: bool, notes: str) -> None:
        if lam is None:
            pred = q["A0_xi_curvature_count"]
        else:
            pred = alpha_from_lambda(q, lam)
        residual = pred - target
        rows.append({
            "stage": stage,
            "name": name,
            "method": method,
            "counted_in_final_formula": counted,
            "lambda_value": "" if lam is None else lam,
            "alpha_inverse_predicted": pred,
            "alpha_inverse_observed": target,
            "absolute_residual": residual,
            # Two relative errors are exported because prior audits used the Xi-QED gap
            # delta_needed as the denominator, while general readers may expect alpha^{-1}.
            "relative_error_vs_alpha_inverse": abs(residual) / abs(target),
            "relative_error_vs_xi_qed_gap": abs(residual) / abs(q["delta_needed"]),
            "relative_error": abs(residual) / abs(q["delta_needed"]),
            "notes": notes,
        })

    # Stage 0: pure Xi curvature count.
    add(0, "Xi curvature count", "A0=2*pi/K2", None, True, "Origin from Xi curvature invariant only.")

    # Stage 1: QED e+mu exact channel without finite muon residual correction.
    pred_qed = q["A0_xi_curvature_count"] + q["Delta2_QED_e_plus_mu"]
    rows.append({
        "stage": 1,
        "name": "Add exact QED e+mu channel",
        "method": "A0+Delta2",
        "counted_in_final_formula": True,
        "lambda_value": "",
        "alpha_inverse_predicted": pred_qed,
        "alpha_inverse_observed": target,
        "absolute_residual": pred_qed - target,
        "relative_error_vs_alpha_inverse": abs(pred_qed - target) / abs(target),
        "relative_error_vs_xi_qed_gap": abs(pred_qed - target) / abs(q["delta_needed"]),
        "relative_error": abs(pred_qed - target) / abs(q["delta_needed"]),
        "notes": "Closed-form spacelike QED kernel at frozen Q2 scale.",
    })

    # Stage 2 and onward: residual correction with increasingly refined damping lambda.
    add(2, "Muon residual undamped", "A0+Delta2+epsilon0", mp.mpf(1), True, "Finite muon threshold residual with no heat damping.")
    add(3, "K2 damping", "lambda=1/(1+K2)", h["lambda_K2"], True, "First heat-scale S-fraction factor.")
    add(4, "R4 quartic factor", "lambda=lambda_K2/(1+b4)", h["lambda_R4"], True, "Quartic heat-curvature correction.")
    add(5, "R6 sextic factor", "lambda=lambda_R4/(1+b6)", h["lambda_R6"], True, "Sextic correction; denominator includes isolated R6 factor 4 anomaly.")
    add(6, "R8 octic factor", "lambda=lambda_R6/(1+b8)", h["lambda_R8"], True, "Fibonacci phase power pi^5 and denominator 2520.")
    add(7, "R10 decic factor", "lambda=lambda_R8/(1+b10)", h["lambda_R10"], True, "Finite R10 heat-curvature lock.")
    add(8, "Fibonacci-binary c0 scale", "c0=8/5-K2/16", h["lambda_c0"], True, "Lower-factor renormalization using Fibonacci ratio and binary quartic scale.")
    add(9, "Final frozen formula", "c=c0+b8/(7+4*b4+K2^2/22)", h["lambda_final"], True, "Best stable frozen near-closure; no target coefficient used.")
    add(10, "Diagnostic exact lambda", "lambda_exact=(delta_needed-Delta2)/epsilon0", q["lambda_exact_diagnostic"], False, "Diagnostic only; uses observed alpha and is not counted in final formula.")
    return rows


# -------------------------------------------------------------------------------------------------
# Validation data section builders
# -------------------------------------------------------------------------------------------------

def axiom_rows() -> List[Dict[str, Any]]:
    """Axioms/theorem assumptions as data so reviewers can see which parts are open."""
    return [
        {
            "id": "A1",
            "axiom_or_theorem": "Xi curvature invariant",
            "statement": "K2=-d^2 log Xi(1/2+iγ)/dγ^2|0 defines the canonical Xi heat scale; A0=2π/K2.",
            "validation_in_script": "K2 and A0 exported; optional low-order Xi derivative recomputation available.",
            "status": "NUMERICALLY_VALIDATED_LOW_ORDER; ANALYTIC_STANDARD_DEFINITION",
        },
        {
            "id": "A2",
            "axiom_or_theorem": "Boundary measure selection",
            "statement": "Xi local heat/phase boundary selects Lebesgue dx, hence QED Feynman measure via x=(1+cosθ)/2.",
            "validation_in_script": "QED x-integral and theta-Jacobian integral are compared numerically.",
            "status": "CONDITIONAL; XI_SELECTION_THEOREM_OPEN",
        },
        {
            "id": "A3",
            "axiom_or_theorem": "Closed-form QED kernel",
            "statement": "Exact spacelike one-loop kernel I(a) is used, not a fitted log approximation.",
            "validation_in_script": "Closed form is compared with numerical quadrature for requested a-list.",
            "status": "NUMERICALLY_VALIDATED",
        },
        {
            "id": "A4",
            "axiom_or_theorem": "Heat-curvature S-fraction law",
            "statement": "Residual damping is generated by K2 and b4,b6,b8,b10 S-fraction factors.",
            "validation_in_script": "All factors and staged improvements are exported.",
            "status": "HIGH_PRECISION_NUMERICAL_STRUCTURE; FULL_ANALYTIC_DERIVATION_OPEN",
        },
        {
            "id": "A5",
            "axiom_or_theorem": "Fibonacci-binary lower-scale refinement",
            "statement": "c=8/5-K2/16+b8/(7+4b4+K2^2/22) renormalizes K2 to K2_eff.",
            "validation_in_script": "Coefficient and final alpha prediction are exported and stress-tested.",
            "status": "STABLE_HIGH_PRECISION_NEAR_CLOSURE; EXACT_THEOREM_OPEN",
        },
    ]


def xi_invariant_rows(c: Dict[str, mp.mpf], q: Dict[str, mp.mpf], h: Dict[str, mp.mpf]) -> List[Dict[str, Any]]:
    return [
        {"quantity": "K2", "value": c["K2"], "definition": "Xi heat curvature", "role": "base scale"},
        {"quantity": "A0=2*pi/K2", "value": q["A0_xi_curvature_count"], "definition": "Xi curvature count", "role": "base alpha inverse approximation"},
        {"quantity": "R4", "value": c["R4"], "definition": "normalized Xi fourth log derivative", "role": "quartic heat correction"},
        {"quantity": "R6", "value": c["R6"], "definition": "normalized Xi sixth log derivative", "role": "sextic heat correction"},
        {"quantity": "R8", "value": c["R8"], "definition": "normalized Xi eighth log derivative", "role": "octic heat correction"},
        {"quantity": "R10", "value": c["R10"], "definition": "normalized Xi tenth log derivative", "role": "decic heat correction"},
        {"quantity": "b4", "value": h["b4"], "definition": "R4/(6*pi^2)", "role": "quartic S-fraction seed"},
        {"quantity": "b6", "value": h["b6"], "definition": "R6/(360*pi^3)", "role": "R6 anomaly corrected seed"},
        {"quantity": "b8", "value": h["b8"], "definition": "R8/(2520*pi^5)", "role": "last-mile residual seed"},
        {"quantity": "b10", "value": h["b10"], "definition": "R10/(113400*pi^8)", "role": "R10 lock seed"},
    ]


def qed_kernel_validation_rows(a_values: Iterable[str]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for text in a_values:
        a = mp.mpf(text)
        closed = qed_kernel_integral_closed(a)
        x_quad = qed_kernel_integral_x_quad(a)
        theta_quad = qed_kernel_integral_theta_jacobian_quad(a)
        rows.append({
            "a": a,
            "closed_form": closed,
            "x_quadrature": x_quad,
            "theta_jacobian_quadrature": theta_quad,
            "closed_vs_x_relerr": abs(closed - x_quad) / max(abs(closed), mp.mpf("1e-80")),
            "closed_vs_theta_relerr": abs(closed - theta_quad) / max(abs(closed), mp.mpf("1e-80")),
            "identity_checked": "x=(1+cos(theta))/2; dx=(sin(theta)/2)dtheta",
        })
    return rows


def qed_channel_rows(q: Dict[str, mp.mpf]) -> List[Dict[str, Any]]:
    return [
        {"quantity": "Q2_MeV", "value": q["Q2_MeV"], "definition": "sqrt(2*pi)*reduced_mass(mu,pi)"},
        {"quantity": "q2_electron", "value": q["q2_electron"], "definition": "electron one-loop contribution at Q2"},
        {"quantity": "q2_muon_exact", "value": q["q2_muon_exact"], "definition": "muon one-loop exact contribution at Q2"},
        {"quantity": "q2_muon_asymptotic", "value": q["q2_muon_asymptotic"], "definition": "muon high-Q asymptotic contribution at Q2"},
        {"quantity": "q2_muon_residual", "value": q["q2_muon_residual"], "definition": "exact minus asymptotic muon residual"},
        {"quantity": "Delta2_QED_e_plus_mu", "value": q["Delta2_QED_e_plus_mu"], "definition": "q2_electron+q2_muon_exact"},
        {"quantity": "epsilon0", "value": q["epsilon0"], "definition": "-q2_muon_residual*K2/(4*pi^2)"},
        {"quantity": "lambda_exact_diagnostic", "value": q["lambda_exact_diagnostic"], "definition": "diagnostic target damping from observed alpha; not counted"},
    ]


def correction_factor_rows(c: Dict[str, mp.mpf], h: Dict[str, mp.mpf]) -> List[Dict[str, Any]]:
    return [
        {
            "factor": "K2",
            "formula": "1/(1+K2)",
            "denominator": "1+K2",
            "phase_power": "",
            "value": 1 / (1 + c["K2"]),
            "origin_status": "canonical heat-scale S-fraction seed",
        },
        {
            "factor": "R4",
            "formula": "1/(1+R4/(6*pi^2))",
            "denominator": 6,
            "phase_power": 2,
            "value": 1 / (1 + h["b4"]),
            "origin_status": "factorial branch and Fibonacci phase power",
        },
        {
            "factor": "R6",
            "formula": "1/(1+R6/(360*pi^3))",
            "denominator": 360,
            "phase_power": 3,
            "value": 1 / (1 + h["b6"]),
            "origin_status": "R6 anomaly factor 4 isolated; exact analytic origin open",
        },
        {
            "factor": "R8",
            "formula": "1/(1+R8/(2520*pi^5))",
            "denominator": 2520,
            "phase_power": 5,
            "value": 1 / (1 + h["b8"]),
            "origin_status": "factorial denominator branch and Fibonacci phase power",
        },
        {
            "factor": "R10",
            "formula": "1/(1+R10/(113400*pi^8))",
            "denominator": 113400,
            "phase_power": 8,
            "value": 1 / (1 + h["b10"]),
            "origin_status": "factorial denominator branch and Fibonacci phase power",
        },
        {
            "factor": "lower_scale_refinement",
            "formula": "K2_eff=K2*(1-b4^2/(2+K2+c*b4))",
            "denominator": "2+K2+c*b4",
            "phase_power": "Fibonacci/binary coefficient c",
            "value": h["K2_eff"],
            "origin_status": "near-closure; exact S-fraction theorem open",
        },
    ]


def coefficient_refinement_rows(c: Dict[str, mp.mpf], h: Dict[str, mp.mpf], q: Dict[str, mp.mpf]) -> List[Dict[str, Any]]:
    # Diagnostic exact coefficient, not used in counted formula.
    heat_product = h["heat_product_R4_to_R10"]
    K2 = c["K2"]
    b4, b8 = h["b4"], h["b8"]
    lambda_exact = q["lambda_exact_diagnostic"]
    s_req = (heat_product / lambda_exact - 1) / K2 - 1
    D_req = -b4 ** 2 / s_req
    eta_req = D_req - (2 + K2)
    c_req = eta_req / b4
    dc_req = c_req - h["c0"]
    D8_req = b8 / dc_req

    candidates = [
        ("c0", "8/5-K2/16", h["c0"], True),
        ("dc", "b8/(7+4*b4+K2^2/22)", h["dc"], True),
        ("c_final", "c0+dc", h["c_final"], True),
        ("c_req_diagnostic", "eta_req/b4", c_req, False),
        ("dc_req_diagnostic", "c_req-c0", dc_req, False),
        ("D8_req_diagnostic", "b8/dc_req", D8_req, False),
    ]
    rows = []
    for name, formula, value, counted in candidates:
        rows.append({
            "name": name,
            "formula": formula,
            "value": value,
            "counted_in_final_formula": counted,
            "relative_error_vs_c_req": "" if name.startswith("d") or name.startswith("D8") else abs(value - c_req) / abs(c_req),
            "notes": "diagnostic uses observed alpha" if not counted else "predeclared counted expression",
        })
    return rows


def precision_rows(args: argparse.Namespace) -> List[Dict[str, Any]]:
    rows = []
    for dps in parse_int_list(args.dps_list):
        mp.mp.dps = dps
        c = constants()
        q = compute_qed_channel(c)
        h = heat_components(c)
        stages = stage_rows(c, q, h)
        final = stages[-2]  # Final frozen formula, before diagnostic exact lambda.
        rows.append({
            "dps": dps,
            "final_alpha_inverse_predicted": final["alpha_inverse_predicted"],
            "final_absolute_residual": final["absolute_residual"],
            "final_relative_error": final["relative_error"],
            "c_final": h["c_final"],
            "K2_eff": h["K2_eff"],
        })
    return rows


def rounding_rows(args: argparse.Namespace) -> List[Dict[str, Any]]:
    mp.mp.dps = int(args.dps)
    cfull = constants()
    qfull = compute_qed_channel(cfull)
    rows = []
    for digits in parse_int_list(args.rounding_digits_list):
        c = dict(cfull)
        for key in ["K2", "R4", "R6", "R8", "R10"]:
            c[key] = round_sig_mpf(c[key], digits)
        # Keep QED channel target fixed, but use rounded K2 in the formula path.
        q = compute_qed_channel(c)
        h = heat_components(c)
        final_pred = alpha_from_lambda(q, h["lambda_final"])
        residual = final_pred - c["alpha_obs_inv"]
        rows.append({
            "significant_digits": digits,
            "K2_rounded": c["K2"],
            "R4_rounded": c["R4"],
            "final_alpha_inverse_predicted": final_pred,
            "final_absolute_residual": residual,
            "final_relative_error_vs_alpha_inverse": abs(residual) / abs(c["alpha_obs_inv"]),
            "final_relative_error_vs_xi_qed_gap": abs(residual) / abs(q["delta_needed"]),
            "final_relative_error": abs(residual) / abs(q["delta_needed"]),
            "relative_difference_vs_full_final": abs(final_pred - alpha_from_lambda(qfull, heat_components(cfull)["lambda_final"])) / abs(alpha_from_lambda(qfull, heat_components(cfull)["lambda_final"])),
        })
    return rows


def random_control_rows(args: argparse.Namespace) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Random local controls around the last coefficient denominator.

    The counted last-mile denominator is D8=7+4b4+K2^2/22.  This control samples nearby
    random denominators in a broad range to see how often a random b8/D8 correction beats the
    counted final formula.  This is not a proof test; it is a falsification/specificity check.
    """
    mp.mp.dps = int(args.dps)
    rnd = random.Random(int(args.seed))
    c = constants()
    q = compute_qed_channel(c)
    h = heat_components(c)
    target = c["alpha_obs_inv"]
    baseline = alpha_from_lambda(q, h["lambda_final"])
    baseline_relerr = abs(baseline - target) / abs(q["delta_needed"])

    b4, b8, K2 = h["b4"], h["b8"], c["K2"]
    c0 = h["c0"]
    heat_product = h["heat_product_R4_to_R10"]
    epsilon0 = q["epsilon0"]
    A0_plus_delta2 = q["A0_xi_curvature_count"] + q["Delta2_QED_e_plus_mu"]

    rows = []
    count_better = 0
    best = None
    for i in range(int(args.random_control_trials)):
        # sample denominator uniformly around the natural range containing 7.
        D8 = mp.mpf(str(rnd.uniform(float(args.random_denom_min), float(args.random_denom_max))))
        coeff = c0 + b8 / D8
        s = -b4 ** 2 / (2 + K2 + coeff * b4)
        K2_eff = K2 * (1 + s)
        lam = heat_product / (1 + K2_eff)
        pred = A0_plus_delta2 + epsilon0 * lam
        relerr = abs(pred - target) / abs(q["delta_needed"])
        relerr_alpha = abs(pred - target) / abs(target)
        if relerr <= baseline_relerr:
            count_better += 1
        row = {
            "trial": i,
            "random_D8": D8,
            "random_coeff": coeff,
            "alpha_inverse_predicted": pred,
            "absolute_residual": pred - target,
            "relative_error_vs_alpha_inverse": relerr_alpha,
            "relative_error_vs_xi_qed_gap": relerr,
            "relative_error": relerr,
        }
        if best is None or relerr < best["relative_error"]:
            best = row
        if i < min(int(args.random_control_trials), int(args.random_sample_rows)):
            rows.append(row)

    summary = {
        "random_control_trials": int(args.random_control_trials),
        "random_denom_min": str(args.random_denom_min),
        "random_denom_max": str(args.random_denom_max),
        "baseline_final_relative_error": baseline_relerr,
        "count_better_or_equal": count_better,
        "p_better_or_equal_with_plus_one": mp.mpf(count_better + 1) / mp.mpf(int(args.random_control_trials) + 1),
        "best_random": best,
    }
    return rows, summary


def proof_status_rows() -> List[Dict[str, Any]]:
    return [
        {
            "item": "Xi derivative definitions",
            "status": "DEFINED_AND_OPTIONALLY_RECOMPUTABLE_LOW_ORDER",
            "what_would_upgrade_status": "Independent symbolic/numerical reproduction of K2,R4,R6,R8,R10 from Xi.",
        },
        {
            "item": "QED Feynman/Jacobian measure identity",
            "status": "NUMERICALLY_VALIDATED_AND_ANALYTIC_CHANGE_OF_VARIABLE_STANDARD",
            "what_would_upgrade_status": "Formal write-up of x=(1+cosθ)/2 pullback and compact moment determinacy.",
        },
        {
            "item": "Xi selects dx boundary measure",
            "status": "OPEN_THEOREM",
            "what_would_upgrade_status": "Prove Xi local heat/phase boundary is max-entropy/Haar in x.",
        },
        {
            "item": "R4/R6/R8/R10 heat-curvature law",
            "status": "HIGH_PRECISION_NUMERICALLY_SUPPORTED; FULL_ANALYTIC_ORIGIN_OPEN",
            "what_would_upgrade_status": "Derive denominators 6,360,2520,113400 and Fibonacci phase powers from an S-fraction theorem.",
        },
        {
            "item": "R6 anomaly factor 4",
            "status": "ISOLATED_BUT_ORIGIN_AMBIGUOUS",
            "what_would_upgrade_status": "Select one analytic source for the factor 4 and prove it is forced.",
        },
        {
            "item": "Fibonacci-binary coefficient refinement",
            "status": "STABLE_HIGH_PRECISION_NEAR_CLOSURE; EXACT_DERIVATION_OPEN",
            "what_would_upgrade_status": "Derive c=8/5-K2/16+b8/(7+4b4+K2^2/22) from Xi heat/S-fraction structure.",
        },
        {
            "item": "Final alpha identity",
            "status": "STABLE_NEAR_IDENTITY_NOT_PROOF",
            "what_would_upgrade_status": "Close both open theorem obligations without using observed alpha as input.",
        },
    ]


# -------------------------------------------------------------------------------------------------
# Main CLI
# -------------------------------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate CSV/JSON validation data for the ROT RH / Xi-QED alpha derivation chain."
    )
    parser.add_argument("--dps", type=int, default=100, help="mpmath decimal precision")
    parser.add_argument("--out-prefix", default="xi_qed_alpha_derivation_validation", help="prefix for data output files")
    parser.add_argument("--a-list", default="0.01,0.1,1,3", help="QED kernel a-values for validation")
    parser.add_argument("--dps-list", default="60,80,100,120,160", help="precision values for stability table")
    parser.add_argument("--rounding-digits-list", default="20,30,40,50,60", help="significant-digit rounding table")
    parser.add_argument("--random-control-trials", type=int, default=5000, help="random denominator controls")
    parser.add_argument("--random-denom-min", default="1", help="minimum denominator for random b8/D control")
    parser.add_argument("--random-denom-max", default="20", help="maximum denominator for random b8/D control")
    parser.add_argument("--random-sample-rows", type=int, default=500, help="number of random-control rows to write")
    parser.add_argument("--seed", type=int, default=314159)
    parser.add_argument("--recompute-xi-low-orders", action="store_true", help="slow: recompute K2,R4,R6,R8 directly from Xi")
    parser.add_argument("--protocol-only", action="store_true", help="write run metadata only and exit")
    args = parser.parse_args()

    mp.mp.dps = int(args.dps)

    protocol = {
        "script": os.path.basename(__file__),
        "time_utc_like": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "dps": args.dps,
        "out_prefix": args.out_prefix,
        "a_list": args.a_list,
        "dps_list": args.dps_list,
        "rounding_digits_list": args.rounding_digits_list,
        "random_control_trials": args.random_control_trials,
        "random_denom_range": [args.random_denom_min, args.random_denom_max],
        "seed": args.seed,
        "claim_type": "validation data generator; no fitting; no manuscript/email output",
        "final_formula": "c=8/5-K2/16+b8/(7+4*b4+K2^2/22); s=-b4^2/(2+K2+c*b4)",
    }
    protocol["protocol_hash"] = protocol_hash(protocol)

    # Always write metadata first.  This is a data file, not a narrative report.
    write_json(args.out_prefix + "_run_metadata.json", protocol)

    if args.protocol_only:
        print(LINE)
        print("ROT RH / Xi-QED alpha-derivation validation protocol written")
        print(LINE)
        print(f"metadata: {args.out_prefix}_run_metadata.json")
        print("No numerical validation run performed because --protocol-only was used.")
        return

    c = constants()
    q = compute_qed_channel(c)
    h = heat_components(c)
    stages = stage_rows(c, q, h)
    final_stage = stages[-2]  # stage 9 final frozen formula

    # Build all validation tables.
    data_tables: Dict[str, List[Dict[str, Any]]] = {
        "axioms": axiom_rows(),
        "xi_invariants": xi_invariant_rows(c, q, h),
        "qed_kernel_validation": qed_kernel_validation_rows([x.strip() for x in args.a_list.split(",") if x.strip()]),
        "qed_channel_components": qed_channel_rows(q),
        "correction_factors": correction_factor_rows(c, h),
        "derivation_stages": stages,
        "coefficient_refinement": coefficient_refinement_rows(c, h, q),
        "precision_stability": precision_rows(args),
        "rounding_stability": rounding_rows(args),
        "proof_status": proof_status_rows(),
    }

    if args.recompute_xi_low_orders:
        data_tables["xi_low_order_recompute"] = recompute_low_xi_invariants(int(args.dps))

    random_rows, random_summary = random_control_rows(args)
    data_tables["random_control_sample"] = random_rows

    # Write CSV data files.
    written_files: List[str] = [args.out_prefix + "_run_metadata.json"]
    for name, rows in data_tables.items():
        path = f"{args.out_prefix}_{name}.csv"
        write_csv(path, rows)
        written_files.append(path)

    # Write JSON data files.
    final_components = {
        "K2": c["K2"],
        "A0_xi_curvature_count": q["A0_xi_curvature_count"],
        "alpha_obs_inv_reference": c["alpha_obs_inv"],
        "Delta2_QED_e_plus_mu": q["Delta2_QED_e_plus_mu"],
        "epsilon0": q["epsilon0"],
        "lambda_final": h["lambda_final"],
        "alpha_inverse_predicted_final": final_stage["alpha_inverse_predicted"],
        "alpha_inverse_residual_final": final_stage["absolute_residual"],
        "alpha_inverse_relative_error_final_vs_gap": final_stage["relative_error_vs_xi_qed_gap"],
        "alpha_inverse_relative_error_final_vs_alpha_inverse": final_stage["relative_error_vs_alpha_inverse"],
        "c_final": h["c_final"],
        "s_final": h["s_final"],
        "K2_eff": h["K2_eff"],
        "protocol_hash": protocol["protocol_hash"],
    }
    write_json(args.out_prefix + "_final_components.json", final_components)
    write_json(args.out_prefix + "_random_control_summary.json", random_summary)

    validation_index = {
        "protocol_hash": protocol["protocol_hash"],
        "final_relative_error_vs_gap": final_stage["relative_error_vs_xi_qed_gap"],
        "final_relative_error_vs_alpha_inverse": final_stage["relative_error_vs_alpha_inverse"],
        "final_absolute_residual": final_stage["absolute_residual"],
        "claim_status": "stable high-precision near-identity; analytic proof obligations remain open",
        "data_files": written_files + [
            args.out_prefix + "_final_components.json",
            args.out_prefix + "_random_control_summary.json",
            args.out_prefix + "_validation_index.json",
        ],
        "no_markdown_or_email_files_generated": True,
    }
    write_json(args.out_prefix + "_validation_index.json", validation_index)
    written_files.extend([
        args.out_prefix + "_final_components.json",
        args.out_prefix + "_random_control_summary.json",
        args.out_prefix + "_validation_index.json",
    ])

    # Console output is concise; the real output is the data files above.
    print(LINE)
    print("ROT RH / Xi-QED — ALPHA DERIVATION VALIDATION DATA GENERATOR")
    print(LINE)
    print(f"protocol hash           : {protocol['protocol_hash']}")
    print(f"out_prefix              : {args.out_prefix}")
    print(f"final alpha^-1 predicted: {mpstr(final_stage['alpha_inverse_predicted'], 70)}")
    print(f"observed alpha^-1 ref   : {mpstr(c['alpha_obs_inv'], 70)}")
    print(f"final residual          : {mpstr(final_stage['absolute_residual'], 70)}")
    print(f"final relative error vs Xi-QED gap : {mpstr(final_stage['relative_error_vs_xi_qed_gap'], 70)}")
    print(f"final relative error vs alpha^-1  : {mpstr(final_stage['relative_error_vs_alpha_inverse'], 70)}")
    print(DASH)
    print("Data files written:")
    for path in written_files:
        print(f"  {path}")
    print(DASH)
    print("Audit status: stable high-precision near-identity; analytic proof obligations remain open.")
    print(LINE)


if __name__ == "__main__":
    main()
