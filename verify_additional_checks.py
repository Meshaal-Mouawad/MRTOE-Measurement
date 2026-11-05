#!/usr/bin/env python3
"""Reproduce offset identities and the quadratic-family condition.

Companion to verify_revision.py for the MRTOE methodological study.
All inputs below are explicitly constructed mathematical controls.
Run: python verify_additional_checks.py > additional_checks_results.json
"""
from __future__ import annotations
import json
from typing import Any
import numpy as np

TAU = 60.0
EPS = 0.05
K = 0.31
N = 600
DT = 0.1
OFFSET = 5.0
MC = K * EPS * np.sqrt(TAU)

def evaluate(x: np.ndarray, t: np.ndarray) -> dict[str, Any]:
    if x.ndim != 1 or t.shape != x.shape or len(x) < 7:
        raise ValueError('Provide equal one-dimensional arrays of at least seven values.')
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(t)) or not np.all(np.diff(t) > 0):
        raise ValueError('Values must be finite and timestamps strictly increasing.')
    w = max(2, min(int(10 * TAU), len(x) // 3))
    w += int(w % 2 == 0)
    p = w // 2
    weights = np.ones(w, dtype=float) / w
    mass = np.convolve(np.ones(len(x)), weights, mode='same')
    smooth = np.convolve(x, weights, mode='same')
    gradient = np.gradient(x, float(np.mean(np.diff(t))))
    gate = np.abs(gradient) > EPS
    integrand = (x - smooth) * gate / TAU
    interior = slice(p, len(x) - p)
    normalized = (x - smooth / mass) * gate / TAU
    shifted = x + OFFSET
    shifted_gate = np.abs(np.gradient(shifted, float(np.mean(np.diff(t))))) > EPS
    if not np.array_equal(gate, shifted_gate):
        raise AssertionError('Offset changes the gate for this input; the identity premise fails.')
    shifted_score = float(np.trapezoid(
        (shifted - np.convolve(shifted, weights, mode='same')) * shifted_gate / TAU, t))
    score = float(np.trapezoid(integrand, t))
    delta_direct = shifted_score - score
    delta_identity = OFFSET / TAU * float(np.trapezoid((1 - mass) * gate, t))
    if abs(delta_direct - delta_identity) >= 1e-10:
        raise AssertionError('Offset identity failed.')
    return {
        'n': len(x), 'window': w, 'half_window': p,
        'original_score': score,
        'normalized_full_score': float(np.trapezoid(normalized, t)),
        'normalized_interior_score': float(np.trapezoid(normalized[interior], t[interior])),
        'interior_gate_all_active': bool(np.all(gate[interior])),
        'strictly_decreasing': bool(np.all(np.diff(x) < 0)),
        'offset_delta_direct': delta_direct,
        'offset_delta_identity': delta_identity,
        'offset_identity_error': abs(delta_direct - delta_identity),
    }

def main() -> None:
    t = np.arange(N, dtype=float) * DT
    end = t[-1]
    curves = {
        'linear': np.linspace(50, 25, N),
        'exponential': 25 + 25 * np.exp(-np.arange(N) / 90),
        'quadratic': 50 - 25 * (t / end) ** 2,
        'constant': np.full(N, 26.0),
    }
    offsets = {label: evaluate(x, t) for label, x in curves.items()}
    p = offsets['quadratic']['half_window']
    mu2 = DT ** 2 * p * (p + 1) / 3
    length = t[-p - 1] - t[p]
    gate_bound = EPS / (2 * t[p])
    score_bound = float(MC * TAU / (mu2 * length))
    joint_bound = max(gate_bound, score_bound)
    coefficients = {
        'gate_pass_score_below_1': 0.003,
        'gate_pass_score_below_2': 0.005,
        'just_below_score_boundary': score_bound * (1 - 1e-4),
        'just_above_score_boundary': score_bound * (1 + 1e-4),
        'original_quadratic_control': 25 / end ** 2,
    }
    family = {}
    for label, b in coefficients.items():
        r = evaluate(50 - b * t ** 2, t)
        analytic = float(b * mu2 * length / TAU)
        assert r['strictly_decreasing'] and r['interior_gate_all_active']
        assert abs(analytic - r['normalized_interior_score']) < 1e-11
        assert (r['normalized_interior_score'] > MC) == (b > score_bound)
        if b > joint_bound:
            assert r['normalized_interior_score'] > MC
        family[label] = {
            'b_C_per_second_squared': float(b),
            'analytic_interior_score': analytic,
            'computed_interior_score': r['normalized_interior_score'],
            'interior_above_threshold': bool(r['normalized_interior_score'] > MC),
            'strictly_decreasing': r['strictly_decreasing'],
            'interior_gate_all_active': r['interior_gate_all_active'],
            'normalized_full_score': r['normalized_full_score'],
        }
    output = {
        'description': 'Read-only mathematical checks; original measured records are checked by verify_revision.py.',
        'parameters': {'n': N, 'dt_seconds': DT, 'tau_seconds': TAU,
                       'epsilon_C_per_second': EPS, 'k_numeric': K, 'threshold': float(MC)},
        'offset_checks': offsets,
        'quadratic_family': {
            'mu2_seconds_squared': float(mu2), 'interior_duration_seconds': float(length),
            'gate_curvature_bound_C_per_second_squared': float(gate_bound),
            'score_curvature_bound_C_per_second_squared': score_bound,
            'sufficient_joint_bound_C_per_second_squared': float(joint_bound),
            'condition': 'b > max(epsilon/(2*t_p), tau*M_c/(mu2*L)); p>=1; L>0; t_i=i*dt',
            'cases': family,
        },
        'status': 'PASS',
    }
    print(json.dumps(output, indent=2))

if __name__ == '__main__':
    main()
