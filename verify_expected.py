#!/usr/bin/env python3
"""Compare self-contained numerical checks to the research supplement outputs."""
from __future__ import annotations
import json
from pathlib import Path
import subprocess
import sys
import math

ROOT = Path(__file__).resolve().parent
TOL = 1e-10

def compare(expected, actual, label='root'):
    if isinstance(expected, dict):
        assert isinstance(actual, dict) and set(actual)==set(expected), f'{label}: key mismatch'
        for k in expected: compare(expected[k], actual[k], f'{label}.{k}')
    elif isinstance(expected, list):
        assert isinstance(actual, list) and len(expected)==len(actual), f'{label}: list mismatch'
        for i,(e,a) in enumerate(zip(expected,actual)): compare(e,a,f'{label}[{i}]')
    elif isinstance(expected, bool):
        assert actual is expected, f'{label}: boolean mismatch'
    elif isinstance(expected, (float,int)) and not isinstance(expected,bool):
        assert isinstance(actual,(float,int)) and math.isclose(expected,actual,rel_tol=0,abs_tol=TOL), f'{label}: expected {expected}, got {actual}'
    else:
        assert expected==actual, f'{label}: expected {expected}, got {actual}'

def main():
    specs = [
      ('verify_controls.py','checks/constructed_controls_results.json'),
      ('verify_additional_checks.py','checks/additional_checks_results.json'),
    ]
    for script,expected_file in specs:
        output = subprocess.check_output([sys.executable,str(ROOT/script)],cwd=ROOT,text=True)
        expected=json.loads((ROOT/expected_file).read_text(encoding='utf-8'))
        actual=json.loads(output)
        compare(expected,actual)
        print(f'PASS: {script} matches {expected_file} within {TOL:g}')
    print('PASS: all public mathematical reproducibility checks')

if __name__=='__main__':
    main()
