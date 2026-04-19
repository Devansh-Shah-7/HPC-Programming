#!/usr/bin/env python3
"""
Compares two mesh output files with a floating-point tolerance.
Usage:  python3 verify.py Mesh.out Test_Mesh.out
"""
import sys, math

def load(path):
    vals = []
    with open(path) as f:
        for line in f:
            vals.extend(float(v) for v in line.split())
    return vals

if len(sys.argv) != 3:
    print("Usage: python3 verify.py <generated> <reference>")
    sys.exit(1)

gen = load(sys.argv[1])
ref = load(sys.argv[2])

if len(gen) != len(ref):
    print(f"FAIL: size mismatch  generated={len(gen)}  reference={len(ref)}")
    sys.exit(1)

tol   = 1e-5
worst = 0.0
fails = 0
for i, (g, r) in enumerate(zip(gen, ref)):
    err = abs(g - r)
    if err > tol:
        fails += 1
        if fails <= 5:
            print(f"  mismatch at index {i}: got {g:.8f}  expected {r:.8f}  diff {err:.2e}")
    worst = max(worst, err)

if fails == 0:
    print(f"PASS  (max error = {worst:.2e}, all {len(gen)} values within tol={tol})")
else:
    print(f"FAIL  {fails}/{len(gen)} values exceed tol={tol}  max_err={worst:.2e}")
    sys.exit(1)
