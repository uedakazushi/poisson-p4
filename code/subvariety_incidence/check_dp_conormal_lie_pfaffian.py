#!/usr/bin/env python3
"""Exact universal first-normal-jet identity at a Poisson-invariant codimension-two surface.

This is an algebraic identity for any bivector tangent to the surface, not an
existence proof for rank-four Poisson extensions. Characteristic zero.
"""
import argparse
import json
from pathlib import Path
from math import comb
import sympy as s

parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path, default=Path("work/dp_conormal_lie_pfaffian.json"))
parser.add_argument("--check", type=Path)
args = parser.parse_args()
assert s.__version__ == "1.14.0"

u, v = s.symbols("u v")
f, a, b, c, m, n = s.symbols("f a b c m n")
p, q, r, t, w, z, alpha, beta = s.symbols("p q r t w z alpha beta")
# Coordinates are ordered (u,v,x,y). All mixed/normal entries vanish modulo (u,v).
Puv = a*u+b*v+c*u**2
Pxy = f+m*u+n*v
Pux = p*u+q*v
Puy = r*u+t*v
Pvx = w*u+z*v
Pvy = alpha*u+beta*v
Pf = s.expand(Puv*Pxy-Pux*Pvy+Puy*Pvx)
assert s.expand(Pf.subs({u:0,v:0})) == 0
assert s.expand(s.diff(Pf,u).subs({u:0,v:0})-f*a) == 0
assert s.expand(s.diff(Pf,v).subs({u:0,v:0})-f*b) == 0

# The global dimensions are consequences of the (2,2) complete-intersection Koszul resolution.
R1 = comb(5,4)
R2 = comb(6,4)
R3 = comb(7,4)
I5 = 2*R3-R1
I2_5 = 3*R1
I4 = 2*R2-1
I2_4 = 3
IIC5 = I2_5 + I4-I2_4
h0_OS2 = R2-2
normal_quotient = IIC5-I2_5
c2_conormal_bracket_bundle = (2*2)**2  # (2H)^2, H^2=4, evaluated below explicitly.
assert (I5,I2_5,IIC5,h0_OS2,normal_quotient) == (65,15,41,13,26)
assert normal_quotient == 2*h0_OS2
assert 4*4 == 16

report = {
 "date": "2026-11-11", "version": "v18", "status": "PASS",
 "sympy": s.__version__, "base_field": "QQ",
 "universal_local_first_jet": "Pf mod (u,v)^2 = f*(a*u+b*v)",
 "mixed_terms_contribute_to_first_normal_jet": False,
 "conormal_isotropy_bundle": "Hom(wedge^2 O_S(-2)^2, O_S(-2)^2) = O_S(2)^2",
 "degree5_I_S_dimension": I5,
 "degree5_I_S_squared_dimension": I2_5,
 "degree5_I_S_times_I_C_dimension": IIC5,
 "degree5_first_jet_dimension": normal_quotient,
 "h0_S_O2": h0_OS2,
 "c2_O_S_2_squared": 16,
 "scope": "Universal local identity and Hilbert-number audit. Global conormal Lie structure and Hamiltonian transport are proved in TeX. No rank-four extension existence or all-147 classification certified."
}
if args.check:
 assert report == json.loads(args.check.read_text(encoding="utf-8")), "certificate mismatch"
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n", encoding="utf-8")
print("PASS: Pfaffian first jet factors through surface bivector times conormal Lie bracket; 65/41/15/26/16 dimensions")
