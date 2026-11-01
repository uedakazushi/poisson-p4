# Generic exclusion of all nonnilpotent modular infinity

Date: 2026-11-01. Exact rational rank stratification, conditional on the archived exhaustive fixed-trace normal forms.

There is a nonempty Zariski open subset of W16 over which no solution of the 134 leading equations has a nonnilpotent modular matrix. The same holds for the 135-leading locus. The exceptional beta image has dimension at most fifteen. This is a generic-base assertion, not an assertion on the whole discriminant complement or on the full affine S61 fiber.

## 1. Complete point coverage

The argument uses point normal forms, not only containment in component closures. The precise historical dependencies are the v03 regular-semisimple theorem, the v04 semisimple-boundary classification, and the v05 mixed-modular theorem. Together they cover every nonnilpotent Jordan type as follows.

1. For regular semisimple modular matrix, every point is symplectically and projectively equivalent to a member of G1 or G2, or one of T1,...,T4. The full resonance analysis, including zero values on the resonant lines, is part of v03. Zero parameters on the exceptional lines are in the log families.
2. For nonregular semisimple trace xy, the reduced fixed fiber has two affine-three planes. Its SL2 action is conjugation on the traceless matrix T=[[-c-d/2,a],[b,c+d/2]], d=0 or1/2. Semisimple T, including zero, gives a log normal form. Nonzero nilpotent T is one orbit for each d; the report includes both representatives a=1,b=0,c=-d/2.
3. For nonregular semisimple trace xy+zw, the v04 chart with U=Qxz semisimple has thirteen explicit curve branches. The three baseline lines are log families. On each of the ten other branches, every nonzero curve parameter is conjugate under the residual diagonal torus to parameter one; all ten representatives, including Kplus and Kminus, are included separately in the report. Parameter zero lies on a baseline. The dual chart with V=Qyw semisimple is obtained by a symplectic transformation and an irrelevant nonzero scalar. The remaining U=V=0 locus is one orbit, represented by Qxy=xy,Qzw=zw, a log point. Nonzero nilpotent U or V is impossible by the archived unit-ideal certificate. Thus no boundary point is omitted.
4. A nonnilpotent, nonsemisimple matrix has one of the two mixed Jordan forms in v05. The equal-nonzero-eigenvalue form has empty fixed fiber. The other form, trace xy+z²/2, has exactly three complete affine lines, all retained below.

A member of G1 with first trace coefficient zero is covered by exchanging Darboux blocks, unless its modular matrix is zero. The missing G2 chart A=0,B!=0 is included separately as N(0,1/2,0,1/2,e). It is precisely the semisimple-T part of the second one-zero trace plane after exchanging Darboux blocks; it cannot be discarded by referring back to the A=1 chart. Thus the three explicit log charts cover every log point with nonzero modular matrix.

## 2. The two log families, including every rank drop

Use

N(A,B,C,D,E;P,R):
Qxy=Axy+Bzw, Qxz=Exz, Qxw=(C-E)xw+Py²,
Qyz=(-B-E)yz+Rx², Qyw=(B-C+E)yw, Qzw=Cxy+Dzw.

The two normalized charts are

G1=N(1,0,0,d,e), G2=N(1,d/2,0,d/2,e).

Their beta-linear leading matrices split into twelve disjoint row/column blocks. Every block has generic rank one, with all 2 by 2 minors identically zero. The gcd of all nonzero entries of each block is a nonzero constant. Each entry factors into affine-linear polynomials in d,e. Therefore each block can vanish only at intersections of two distinct factor lines, and enumeration of those intersections is exhaustive over C. The code evaluates the entries at every candidate and keeps exactly the common zeros.

The union of exceptional parameter points consists of 27 rational points for G1 and 23 for G2. Outside those finite lists the total rank is twelve. At every point on the lists it is at least ten; the entire lists and ranks are saved, not inferred from a grid sample.

The diagonal two-dimensional Cartan torus fixes every log tensor, including exceptional parameters. Its Sp4 orbit therefore has dimension at most eight. The incidence over the open parameter plane has dimension at most

2+8+(16-12)=14.

At any of the finitely many exceptional parameter points the bound is

8+(16-10)=14.

Thus those two log-family beta images have proper closed closure.

The additional boundary chart N(0,1/2,0,1/2,e) has twelve generic rank-one blocks and two identically zero singleton columns. The entry gcds show that its only rank drops are e=-1,-1/2,0,1/2; the rank is eleven at each and twelve elsewhere. The same Cartan torus fixes every point. Its incidence bound is therefore 1+8+4=13 off the exceptional values and 8+5=13 at those values. The corrected certificate includes the complete chart and all four exceptional points.

## 3. Finite exceptional orbits

The report computes exact rational minors for sixteen representatives: T1,...,T4, the two nonzero nilpotent SL2 orbits in the one-zero trace case, and all ten nonzero exceptional branches of the equal-eigenvalue chart. Every rank is at least eleven. Even the coarse map from Sp4 times the beta kernel has dimension at most

10+(16-11)=15.

Scalar multiplication of Q is irrelevant both projectively and for the zero set of the quadratic leading forms. No additional scalar parameter is needed in this incidence bound.

## 4. Complete mixed-Jordan lines

The three lines have brackets

Qxy=(1-b)xy+(1/2-c)z², Qxz=(b/2)xz,
Qxw=(1-a-2c)xz+(b/2)xw, Qyz=-(b/2)yz,
Qyw=a yz-(b/2)yw, Qzw=bxy+cz²,

with a arbitrary and (b,c) respectively (0,1/4),(0,1/2),(1/2,1/2). The exact beta minors are constants in a, of ranks eleven, eleven and twelve. There are no unexamined exceptional parameter values.

The nonzero modular field belongs to sp4 and preserves Q; the code checks both identities universally in a. Hence each Q has a stabilizer of dimension at least one, and its Sp4 orbit has dimension at most nine. The whole incidence therefore has dimension at most

1+9+(16-11)=15.

This covers the complete three lines, including values where other invariants degenerate.

## 5. Conclusion and scope

The finite union of the above constructible beta images has closure of dimension at most fifteen. Removing that closure gives the asserted nonempty base open. Properness of the incidence is not needed: the dimension of the closure of a constructible image equals its dimension. This argument does not use the invalid principle that a single finite special affine fiber controls every generic component.

Combining with the separately certified Jordan-(2,2), regular-nilpotent, and minimal-open exclusions leaves only M0,M1 among the nonzero-modular types over sufficiently general beta. Both actually survive, and their modular image lines lie on E_tau. Their exceptional fibers still require the separate analysis described in the minimal-nilpotent notes.

## Reproduction

```
.venv/bin/python code/s61_nonnilpotent_infinity_exclusion_20261101.py --output work/nonnilpotent_check.json
```

Frozen report: `reports/S61_NONNILPOTENT_INFINITY_EXCLUSION_20261101.json`. It records all block matrices and factor lines, all 50 two-parameter exceptional points and the four boundary-chart exceptions, all constant representatives and minors, the universal mixed-line minors, and hashes of the three historical point-exhaustion sources. No random sampling, finite-field inference, or generic-rank-only exclusion is used.

Review correction: an independent audit identified the missing A=0 chart in the first draft. The current code, report and argument explicitly include it; the dimension bound and theorem remain unchanged.
