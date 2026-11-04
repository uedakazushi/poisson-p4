# Ten actual pure cubic source kernels on the whole P100 curve

Date: 2026-11-04. Exact primary certificate; independent review pending.
The source-only preflight gives six pivot sources, the original-frame sigma
lifts and ten formal source-kernel vectors. This job evaluates their actual
pure84 coefficients on the entire P100 unit ring.

The pivot genuine indices are `[1,6,8,9,11,14]`; their504 accepted values
are reused. The nonpivots are `[0,2,3,4,5,7,10,12,13,15]`, actual135 labels
`[18,23,39,40,43,68,92,115,120,134]`. G18's84 values are also reused.
All other nine complete normalized source rows are streamed from the
immutable full57-by91 matrix. Each core coefficient is discarded after
specialization; the139746-term source JSON is never loaded. All112 full
H-linear coefficients are checked against this exact normalized matrix.

The saved native Tn normalization already divides by the accepted row unit.
It is not applied again. The six cached pivot bindings record their actual
135 labels, basis positions, coefficient-root proof and common pure84 order;
the504 coefficients are not reevaluated. No restricted-ZW source is used.

For each nonpivot i the job constructs

```
alpha_i,p = sum_j H_i[j]*sigma_j,p,
K_i = C_i-sum_p alpha_i,p*C_p.
```

It verifies all16 H_i.eta identities, all49 sigma lift coefficients and all36
pivot projection identities. Both the original and actual delta7 source
cofactor frames are saved. All140 H-kernel identities and all100 entries of
the ten-by-ten nonpivot identity block pass. Every pure coefficient is checked
in the reversed original-frame sum and again using the actual delta7 frame:
1680 literal coefficient identities in total. Denominators use only the five
known P100 curve units.

The pivot-first source permutation has sign+1; the lower block triangular
source change has determinant1. The actual delta7-frame determinant is
1/delta, where delta=-64000*C^2 is a whole-curve unit. Hence these ten vectors
are a basis of the formal source H-kernel on the whole curve. Their actual
pure cubic images can still have relations modulo the quadratic ideal;
neither image rank10 nor a later70-row quartic rank is asserted.

The ten pure84 rows have RF term counts
`[1258,1635,1258,1631,1720,1104,1972,1804,1888,1972]`.
Seven have84 nonzero coefficients; the three literal pure cubics have77.
Each output includes its full16 source cofactors in both frames, so its
provenance remains directly reusable for actual polynomial products.

## Bound and evidence

The proved live bound was58998 RF terms; the observed peak was32359,
including a16797-term stream/literal reserve. The input streams79559 core
terms, but only one core coefficient remains live during its evaluation.
The job serialized16766 terms and took5.981 seconds including compilation;
native arithmetic took0.413 seconds. No cap was raised or job repeated.

Code: `code/s61_nonsquare_P100_pure_cubic_kernel_20261104_v63.py` and `.cpp`.
Reports: `S61_NONSQUARE_P100_PURE_CUBIC_KERNEL_20261104_v63`, with exact
`source-binding.json`, live `preflight.json`, `audit.json`, and ten
`pure-cubic84-i.json` files. Each saves the old and actual delta7 cofactor
vectors together with its84 pure values.

No70 quartic products, NF117 evaluation,15 mixed Koszul/four overlap vector
construction, matrix rank, new oracle value, point, prime or GB is included.
The original quartic rank bounds are unchanged by this coefficient receipt.
All primary bytes are immutable.
