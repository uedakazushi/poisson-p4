# Five graph cubic coefficients forced zero by the actual cover

The fixed known initial module cover of the original comparison C3 is
`M=A eA ⊕ B eB ⊕ H eH`, where
`A=QQ[r,ell]/(gA3,ell³)`, `B=QQ[r,ell]/(gB3,ell³)` and
`H=QQ[r,ell]/(r,ell)²`. This is a cover of the associated graded module;
it does not split the actual comparison. Its Hilbert function is
36642 and its length21. The first21 and complete quadratic source
proofs fix the raw eA/eB sectors. The quadratic H relations do not kill
eA r ell² or its r-multiple eA r²ell².

The checker constructs the full degree3 NF map to
`(A r²ell,A r ell²,B r²ell,B r ell²)` and degree4 NF map to
`(A r²ell²,B r²ell²)`. It verifies the actual source cubics gA3/gB3,
all ell³ multiples, and the free unit columns. In particular each
r ell² vector has literal NF(0,1) in its block, and each r²ell²
vector has literal NF1, independent of any unknown coefficient.

For each physical row63,64,68,73,74, the accepted tauA constant kernel
contains the same free unit source e_j: bB_j,2 has zero plane part.
The full tauA construction therefore gives A cubic delta_j3 while
retaining arbitrary B cubic carry. The tauB column has zero first.
Its delta2 eB terms are only zz; its actual mixed bA2 H terms are
retained and cleared with r/ell times RH. These clears have zero B
plane cubic contribution because RH2_B,plane=delta_i2,plane=0. Thus
the tauB source gives B cubic delta_j3 and no A cubic component.

Existing C2 readbacks give delta_j3 r³=r²ell=0. Modulo ell³ write
delta_j3=d r ell²; the ell³ coefficient remains unknown. The two
actual source cubic NF vectors are then `(0,d,u,v)` and `(0,0,0,d)`.
Their fixed A2/B2 minor is d². Their r-multiples have quartic NF
`(d,-(gB_r²ell/gB_r³)u+v)` and `(0,d)`, again with determinant d².
Every beta2 and mixed source correction is preserved in u,v.

If d were nonzero, two independent degree3 and two independent degree4
initial relations would remove at least four dimensions from M. This
would give length(C3)≤17. The independently accepted actual DVR rank4,
length(C/ell)=8 and length(C/ell²)=14 give torsion exponents
`(1,1,e3≥2,e4≥2)`, hence length(C3)≥18. Therefore d=0 for all five
rows. No graph-root readback, point, implicit solve or new jet is used.

This excludes exactly the five graph cubic classes. It does not prove
that every cubic source vanishes, determine the actual C3 length, or
settle original H0 length or the full nil cube.
