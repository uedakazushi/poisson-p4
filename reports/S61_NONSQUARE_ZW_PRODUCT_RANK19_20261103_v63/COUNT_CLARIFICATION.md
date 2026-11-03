# Cofactor support count clarification

The frozen producer reports 117 terms by summing SymPy `Poly.terms()` lengths.
One of the five cofactor polynomials is zero, whose SymPy representation contributes
one serialized zero term. The native witness support lengths are
`[38,2,38,0,38]`, so there are **116 nonzero cofactor terms**.

The frozen source, JSON certificate, native witness, literal QQ sum=1, actual
unit pivots and all rank/Hilbert conclusions are unchanged. The archive hash
continues to identify the original frozen bytes. Independent review identified
this bookkeeping convention; it has no algebraic effect.
