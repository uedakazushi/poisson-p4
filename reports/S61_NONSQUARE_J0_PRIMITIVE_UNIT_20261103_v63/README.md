# New global J0 four-constant unit pivot

For the pinned23x57 residual over
QQ[a,e,1/(a*F*Bs*C*(3e+2))], the exact row7 identity is

    M[7,4] - M[7,15] - M[7,25] + M[7,30] = c*a*C*(3e+2)^3,

where c is the explicit nonzero rational in report.json. All indices are
zero-based. The four source polynomials and all10 terms of their difference
are retained in identity.json.gz. Independent QQ coefficient arithmetic
verifies the entire identity. No unknown factor or J was inverted.

The constant column replacement is unimodular, and an exact native pivot
replay gives22x56. This is one new safe pivot,35total from the original57
basis. It does not prove fullrank57, fullcokerunit, defectfiniteness or
ambient nilpotence. The new39.5MB residual is reproducible in about5seconds;
its exacthash is retained rather than duplicating the regenerated file.

Replay the compact identity in any checkout:

```sh
python code/s61_nonsquare_J0_primitive_identity_replay_20261103_v63.py
```

The source-matrix binding uses the pinned23x57 hash from the previous J0
checkpoint. The saved primitive-search.cpp and row7.witness reproduce the
original-row witness; polynomial_unit_span_pivots_20261103_v62.cpp applies
the saved witness with metadata=1,colclean,target-row7. See the note for
the unchanged known-unit set and the limitations of the earlier refusals.

Proof: research/notes/S61_NONSQUARE_J0_PRIMITIVE_UNIT_20261103_v63.md.
