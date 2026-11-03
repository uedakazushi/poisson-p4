The archived `replay-fixed.sing` repairs the target dimensions but its check is
still invalid: Singular's `size` of a zero matrix counts its entries. A fresh
zero-matrix fixture gives `size(Z)=2` and `size(module(Z))=0` for a 2-by-1 zero
matrix. Both archived wrapper checks must therefore be treated as INVALID.

The original and primitive cofactor certificates remain valid: separate FLINT
multiplication verified all twenty-four literal polynomial identities. Their
acceptance does not use either invalid wrapper check. The subsequent v48 source
uses `size(module(Q*L-RHS))` and retains both earlier failed-wrapper receipts.
No archived data or hash has been replaced.
