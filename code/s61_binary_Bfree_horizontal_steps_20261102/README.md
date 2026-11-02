# B-free horizontal certificate steps

These small deterministic steps reconstruct the three Schur obstructions from the original 135 equations. The public producer executes them in a private directory. Literal `/tmp/` destinations are remapped to that directory; `work/` paths refer only to outputs generated earlier in the same run. No pre-existing ignored work file or large cache is an input. The mathematical expressions and assertions are retained from the separately checked derivation. Stage scope strings describe intermediate outputs, not the final theorem.

The common coordinate convention is `old_light = D Ti(new_light), old_heavy = H`, equivalently `D^w f(Ti(new_light), H/D^2)` in weighted degree w. All row contents used in this directory are nonzero elements of QQ(a,e); no conclusion at their zeros is made.
