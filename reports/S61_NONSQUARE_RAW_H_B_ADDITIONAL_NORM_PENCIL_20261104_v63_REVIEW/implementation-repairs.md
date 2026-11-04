The frozen successful report pins both preserved failure directories.

The first new-checker attempt failed at compilation because the helper
name `plus` conflicted with `std::plus`. No native determinant or GCD
was run in that attempt. Renaming the helper fixed only compilation.

The next native attempt passed the additional determinant and full
sixteen-entry adjugate checks but failed at its first Laurent source
membership comparison. The checker's manual a-exponent shift had
sorted FLINT terms without reducing the rational coefficient content.
This was a checker canonicalization defect. A separate reviewer-owned
Fraction dictionary multiplication verified the entire source identity:
two 4,862-term products cancel to the saved 575-term norm, with zero
difference. The repair canonicalizes each reconstructed shift and
adds FLINT canonicality guards throughout the checker.

The final complete fixed replay passes 309 coefficient identities,
nine exact GCD checks and 33 last-row Laplace determinants in 0.671
native seconds (2.210 seconds for the final compile/driver receipt).
No author arithmetic identity changed, no PASS was claimed before
that final replay, and no GB, resultant, parameter point or prime
was used in the repairs. All failed source bytes and logs remain
under `failed-compile/` and `failed-native/`.
