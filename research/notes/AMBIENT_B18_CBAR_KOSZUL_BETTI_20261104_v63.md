# Actual Cbar: finite Koszul Betti numbers

This packet uses only the independently accepted actual eight coordinate actions on
Cbar=C/ell C. The coefficient ring is the regular eight-variable completed CA-F9
ring, with coordinates zA_0,...,zA_6,r_CA. The common r is r_CA, not r_CB.

The sparse Koszul maps have dimensions 8,64,224,448,560,448,224,64,8 and a total of
1024 nonzero entries. Fixed rational column elimination records every operation,
normalized pivot and zero remainder. Direct sparse compositions verify d²=0 in
every column. The differential ranks are 5,35,105,175,175,105,35,5. Therefore the
minimal Betti numbers are

    3,24,84,168,210,168,84,24,3.

Indeed tensoring the Koszul resolution of the residue field with the actual Cbar
computes Tor over this regular coefficient ring; these dimensions are the ranks
of its minimal free resolution. This does not identify a resolution of the
comparison kernel K or the original D0.

The same eight actions give m-adic filtration dimensions 8,5,2,0, Loewy Hilbert
function (3,3,2), minimal generator count 3, and socle dimension 3. Cbar is treated
as a module; no unital algebra structure is asserted. No generic-DVR rank,
comparison-kernel torsion, or D0 depth follows from this finite calculation.

All source values are reused unchanged. No source solve, Hensel expansion,
Groebner basis, field point or finite prime is used. The report pins the source,
prior independent review and all sparse elimination receipts.
