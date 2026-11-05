# Five active variables suffice for the constant d3 source rank

The actual C2 actions of zA2,zA3,zA6 vanish identically. Splitting the
periodic/Koszul complex by the fixed wedge subset of these three zero
variables leaves the same complex in the five active variables
zA0,zA1,zA4,zA5,rA. Put active wedge factors first and zero factors last;
the stored signed permutation accounts for every exterior orientation.
All304 full d2 and1137 full d3 coefficients are reconstructed exactly.

The active d1,d2,d3 shapes are14-by84,84-by224,224-by364, with37,193,447
nonzero entries. Full d1 rank11 and full d2 rank89 imply active d2 rank
56. Thus the active d3 rank is at most224-56=168. Full d3 rank is its
rank plus3*56+3*11, and the actual second C2 Betti number is228 minus
the active rank.

The one proposed rank job uses at most37632 independent-basis slots,
61152 factorization slots,14028 pivot-operation slots,447 input slots
and447 sparse-minor slots, plus10000 metadata slots. It saves a sparse
original selected minor and an exact factored nonzero determinant;
it does not duplicate a dense168-square inverse. Hard384bit coefficient,
150000 slot,32MiB,60s and2GB guards apply. No physical source, jet,
homotopy coefficient, full518-by1302 rank or T2 length is computed.
