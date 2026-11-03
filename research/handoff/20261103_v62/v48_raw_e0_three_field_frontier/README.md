# One raw-coordinate lift at the same three fields: frozen partial frontier

Extract after v36, v45 and v47, including v47/CORRECTION.md. No new prime values
are used. The same redSB/liftstd protocol computes ONLY h*e0, h=a^5J³, over the
existing fields2147483647 and2147483629. Their complete160-column canonical
bases are byte-for-byte identical to the previously saved slimgb bases in those
same fields. No existing canonical basis is overwritten.

The first large-field wrapper initially used `size` of a matrix as a zero test.
Singular counts matrix entries even for a zero matrix. The original25.003-second
wrapper/stdout is retained as INVALID. The subsequent check correctly casts the
matrix difference to a module. That run takes27.932seconds, so the first field's
cumulative solver time is52.935seconds, below its60-second budget. The second
field takes24.545seconds. Both use at most441MiB and no workers. A fresh2-by1
zero-matrix fixture records size(Z)=2 and size(module(Z))=0. Both old v47 wrappers
are also explicitly marked INVALID; their cofactor data were validated only by
independent direct polynomial identities.

For EACH of the three fields, separate FLINT multiplication against the literal
QQ raw source reduced at that prime verifies ALL24 equations M*L=a^5J³*e0.
Exact common-content division strips a^4J³, and another ALL24 direct check gives
M*Lprimitive=a*e0. The primitive cofactor has maximum degree11,8919 terms in both
large fields and8918 at32003. Their full supports agree after explicitly padding
ONE vanished32003 coefficient. These are complete raw-matrix field certificates,
not a full24-target inverse and not a characteristic-zero result.

The bounded three-field77-bit CRT has ZERO equal balanced rational reconstructions
between adjacent prefixes on8919 positions. Farey reconstructions cover5427 and
5365 positions, which is only the random coverage baseline. There is no complete
QQ candidate or QQ polynomial identity. Matching protocol/support does not itself
prove the lifts are reductions of a common rational cofactor. No QQ support solve,
fourth field, degree5 coefficient test or further CRT extension is performed.

All field cofactors, primitive data, canonical bases, source pins, corrected and
INVALID wrappers/receipts, independent verification sources/logs, support/CRT
scripts and compact CRT state are saved. No ELF is included. The already proved
47 QQ seed identities and degree-at-most4 obstructions are separate results;
they do not turn this partial cofactor reconstruction into a QQ rank41 proof.
