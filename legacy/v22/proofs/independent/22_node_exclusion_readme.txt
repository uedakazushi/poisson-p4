CERTIFIED RESULT (subject to integrating root's E8 local-normal-form lemma):
The 22_node configuration cannot be the reduced singular support of a normal
actual quintic degeneracy divisor of a generically rank4 Poisson structure onP4.

PROOF
22_node_exclusion_proof.tex gives the complete reduction. The root's pure conic
jet lemma excludes z30=0 orz14=0; this has been independently audited here.
The open branch z30*z14 nonzero is excluded by four exact certificates:
1. 22_node_NN_heis_both_slopes_certificate.json: slopes equal; λ denominators.
2. 22_node_SS_heis_both_SS1_certificate.json: common nonzero slope; unit; λ only.
3. 22_node_NN_heis_both_certificate.json: both nilpotent; unit; λ andλ+1.
4. 22_node_NN_spectral_exception_certificate.json: λ=-1; E8 spectral unit; QQ only.

MODEL DERIVATION
heisenberg.py and spectral.py build conormal bracket/tangent-action matrices
from the archived full150 coefficient model. H convention: row=input covector.
The second conic normal frame is(q2,x1,x2), whereas root generator uses(x1,x2,q2).
With cyclic permutation (1,2,0) and transposeH, all36 matrix comparisons PASS.

CAS AUDIT
PYTHONPATH=/tmp/p4_v6_math python3 audit_certificates.py
reads the original v21 Jacobi model and exact jet-coefficient data, reconstructs
each stage from the preceding equations and substitutions, checks all rational
linear-combination identities, checks substitution equations, and records all
irreducible denominator factors. All176 checks passed.
certificate_audit.json and independent_jet_comparison.json are finalreports.

MATERIAL TO INCLUDE
The four named certificates,22_node_heisenberg.json,22_node_spectral.json,
prooftex,readme,audit script,and auditreports.
Other *_reduced.json files and reduce_* scripts are exploratory,unneeded.
31_external Heisenberg/spectral data were generated but NOT used for a complete
TC+line exclusion. No claim about31 follows from this subtask.

INDEPENDENT STDLIB AUDIT
python3 audit_certificates_stdlib.py --model PATH_TO_V21_22_MODEL
No SymPy or other CAS; only ast/Fraction sparse polynomial arithmetic.
Reconstructs P->B,H, all necessary jet-equation coefficients, Hessian identities,
all7 handwrittenrelations and2 derivationidentities, all4 certificatechains.
All228 checks PASS. See certificate_audit_stdlib.json.
Denominators handled exactly inQ[lambda,lambda^-1,(lambda+1)^-1]; no sampling.

ADDITIONAL211 AUDIT
audit_211_certificate_stdlib.py imports the same stdlib arithmetic framework
and audits root211_triple certificate, including independent reconstruction
ofall three components and allE8 equations. Geometricexhaustiveness ofits
initial substitution remains separately proved byroot/junction agent.
