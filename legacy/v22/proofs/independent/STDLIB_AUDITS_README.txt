P4 V22 INDEPENDENT EXACT AUDIT PACKET

All verifiers below need only Python's standard library. No SymPy, CAS, numerical
sampling, or network. Keep these five files in the same code directory:
  audit_certificates_stdlib.py
  audit_211_certificate_stdlib.py
  audit_31_ss_stdlib.py
  univariate_fraction_stdlib.py
  audit_31_nn_stdlib.py

22 conic-pair:
python3 audit_certificates_stdlib.py --model V21_TREE_JACOBI/22_node.json \
  --data-directory TWOCOMPONENT_CERTIFICATES --output REPORT.json
228 checks. Data: four22certificates namedin22_node_exclusion_readme.txt,
22_node_heisenberg.json,22_node_spectral.json. Reconstructs matricesfromP.

211 triple:
python3 audit_211_certificate_stdlib.py --model V21_TREE_JACOBI/211_triple.json \
  --jets V22_E8_JETS/211_triple.json --certificate 211_triple_certificate.json \
  --output REPORT.json
110 checks. Reconstructs allthreecomponentjet matricesandallnecessaryE8equations.

31 semisimple cases:
python3 audit_31_ss_stdlib.py --model V21_TREE_JACOBI/31_external.json \
  --jets V22_E8_JETS/31_external.json --certificate-directory SS_CERTIFICATES \
  --output REPORT.json
1243 checks. Seven31_D/E2/E3_SS_*_certificate.json files aslisted
in31_ss_readme.txt. Includes fullP->allE8jets reconstruction. Sixunitcertificates;
one exactremaining hyperbola z9*z21=1/9 withquadratic*cubic actualquintic.

31 nilpotent cases:
python3 audit_31_nn_stdlib.py --model V21_TREE_JACOBI/31_external.json \
  --jets V22_E8_JETS/31_external.json --certificate-directory NN_CERTIFICATES \
  --output REPORT.json
658 checks. Filesfromchains:
  31_C_NN_k3_certificate.json
  31_NN_k3_final_certificate.json
  31_NN_k1_II_certificate.json
  31_NN_k1_IA_certificate.json
  31_NN_k1_IB_certificate.json
  31_NN_k1_IA_rational_certificate.json
  31_NN_k1_IB_quadratic_field_certificate.json
  31_NN_k1_IB_rational_certificate.json
  31_NN_partition_certificate.json
Generic coefficient fields Q(z17) orQ(z21) checkedexactly byunivariate
Euclidean polynomial arithmetic. Everyextra denominatorzero exceptionhas
its owncertificate. Quadraticexception ischeckedinQ[theta]/(theta^2+13511),
coveringbothembeddings. Rationalfamily hasexactgraphrelations
z11=3*z18*z20,z1=3*z20^2; itsactualF liesin(x2,x3)^2 andisnonzero.

GEOMETRIC DEPENDENCIES
These verifiers certify allalgebraic identities, sourceprovenance, eliminations,
denominatorcoverage, andnonnormality offinalfamilies. TheydoNOTreplaceproofs
ofgeometricexhaustiveness ofinitialsubstitutions. Thoseareseparatelydocumented
inroot/junction/chainsmanuscripts. Square-zero conclusions areclosed-point
arguments, neveridentificationofanidealwithitsradical.
