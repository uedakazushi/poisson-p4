31 SEMISIMPLE NODE EXCLUSION
See31_ss_exclusion_proof.tex. Geometric prerequisite: both residue slopes nonzero,
TCB(p)=0; eigenvalue/vanishing-order arguments reduce to7 slices below.
This directory handles algebraic exclusion ofall7; geometric exhaustiveness
isprovedseparatelybyroot/junction, notinferredfrom coefficient sampling.

CERTIFICATES TO INCLUDE
31_D_SS_m1_certificate.json
31_D_SS_m3_certificate.json
31_D_SS_3_certificate.json
31_D_SS_m3_2_certificate.json
31_D_SS_3_2_certificate.json
31_E2_SS_m3_2_certificate.json
31_E3_SS_m2_certificate.json
31_D_SS_m1_final_family.json

Firstcertificate isNOTunit: fullremainingideal=(z9*z21-1/9).
Itparameterizes anactualJacobi/E8-first-jet family withquinticF quadratic*cubic.
Normality fails alongsurfaceintersectionofbothfactors.
Allother6certificates areunit. No parameterdivisionusedincertificatechains.

STDLIB AUDIT
python3 audit_31_ss_stdlib.py --model PATH_TO_V21_31_JSON --jets PATH_TO_V22_31_E8_JSON
Dependencies: standardPythononly; imports audit_certificates_stdlib.py.
Report31_ss_audit_stdlib.json: all1243checksPASS.
Audit reconstructsP->H/B, all186necessaryE8eqs, allcertificatechains,
remaininghyperbolaideal,andFfactorization. DirectJacobi+Pfaffiangradient
verification includedforactualremainingfamily.

Old31_A/31_B/31_Cgeneric outputs are exploratory; doNOTinclude/needthem.
The nonzero-TCB(p) branch isexcludedbygeometriccenterdegree+flowinvariance.
Mixed-residuebranches areexcludedbygeometriczero-slope-nilpotencylemma.
NN ishandledbychainsagent, separatecertificatepacket.
