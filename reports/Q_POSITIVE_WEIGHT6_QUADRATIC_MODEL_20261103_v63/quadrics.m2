R=QQ[z0,z1,z2,z3,z4,z5,z6,z7];
I=ideal(z1^2 + 120*z5*z7 - 14*z7^2,(15*z1*z3 - 5*z1*z6 + 12*z3*z7 + 8*z4*z5 - 2*z4*z7 + 16*z5*z6 - 8*z6*z7)/15,(3*z1*z4 + 2*z1*z6 + 36*z3*z7 - 16*z4*z5 + 12*z4*z7 - 32*z5*z6 + 8*z6*z7)/3,(4*z1*z5 + 36*z5*z7 - z7^2)/4,z7*(z1 - 20*z5 + 5*z7));
G=gb I;
print("DIM="|toString dim(R/I));print("DEG="|toString degree(R/I));print("HS="|toString hilbertSeries(R/I));
K=ideal(z1,z7,z5*(z4+2*z6));assert(gens I%gb K==0);
assert(z1^3%G==0);assert(z7^3%G==0);assert((z5*(z4+2*z6))^2%G==0);
assert(dim(R/I)==5);assert(degree(R/I)==3);print("EXACT_RADICAL_TWO_FIVE_PLANES_PASS");
