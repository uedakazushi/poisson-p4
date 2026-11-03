R=QQ[x0,x1,x2,x3,x4];
F=-8*x0^5 - 720*x0^3*x2*x4 + 360*x0^3*x3^2 + 540*x0^2*x1^2*x4 - 2160*x0^2*x1*x2*x3 + 270*x0^2*x2^3 - 2700*x0^2*x4^3 + 3240*x0*x1^3*x3 + 810*x0*x1^2*x2^2 - 16200*x0*x1*x3*x4^2 - 12150*x0*x2^2*x4^2 - 2430*x1^4*x2 + 24300*x1^2*x2*x4^2 + 24300*x1^2*x3^2*x4 - 24300*x1*x2^2*x3*x4 + 6075*x2^4*x4 - 60750*x2*x4^4 - 121500*x3^2*x4^3;
J=ideal(diff(vars R,F));G=gb J;
print("JAC_DIM="|toString dim(R/J));
print("JAC_DEG="|toString degree(R/J));
print("HS="|toString hilbertSeries(R/J,Reduce=>true));
print("GB_GENERATORS="|toString numColumns gens G);
assert(dim(R/J)==2);assert(degree(R/J)==40);
print "RESONANCE ACTUAL DIVISOR NORMAL PASS";exit 0;
