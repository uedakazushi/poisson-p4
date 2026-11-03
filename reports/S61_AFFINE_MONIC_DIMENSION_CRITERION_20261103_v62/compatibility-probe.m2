A=QQ[ss1,ss2,ss3,ss4,ss5,cc];K=frac A;R=K[x,y];
F=cc^2-x*cc-y;I=ideal(x*x-2,y-3,F);G=gb(I,ChangeMatrix=>true);assert(1_R%G==0);C=matrix{{1_R}}//gens I;assert((gens I)*C==matrix{{1_R}});
J=ideal(x*y-1,F);H=gb(J,ChangeMatrix=>true);assert(1_R%H!=0);D=(gens H)//gens J;assert((gens J)*D==gens H);assert(numColumns basis(R/J)==2);
N=ideal(x*x,y,F);P=gb(N,ChangeMatrix=>true);assert(1_R%P==0);E=matrix{{1_R}}//gens N;assert((gens N)*E==matrix{{1_R}});
print "MONIC_DETECTOR_FINITE_CURVE_AND_NILPOTENT_TRACES_PASS";stdio << flush;exit 0;
