S=QQ[p,t,r,c,a,b,d,MonomialOrder=>Eliminate 1];
J=ideal(p^2-a,p*r-b,p*t-d,2*c*p-r^2+25*t^2);
K=eliminate({p},J);
T=QQ[t,r,c,a,b,d];
KK=sub(K,T);
I=ideal(2*c*d - r^2*t + 25*t^3,2*b*c - r^3 + 25*r*t^2,a*t^2 - d^2,a*r*t - b*d,2*a^2*c*t - b^2*d + 25*d^3,b*t - d*r,2*a*c - b*r + 25*d*t,a*r^2 - b^2,2*a^2*c*r - b^3 + 25*b*d^2,-2*a*b*c + b^2*r - 25*d^2*r,4*a^3*c^2 - b^4 + 50*b^2*d^2 - 625*d^4);
assert(isSubset(I,KK) and isSubset(KK,I));
assert(dim(T/I)==3);
print "INDEPENDENT_B16_ENTIRE_ELIMINATION_PASS";
