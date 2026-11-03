R=QQ[a,b,c,d,e,f,w,u,v];
I=ideal(a,b,c,d,e,f,w*u,w*v,u*v,u^3);
assert(numgens trim I==10);
S=QQ[a,b,c,d,e,f,w,v];
J=ideal(a,b,c,d,e,f,w*v,v^2);
assert(numgens trim J==8);
print "ODD_DISCRIMINANT_MINIMAL_GENERATORS_PASS";
exit 0;
