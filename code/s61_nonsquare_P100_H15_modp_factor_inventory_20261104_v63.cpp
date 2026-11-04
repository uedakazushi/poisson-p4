// ONE existing-prime factor inventory; deterministic Rabin certificates.
#include <flint/flint.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
const ulong P=32003;
struct Poly{nmod_poly_t v;Poly(){nmod_poly_init(v,P);}Poly(const Poly&a):Poly(){nmod_poly_set(v,a.v);}~Poly(){nmod_poly_clear(v);}};
void jp(ostream&o,const Poly&a){o<<'[';for(long i=0;i<nmod_poly_length(a.v);i++){if(i)o<<',';o<<nmod_poly_get_coeff_ui(a.v,i);}o<<']';}
int main(int argc,char**argv){assert(argc==3);ifstream in(argv[1]);long n;in>>n;assert(n==16);Poly H;for(long i=0;i<n;i++){ulong a;in>>a;assert(a<P);nmod_poly_set_coeff_ui(H.v,i,a);}assert(nmod_poly_degree(H.v)==15&&nmod_poly_get_coeff_ui(H.v,15)==1);
 nmod_poly_factor_t fac;nmod_poly_factor_init(fac);ulong c=nmod_poly_factor(fac,H.v);assert(c==1);Poly product;nmod_poly_one(product.v);ofstream o(argv[2]);o<<"{\"prime\":32003,\"factor_algorithm\":\"nmod_poly_factor\",\"factor_scalar\":"<<c<<",\"Hbar\":";jp(o,H);o<<",\"factors\":[";
 for(long i=0;i<fac->num;i++){if(i)o<<',';Poly f;nmod_poly_factor_get_poly(f.v,fac,i);long deg=nmod_poly_degree(f.v),mult=fac->exp[i];assert(deg>=1);bool flag=nmod_poly_is_irreducible_rabin(f.v);assert(flag);Poly pwr;nmod_poly_pow(pwr.v,f.v,mult);nmod_poly_mul(product.v,product.v,pwr.v);
  Poly x,X;nmod_poly_set_coeff_ui(x.v,1,1);nmod_poly_rem(x.v,x.v,f.v);nmod_poly_set(X.v,x.v);vector<long> tests;long r=deg;for(long d=2;d<=r;d++)if(r%d==0){tests.push_back(deg/d);while(r%d==0)r/=d;}
  o<<"{\"degree\":"<<deg<<",\"multiplicity\":"<<mult<<",\"coefficients\":";jp(o,f);o<<",\"library_rabin_irreducible\":true,\"Frobenius_chain\":[";
  for(long k=1;k<=deg;k++){if(k>1)o<<',';nmod_poly_powmod_ui_binexp(X.v,X.v,P,f.v);jp(o,X);}
  assert(nmod_poly_equal(X.v,x.v));o<<"],\"terminal_x_p_degree_equals_x\":true,\"proper_subfield_tests\":[";
  for(long ti=0;ti<(long)tests.size();ti++){if(ti)o<<',';long k=tests[ti];Poly T,B,g,s,t,check,w;nmod_poly_set(T.v,x.v);for(long j=0;j<k;j++)nmod_poly_powmod_ui_binexp(T.v,T.v,P,f.v);nmod_poly_sub(B.v,T.v,x.v);nmod_poly_xgcd(g.v,s.v,t.v,B.v,f.v);assert(nmod_poly_is_one(g.v));nmod_poly_mul(check.v,s.v,B.v);nmod_poly_mul(w.v,t.v,f.v);nmod_poly_add(check.v,check.v,w.v);assert(nmod_poly_is_one(check.v));o<<"{\"frobenius_exponent\":"<<k<<",\"difference\":";jp(o,B);o<<",\"Bezout_difference\":";jp(o,s);o<<",\"Bezout_factor\":";jp(o,t);o<<",\"gcd_is1\":true}";}
  o<<"]}";
 }
 assert(nmod_poly_equal(product.v,H.v));Poly der,g;nmod_poly_derivative(der.v,H.v);nmod_poly_gcd(g.v,H.v,der.v);o<<"],\"full_product\":";jp(o,product);o<<",\"full_product_equals_Hbar\":true,\"derivative_gcd\":";jp(o,g);o<<",\"Hbar_squarefree\":"<<(nmod_poly_is_one(g.v)?"true":"false")<<",\"all_passed\":true}\n";nmod_poly_factor_clear(fac);flint_cleanup();}
