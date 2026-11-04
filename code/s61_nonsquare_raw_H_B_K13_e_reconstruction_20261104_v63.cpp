// ONE companion eigenvalue reconstruction, with certified unit pivot.
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
struct P{fmpz_poly_t v;P(){fmpz_poly_init(v);}P(const P&x):P(){fmpz_poly_set(v,x.v);}P(P&&x):P(){fmpz_poly_swap(v,x.v);}~P(){fmpz_poly_clear(v);}P&operator=(P x){fmpz_poly_swap(v,x.v);return*this;}};
struct Z{fmpz_t v;Z(){fmpz_init(v);}~Z(){fmpz_clear(v);}};
void readp(istream&in,P&p){long n;in>>n;assert(in&&n>=0&&n<150000);Z z;for(long i=0;i<n;i++){long k;string c;in>>k>>c;assert(k>=0&&!fmpz_set_str(z.v,c.c_str(),10));fmpz_poly_set_coeff_fmpz(p.v,k,z.v);}}
void jp(ostream&o,const P&p){o<<'[';bool first=true;Z z;for(long i=fmpz_poly_degree(p.v);i>=0;i--){fmpz_poly_get_coeff_fmpz(z.v,p.v,i);if(fmpz_is_zero(z.v))continue;if(!first)o<<',';first=false;char*s=fmpz_get_str(nullptr,10,z.v);o<<'['<<i<<",\""<<s<<"\"]";flint_free(s);assert(o.tellp()<24*1024*1024);}o<<']';}
void jz(ostream&o,const fmpz_t z){char*s=fmpz_get_str(nullptr,10,z);o<<'"'<<s<<'"';flint_free(s);}
int main(int argc,char**argv){
 assert(argc==3);auto start=chrono::steady_clock::now();ifstream in(argv[1]);P g;readp(in,g);assert(fmpz_poly_degree(g.v)==13);vector<P>V(16),H(16);for(auto&p:V)readp(in,p);for(auto&p:H)readp(in,p);Z DH,DV;string s;in>>s;assert(!fmpz_set_str(DH.v,s.c_str(),10));in>>s;assert(!fmpz_set_str(DV.v,s.c_str(),10));
 int selected=-1;long minbits=100000000;vector<P>prim(16);vector<Z*>contents;
 for(int i=0;i<16;i++){auto*c=new Z;contents.push_back(c);fmpz_poly_content(c->v,V[i].v);if(fmpz_is_zero(c->v))continue;fmpz_poly_scalar_divexact_fmpz(prim[i].v,V[i].v,c->v);long b=labs((long)fmpz_poly_max_bits(prim[i].v));if(b<minbits){minbits=b;selected=i;}}
 assert(selected>=0);P unitg;fmpz_poly_gcd_subresultant(unitg.v,prim[selected].v,g.v);
 if(fmpz_poly_degree(unitg.v)>0){ofstream o(argv[2]);o<<"{\"pivot_unit\":false,\"selected\":"<<selected<<",\"proper_pivot_gcd\":";jp(o,unitg);o<<"}\n";for(auto*c:contents)delete c;return 0;}
 P S,T;Z R;fmpz_poly_xgcd_modular(R.v,T.v,S.v,g.v,prim[selected].v);assert(!fmpz_is_zero(R.v));P x,y,r;fmpz_poly_mul(x.v,S.v,prim[selected].v);fmpz_poly_mul(y.v,T.v,g.v);fmpz_poly_add(x.v,x.v,y.v);fmpz_poly_set_fmpz(r.v,R.v);assert(fmpz_poly_equal(x.v,r.v));
 vector<P>EVraw(16),EV(16),EVq(16);vector<ulong>powers(16);Z L;fmpz_set(L.v,fmpz_poly_lead(g.v));
 for(int i=0;i<16;i++){if(i)fmpz_poly_scalar_mul_fmpz(EVraw[i].v,V[i-1].v,DH.v);fmpz_poly_mul(x.v,H[i].v,V[15].v);fmpz_poly_sub(EVraw[i].v,EVraw[i].v,x.v);fmpz_poly_pseudo_divrem(EVq[i].v,EV[i].v,&powers[i],EVraw[i].v,g.v);Z f;fmpz_pow_ui(f.v,L.v,powers[i]);fmpz_poly_scalar_mul_fmpz(x.v,EVraw[i].v,f.v);fmpz_poly_mul(y.v,g.v,EVq[i].v);fmpz_poly_add(y.v,y.v,EV[i].v);assert(fmpz_poly_equal(x.v,y.v));}
 P traw,trem,tq;fmpz_poly_mul(traw.v,EV[selected].v,S.v);ulong tp;fmpz_poly_pseudo_divrem(tq.v,trem.v,&tp,traw.v,g.v);Z td,lf;fmpz_pow_ui(lf.v,L.v,powers[selected]+tp);fmpz_mul(td.v,DH.v,lf.v);fmpz_mul(td.v,td.v,contents[selected]->v);fmpz_mul(td.v,td.v,R.v);
 fmpq_poly_t t;fmpq_poly_init(t);fmpq_poly_set_fmpz_poly(t,trem.v);fmpq_poly_scalar_div_fmpz(t,t,td.v);fmpq_poly_canonicalise(t);P tn;Z den;fmpq_poly_get_numerator(tn.v,t);fmpq_poly_get_denominator(den.v,t);fmpq_poly_clear(t);
 vector<P>identityq(16);for(int i=0;i<16;i++){fmpz_pow_ui(lf.v,L.v,powers[i]);fmpz_mul(lf.v,lf.v,DH.v);fmpz_poly_mul(x.v,tn.v,V[i].v);fmpz_poly_scalar_mul_fmpz(x.v,x.v,lf.v);fmpz_poly_scalar_mul_fmpz(y.v,EV[i].v,den.v);fmpz_poly_sub(x.v,x.v,y.v);assert(fmpz_poly_divides(identityq[i].v,x.v,g.v));fmpz_poly_mul(y.v,identityq[i].v,g.v);assert(fmpz_poly_equal(x.v,y.v));}
 ofstream o(argv[2]);o<<"{\"pivot_unit\":true,\"selected\":"<<selected<<",\"pivot_primitive_bits\":"<<minbits<<",\"pivot_content\":";jz(o,contents[selected]->v);o<<",\"pivot_primitive\":";jp(o,prim[selected]);o<<",\"unit_S\":";jp(o,S);o<<",\"unit_T\":";jp(o,T);o<<",\"unit_constant\":";jz(o,R.v);o<<",\"e_numerator\":";jp(o,tn);o<<",\"e_denominator\":";jz(o,den.v);
 o<<",\"e_raw_numerator\":";jp(o,traw);o<<",\"e_reduction_quotient\":";jp(o,tq);o<<",\"e_reduction_remainder\":";jp(o,trem);o<<",\"e_reduction_L_power\":"<<tp<<",\"e_raw_denominator\":";jz(o,td.v);o<<",\"EV_rows\":[";
 for(int i=0;i<16;i++){if(i)o<<',';o<<"{\"raw\":";jp(o,EVraw[i]);o<<",\"quotient\":";jp(o,EVq[i]);o<<",\"remainder\":";jp(o,EV[i]);o<<",\"L_power\":"<<powers[i]<<",\"EV_equals_tV_quotient\":";jp(o,identityq[i]);o<<'}';}
 o<<"],\"all16_EV_equals_tV\":true,\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";assert(o.tellp()<24*1024*1024);for(auto*c:contents)delete c;cerr<<"UNIT_PIVOT_AND_ALL16_COMPANION_EIGEN_IDENTITIES_PASS\n";
}
