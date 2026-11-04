// ONE saved matrix-vector product over QQ[a]/g13; integer pseudo-division.
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
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
void readp(istream&in,P&p){long n;in>>n;assert(in&&n>=0&&n<=150000);Z z;vector<long>seen;for(long j=0;j<n;j++){long k;string c;in>>k>>c;assert(in&&k>=0&&find(seen.begin(),seen.end(),k)==seen.end());seen.push_back(k);assert(!fmpz_set_str(z.v,c.c_str(),10)&&!fmpz_is_zero(z.v));fmpz_poly_set_coeff_fmpz(p.v,k,z.v);}}
void jp(ostream&o,const P&p){o<<'[';bool first=true;Z z;for(long i=fmpz_poly_degree(p.v);i>=0;i--){fmpz_poly_get_coeff_fmpz(z.v,p.v,i);if(fmpz_is_zero(z.v))continue;if(!first)o<<',';first=false;char*s=fmpz_get_str(nullptr,10,z.v);o<<'['<<i<<",\""<<s<<"\"]";flint_free(s);assert(o.tellp()<24*1024*1024);}o<<']';}
void jz(ostream&o,const fmpz_t z){char*s=fmpz_get_str(nullptr,10,z);o<<'"'<<s<<'"';flint_free(s);}
int main(int argc,char**argv){
 assert(argc==3);auto start=chrono::steady_clock::now();ifstream in(argv[1]);vector<vector<P>>N(16,vector<P>(16));vector<P>V(16),raw(16),quo(16),rem(16);P modulus;readp(in,modulus);assert(fmpz_poly_degree(modulus.v)==13);for(auto&r:N)for(auto&p:r){readp(in,p);assert(fmpz_poly_degree(p.v)<13);}for(auto&p:V){readp(in,p);assert(fmpz_poly_degree(p.v)<13);}string extra;assert(!(in>>extra));
 vector<ulong>powers(16);long peakbits=0,nonzero=0;Z lead;fmpz_set(lead.v,fmpz_poly_lead(modulus.v));
 for(int i=0;i<16;i++){
  P t;for(int j=0;j<16;j++){fmpz_poly_mul(t.v,N[i][j].v,V[j].v);fmpz_poly_add(raw[i].v,raw[i].v,t.v);}
  fmpz_poly_pseudo_divrem(quo[i].v,rem[i].v,&powers[i],raw[i].v,modulus.v);assert(powers[i]<=12&&fmpz_poly_degree(rem[i].v)<13);Z factor;fmpz_pow_ui(factor.v,lead.v,powers[i]);P lhs,rhs;fmpz_poly_scalar_mul_fmpz(lhs.v,raw[i].v,factor.v);fmpz_poly_mul(rhs.v,modulus.v,quo[i].v);fmpz_poly_add(rhs.v,rhs.v,rem[i].v);assert(fmpz_poly_equal(lhs.v,rhs.v));
  peakbits=max(peakbits,labs((long)fmpz_poly_max_bits(rem[i].v)));nonzero+=!fmpz_poly_is_zero(rem[i].v);assert(peakbits<=12629);
 }
 {ofstream o(string(argv[2])+".W-checkpoint.json");o<<"{\"all16_integer_pseudoreductions_checked\":true,\"rows\":[";for(int i=0;i<16;i++){if(i)o<<',';o<<"{\"raw_integer_product\":";jp(o,raw[i]);o<<",\"g13_quotient\":";jp(o,quo[i]);o<<",\"remainder\":";jp(o,rem[i]);o<<",\"leading_coefficient_power\":"<<powers[i]<<'}';}o<<"],\"nonzero_remainder_components\":"<<nonzero<<",\"max_remainder_coefficient_bits\":"<<peakbits<<"}\n";}
 P g=modulus;ofstream o(argv[2]);o<<"{\"gcd_chain\":[";bool first=true;int steps=0;
 for(int i=0;i<16;i++){
  if(fmpz_poly_is_zero(rem[i].v))continue;assert(chrono::duration<double>(chrono::steady_clock::now()-start).count()<52);P gn,qg,qr,S,T;Z R;fmpz_poly_gcd_subresultant(gn.v,g.v,rem[i].v);if(fmpz_sgn(fmpz_poly_lead(gn.v))<0)fmpz_poly_neg(gn.v,gn.v);assert(fmpz_poly_divides(qg.v,g.v,gn.v)&&fmpz_poly_divides(qr.v,rem[i].v,gn.v));
  if(fmpz_poly_degree(qg.v)==0){fmpz_poly_get_coeff_fmpz(R.v,qg.v,0);fmpz_poly_one(S.v);}
  else if(fmpz_poly_degree(qr.v)==0){fmpz_poly_get_coeff_fmpz(R.v,qr.v,0);fmpz_poly_one(T.v);}
  else fmpz_poly_xgcd_modular(R.v,S.v,T.v,qg.v,qr.v);
  assert(!fmpz_is_zero(R.v));P a,b,c;fmpz_poly_mul(a.v,S.v,g.v);fmpz_poly_mul(b.v,T.v,rem[i].v);fmpz_poly_add(a.v,a.v,b.v);fmpz_poly_scalar_mul_fmpz(c.v,gn.v,R.v);assert(fmpz_poly_equal(a.v,c.v));
  if(!first)o<<',';first=false;o<<"{\"component\":"<<i<<",\"previous_gcd\":";jp(o,g);o<<",\"new_gcd\":";jp(o,gn);o<<",\"previous_quotient\":";jp(o,qg);o<<",\"component_quotient\":";jp(o,qr);o<<",\"Bezout_S\":";jp(o,S);o<<",\"Bezout_T\":";jp(o,T);o<<",\"Bezout_constant\":";jz(o,R.v);o<<'}';steps++;g=gn;if(fmpz_poly_degree(g.v)==0)break;
 }
 o<<"],\"final_gcd\":";jp(o,g);o<<",\"gcd_is_QQunit\":"<<(fmpz_poly_degree(g.v)==0?"true":"false")<<",\"all_W_components_zero\":"<<(nonzero==0?"true":"false")<<",\"nonzero_W_components\":"<<nonzero<<",\"gcd_steps\":"<<steps<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";assert(o.tellp()<24*1024*1024);cerr<<"ALL16_W_PSEUDOREDUCTIONS_AND_GCD_CHAIN_PASS\n";
}
