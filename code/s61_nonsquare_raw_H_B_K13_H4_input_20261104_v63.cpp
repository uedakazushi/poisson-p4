// ONE H4 coefficient preparation. Integer substitutions and exact remainders.
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
void readp(istream&in,P&p){long n;in>>n;assert(n>=0&&n<150000);Z z;for(long i=0;i<n;i++){long d;string s;in>>d>>s;assert(in&&d>=0&&!fmpz_set_str(z.v,s.c_str(),10));fmpz_poly_set_coeff_fmpz(p.v,d,z.v);}}
void jp(ostream&o,const P&p){o<<'[';bool first=true;Z z;for(long i=fmpz_poly_degree(p.v);i>=0;i--){fmpz_poly_get_coeff_fmpz(z.v,p.v,i);if(fmpz_is_zero(z.v))continue;if(!first)o<<',';first=false;char*s=fmpz_get_str(nullptr,10,z.v);o<<'['<<i<<",\""<<s<<"\"]";flint_free(s);assert(o.tellp()<28*1024*1024);}o<<']';}
void jz(ostream&o,const fmpz_t z){char*s=fmpz_get_str(nullptr,10,z);o<<'"'<<s<<'"';flint_free(s);}
int main(int argc,char**argv){
 assert(argc==3);auto start=chrono::steady_clock::now();ifstream in(argv[1]);P g,tn;readp(in,g);readp(in,tn);assert(fmpz_poly_degree(g.v)==13&&fmpz_poly_degree(tn.v)<13);Z td,DH,g0,L;string s;in>>s;assert(!fmpz_set_str(td.v,s.c_str(),10));in>>s;assert(!fmpz_set_str(DH.v,s.c_str(),10));fmpz_poly_get_coeff_fmpz(g0.v,g.v,0);assert(!fmpz_is_zero(g0.v));fmpz_set(L.v,fmpz_poly_lead(g.v));
 P ain;fmpz_poly_shift_right(ain.v,g.v,1);fmpz_poly_neg(ain.v,ain.v);vector<P>ap(6),tp(9),raw(5),q(5),r(5),normal(5);fmpz_poly_one(ap[0].v);fmpz_poly_one(tp[0].v);for(int i=1;i<=5;i++)fmpz_poly_mul(ap[i].v,ap[i-1].v,ain.v);for(int i=1;i<=8;i++)fmpz_poly_mul(tp[i].v,tp[i-1].v,tn.v);
 long count;in>>count;assert(count==74);for(long j=0;j<count;j++){int k,e,d;in>>k>>e>>d>>s;assert(in&&0<=k&&k<=5&&0<=e&&e<=8&&0<=d&&d<=4);Z c,z1,z2;assert(!fmpz_set_str(c.v,s.c_str(),10));fmpz_pow_ui(z1.v,g0.v,5-k);fmpz_pow_ui(z2.v,td.v,8-e);fmpz_mul(c.v,c.v,z1.v);fmpz_mul(c.v,c.v,z2.v);P term;fmpz_poly_mul(term.v,ap[k].v,tp[e].v);fmpz_poly_scalar_mul_fmpz(term.v,term.v,c.v);fmpz_poly_add(raw[d].v,raw[d].v,term.v);}
 Z common,z1,z2;fmpz_pow_ui(z1.v,g0.v,5);fmpz_pow_ui(z2.v,td.v,8);fmpz_mul(common.v,DH.v,z1.v);fmpz_mul(common.v,common.v,z2.v);vector<ulong>powers(5);vector<Z*>den;long maxbits=0;
 for(int i=0;i<5;i++){fmpz_poly_pseudo_divrem(q[i].v,r[i].v,&powers[i],raw[i].v,g.v);assert(powers[i]<=144&&fmpz_poly_degree(r[i].v)<13);Z lf;fmpz_pow_ui(lf.v,L.v,powers[i]);P lhs,rhs;fmpz_poly_scalar_mul_fmpz(lhs.v,raw[i].v,lf.v);fmpz_poly_mul(rhs.v,g.v,q[i].v);fmpz_poly_add(rhs.v,rhs.v,r[i].v);assert(fmpz_poly_equal(lhs.v,rhs.v));auto*dd=new Z;den.push_back(dd);Z rd;fmpz_mul(rd.v,common.v,lf.v);fmpq_poly_t pp;fmpq_poly_init(pp);fmpq_poly_set_fmpz_poly(pp,r[i].v);fmpq_poly_scalar_div_fmpz(pp,pp,rd.v);fmpq_poly_get_numerator(normal[i].v,pp);fmpq_poly_get_denominator(dd->v,pp);fmpq_poly_clear(pp);fmpz_poly_scalar_mul_fmpz(lhs.v,normal[i].v,rd.v);fmpz_poly_scalar_mul_fmpz(rhs.v,r[i].v,dd->v);assert(fmpz_poly_equal(lhs.v,rhs.v));maxbits=max(maxbits,labs((long)fmpz_poly_max_bits(normal[i].v)));}
 assert(fmpz_poly_is_one(normal[4].v)&&fmpz_is_one(den[4]->v));
 ofstream o(argv[2]);o<<"{\"common_source_denominator\":";jz(o,common.v);o<<",\"a_inverse_numerator\":";jp(o,ain);o<<",\"a_inverse_denominator\":";jz(o,g0.v);o<<",\"H4_coefficients\":[";
 for(int i=0;i<5;i++){if(i)o<<',';o<<"{\"d_degree\":"<<i<<",\"raw_common_numerator\":";jp(o,raw[i]);o<<",\"g13_quotient\":";jp(o,q[i]);o<<",\"remainder\":";jp(o,r[i]);o<<",\"L_power\":"<<powers[i]<<",\"numerator\":";jp(o,normal[i]);o<<",\"denominator\":";jz(o,den[i]->v);o<<'}';}
 o<<"],\"monic_H4_verified\":true,\"all5_source_reductions_and_normalizations_verified\":true,\"maximum_canonical_numerator_bits\":"<<maxbits<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";assert(o.tellp()<28*1024*1024);for(auto*d:den)delete d;cerr<<"H4_ONLY_EXACT_SOURCE_INPUT_PASS\n";
}
