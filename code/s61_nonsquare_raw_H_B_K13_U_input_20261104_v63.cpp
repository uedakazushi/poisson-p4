// ONE FULLU coefficient preparation, no H/U gcd and no matrix solve.
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
void readp(istream&in,P&p){long n;in>>n;assert(n>=0&&n<150000);Z z;for(long j=0;j<n;j++){long a;string s;in>>a>>s;assert(in&&a>=0&&!fmpz_set_str(z.v,s.c_str(),10));fmpz_poly_set_coeff_fmpz(p.v,a,z.v);}}
void jp(ostream&o,const P&p){o<<'[';bool first=true;Z z;for(long a=fmpz_poly_degree(p.v);a>=0;a--){fmpz_poly_get_coeff_fmpz(z.v,p.v,a);if(fmpz_is_zero(z.v))continue;if(!first)o<<',';first=false;char*s=fmpz_get_str(nullptr,10,z.v);o<<'['<<a<<",\""<<s<<"\"]";flint_free(s);assert(o.tellp()<28*1024*1024);}o<<']';}
void jz(ostream&o,const fmpz_t z){char*s=fmpz_get_str(nullptr,10,z);o<<'"'<<s<<'"';flint_free(s);}
int main(int argc,char**argv){
 assert(argc==3);auto start=chrono::steady_clock::now();ifstream in(argv[1]);P g,tn;readp(in,g);readp(in,tn);Z td,DU,L;string s;in>>s;assert(!fmpz_set_str(td.v,s.c_str(),10));in>>s;assert(!fmpz_set_str(DU.v,s.c_str(),10));fmpz_set(L.v,fmpz_poly_lead(g.v));vector<P>power(8),raw(4),q(4),r(4),normal(4);fmpz_poly_one(power[0].v);for(int j=1;j<=7;j++)fmpz_poly_mul(power[j].v,power[j-1].v,tn.v);
 long count;in>>count;assert(count==49);for(long j=0;j<count;j++){long a,e,d;in>>a>>e>>d>>s;assert(in&&a>=0&&a<=5&&e>=0&&e<=7&&d>=0&&d<=3);Z c,factor;assert(!fmpz_set_str(c.v,s.c_str(),10));fmpz_pow_ui(factor.v,td.v,7-e);fmpz_mul(c.v,c.v,factor.v);P term;fmpz_poly_shift_left(term.v,power[e].v,a);fmpz_poly_scalar_mul_fmpz(term.v,term.v,c.v);fmpz_poly_add(raw[d].v,raw[d].v,term.v);}
 Z common;fmpz_pow_ui(common.v,td.v,7);fmpz_mul(common.v,common.v,DU.v);vector<ulong>powers(4);vector<Z*>den;bool zero=true;
 for(int i=0;i<4;i++){fmpz_poly_pseudo_divrem(q[i].v,r[i].v,&powers[i],raw[i].v,g.v);assert(powers[i]<=77&&fmpz_poly_degree(r[i].v)<13);Z f,rd;fmpz_pow_ui(f.v,L.v,powers[i]);fmpz_mul(rd.v,common.v,f.v);P x,y;fmpz_poly_scalar_mul_fmpz(x.v,raw[i].v,f.v);fmpz_poly_mul(y.v,g.v,q[i].v);fmpz_poly_add(y.v,y.v,r[i].v);assert(fmpz_poly_equal(x.v,y.v));auto*d=new Z;den.push_back(d);fmpq_poly_t p;fmpq_poly_init(p);fmpq_poly_set_fmpz_poly(p,r[i].v);fmpq_poly_scalar_div_fmpz(p,p,rd.v);fmpq_poly_get_numerator(normal[i].v,p);fmpq_poly_get_denominator(d->v,p);fmpq_poly_clear(p);fmpz_poly_scalar_mul_fmpz(x.v,normal[i].v,rd.v);fmpz_poly_scalar_mul_fmpz(y.v,r[i].v,d->v);assert(fmpz_poly_equal(x.v,y.v));zero=zero&&fmpz_poly_is_zero(normal[i].v);}
 ofstream o(argv[2]);o<<"{\"common_source_denominator\":";jz(o,common.v);o<<",\"U_coefficients\":[";
 for(int i=0;i<4;i++){if(i)o<<',';o<<"{\"d_degree\":"<<i<<",\"raw_common_numerator\":";jp(o,raw[i]);o<<",\"g13_quotient\":";jp(o,q[i]);o<<",\"remainder\":";jp(o,r[i]);o<<",\"L_power\":"<<powers[i]<<",\"numerator\":";jp(o,normal[i]);o<<",\"denominator\":";jz(o,den[i]->v);o<<'}';}
 o<<"],\"all_U_coefficients_zero\":"<<(zero?"true":"false")<<",\"no_H4gcd_or_Gaussian\":true,\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";assert(o.tellp()<28*1024*1024);for(auto*d:den)delete d;cerr<<"FULLU_COEFFICIENT_INPUT_ONLY_PASS\n";
}
