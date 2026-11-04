// ONE m2/m3 scalar job; optional B21 coefficients only when both zero.
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
void readp(istream&f,P&p){long n;f>>n;assert(f&&n>=0&&n<150000);Z z;for(long j=0;j<n;j++){long a;string s;f>>a>>s;assert(f&&a>=0&&!fmpz_set_str(z.v,s.c_str(),10));fmpz_poly_set_coeff_fmpz(p.v,a,z.v);}}
void jp(ostream&o,const P&p){o<<'[';bool first=true;Z z;for(long a=fmpz_poly_degree(p.v);a>=0;a--){fmpz_poly_get_coeff_fmpz(z.v,p.v,a);if(fmpz_is_zero(z.v))continue;if(!first)o<<',';first=false;char*s=fmpz_get_str(nullptr,10,z.v);o<<'['<<a<<",\""<<s<<"\"]";flint_free(s);assert(o.tellp()<24*1024*1024);}o<<']';}
void jz(ostream&o,const fmpz_t z){char*s=fmpz_get_str(nullptr,10,z);o<<'"'<<s<<'"';flint_free(s);}
void normal(P&num,fmpz_t den,const P&r,const fmpz_t rawden){fmpq_poly_t q;fmpq_poly_init(q);fmpq_poly_set_fmpz_poly(q,r.v);fmpq_poly_scalar_div_fmpz(q,q,rawden);fmpq_poly_get_numerator(num.v,q);fmpq_poly_get_denominator(den,q);fmpq_poly_clear(q);}
int main(int argc,char**argv){
 assert(argc==3);auto start=chrono::steady_clock::now();ifstream in(argv[1]);P g,tn;readp(in,g);readp(in,tn);assert(fmpz_poly_degree(g.v)==13&&fmpz_poly_degree(tn.v)<13);Z td,L;string s;in>>s;assert(!fmpz_set_str(td.v,s.c_str(),10));fmpz_set(L.v,fmpz_poly_lead(g.v));vector<P>power(24);fmpz_poly_one(power[0].v);for(int j=1;j<24;j++)fmpz_poly_mul(power[j].v,power[j-1].v,tn.v);
 ofstream o(argv[2]);o<<"{\"scalars\":[";bool first=true,bothzero=true,unit=false;int evaluated=0;
 for(int index=2;index<=3;index++){
  int count,ee;in>>count>>ee>>s;Z D;assert(in&&!fmpz_set_str(D.v,s.c_str(),10));P raw;
  for(int j=0;j<count;j++){long a,e;in>>a>>e>>s;assert(in&&a>=0&&e>=0&&e<=ee);Z c,factor;assert(!fmpz_set_str(c.v,s.c_str(),10));fmpz_pow_ui(factor.v,td.v,ee-e);fmpz_mul(c.v,c.v,factor.v);P term;fmpz_poly_shift_left(term.v,power[e].v,a);fmpz_poly_scalar_mul_fmpz(term.v,term.v,c.v);fmpz_poly_add(raw.v,raw.v,term.v);}
  Z common,lf,rd;fmpz_pow_ui(common.v,td.v,ee);fmpz_mul(common.v,common.v,D.v);P q,r,numer;ulong lp;fmpz_poly_pseudo_divrem(q.v,r.v,&lp,raw.v,g.v);fmpz_pow_ui(lf.v,L.v,lp);P lhs,rhs;fmpz_poly_scalar_mul_fmpz(lhs.v,raw.v,lf.v);fmpz_poly_mul(rhs.v,g.v,q.v);fmpz_poly_add(rhs.v,rhs.v,r.v);assert(fmpz_poly_equal(lhs.v,rhs.v));fmpz_mul(rd.v,common.v,lf.v);Z den;normal(numer,den.v,r,rd.v);assert(fmpz_poly_degree(numer.v)<13);P gg,S,T,prim;Z cnum,R;
  bool zero=fmpz_poly_is_zero(numer.v);bothzero=bothzero&&zero;
  if(!zero){fmpz_poly_content(cnum.v,numer.v);fmpz_poly_scalar_divexact_fmpz(prim.v,numer.v,cnum.v);fmpz_poly_gcd_subresultant(gg.v,prim.v,g.v);if(fmpz_poly_degree(gg.v)==0){fmpz_poly_xgcd_modular(R.v,T.v,S.v,g.v,prim.v);assert(!fmpz_is_zero(R.v));fmpz_poly_mul(lhs.v,S.v,prim.v);fmpz_poly_mul(rhs.v,T.v,g.v);fmpz_poly_add(lhs.v,lhs.v,rhs.v);fmpz_poly_set_fmpz(rhs.v,R.v);assert(fmpz_poly_equal(lhs.v,rhs.v));unit=true;}}
  if(!first)o<<',';first=false;evaluated++;o<<"{\"mixed_index\":"<<index<<",\"raw_numerator\":";jp(o,raw);o<<",\"common_source_denominator\":";jz(o,common.v);o<<",\"g13_quotient\":";jp(o,q);o<<",\"remainder\":";jp(o,r);o<<",\"L_power\":"<<lp<<",\"numerator\":";jp(o,numer);o<<",\"denominator\":";jz(o,den.v);o<<",\"scalar_zero\":"<<(zero?"true":"false")<<",\"unit_modg13\":"<<(unit?"true":"false");
  if(!zero){o<<",\"numerator_content\":";jz(o,cnum.v);o<<",\"primitive_numerator\":";jp(o,prim);o<<",\"gcd_with_g13\":";jp(o,gg);}
  if(unit){o<<",\"unit_S\":";jp(o,S);o<<",\"unit_T\":";jp(o,T);o<<",\"unit_constant\":";jz(o,R.v);}
  o<<'}';assert(chrono::duration<double>(chrono::steady_clock::now()-start).count()<50);if(unit)break;
 }
 o<<"],\"any_scalar_unit\":"<<(unit?"true":"false")<<",\"both_scalars_zero\":"<<((evaluated==2&&bothzero)?"true":"false")<<",\"scalar_count\":"<<evaluated;
 // B input follows both scalar records; never read/use it after early unit stop.
 bool prep=false;long scalarbytes=(long)o.tellp();
 if(evaluated==2&&bothzero&&scalarbytes+8165976+100000<32*1024*1024&&chrono::duration<double>(chrono::steady_clock::now()-start).count()<40){
  in>>s;Z DB,common;assert(in&&!fmpz_set_str(DB.v,s.c_str(),10));fmpz_pow_ui(common.v,td.v,4);fmpz_mul(common.v,common.v,DB.v);
  ofstream b(string(argv[2])+".B21.json");b<<"{\"source_prepared_B7x3\":true,\"d_degree3_noH4reduction\":true,\"entries\":[";
  for(int entry=0;entry<21;entry++){long count;in>>count;assert(in&&count>=0);vector<P>raw(4);for(long j=0;j<count;j++){long a,e,d;in>>a>>e>>d>>s;assert(in&&a>=0&&e>=0&&e<=4&&d>=0&&d<=3);Z c,factor;assert(!fmpz_set_str(c.v,s.c_str(),10));fmpz_pow_ui(factor.v,td.v,4-e);fmpz_mul(c.v,c.v,factor.v);P term;fmpz_poly_shift_left(term.v,power[e].v,a);fmpz_poly_scalar_mul_fmpz(term.v,term.v,c.v);fmpz_poly_add(raw[d].v,raw[d].v,term.v);}
   if(entry)b<<',';b<<'[';for(int d=0;d<4;d++){P q,r,num;ulong lp;fmpz_poly_pseudo_divrem(q.v,r.v,&lp,raw[d].v,g.v);Z f,rd,den;fmpz_pow_ui(f.v,L.v,lp);fmpz_mul(rd.v,common.v,f.v);P x,y;fmpz_poly_scalar_mul_fmpz(x.v,raw[d].v,f.v);fmpz_poly_mul(y.v,g.v,q.v);fmpz_poly_add(y.v,y.v,r.v);assert(fmpz_poly_equal(x.v,y.v));normal(num,den.v,r,rd.v);if(d)b<<',';b<<"{\"numerator\":";jp(b,num);b<<",\"denominator\":";jz(b,den.v);b<<'}';}b<<']';assert(b.tellp()+scalarbytes+100000<32*1024*1024);}
  b<<"],\"matrix_Gaussian_NOT_computed\":true}\n";prep=true;
 }
 o<<",\"B21_coefficients_prepared\":"<<(prep?"true":"false")<<",\"no_H4_reduction_Gaussian_or_rightinverse\":true,\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";assert(o.tellp()<24*1024*1024);cerr<<"MIXED23_SCALARS_AND_OPTIONAL_B21_INPUT_ONLY_PASS\n";
}
