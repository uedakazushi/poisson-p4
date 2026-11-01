#pragma once
// Exact QQ(e,d) arithmetic used by the quartic certificate.
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cassert>
#include <flint/fmpq_mpoly.h>
using namespace std;fmpq_mpoly_ctx_t ctx;
struct Rat{fmpq_mpoly_struct n[1],d[1];Rat(){fmpq_mpoly_init(n,ctx);fmpq_mpoly_init(d,ctx);fmpq_mpoly_one(d,ctx);}Rat(const Rat&o):Rat(){fmpq_mpoly_set(n,o.n,ctx);fmpq_mpoly_set(d,o.d,ctx);}Rat(Rat&&o):Rat(){fmpq_mpoly_swap(n,o.n,ctx);fmpq_mpoly_swap(d,o.d,ctx);}~Rat(){fmpq_mpoly_clear(n,ctx);fmpq_mpoly_clear(d,ctx);}Rat&operator=(Rat o){fmpq_mpoly_swap(n,o.n,ctx);fmpq_mpoly_swap(d,o.d,ctx);return *this;}bool zero()const{return fmpq_mpoly_is_zero(n,ctx);}void reduce(){if(zero()){fmpq_mpoly_one(d,ctx);return;}fmpq_mpoly_t g,a,b;fmpq_mpoly_init(g,ctx);fmpq_mpoly_init(a,ctx);fmpq_mpoly_init(b,ctx);assert(fmpq_mpoly_gcd_cofactors(g,a,b,n,d,ctx));fmpq_mpoly_swap(n,a,ctx);fmpq_mpoly_swap(d,b,ctx);fmpq_t lc;fmpq_init(lc);fmpq_mpoly_get_term_coeff_fmpq(lc,d,0,ctx);fmpq_mpoly_scalar_div_fmpq(n,n,lc,ctx);fmpq_mpoly_scalar_div_fmpq(d,d,lc,ctx);fmpq_clear(lc);fmpq_mpoly_clear(g,ctx);fmpq_mpoly_clear(a,ctx);fmpq_mpoly_clear(b,ctx);}};
void monic(Rat&c){if(c.zero()){fmpq_mpoly_one(c.d,ctx);return;}fmpq_t lc;fmpq_init(lc);fmpq_mpoly_get_term_coeff_fmpq(lc,c.d,0,ctx);fmpq_mpoly_scalar_div_fmpq(c.n,c.n,lc,ctx);fmpq_mpoly_scalar_div_fmpq(c.d,c.d,lc,ctx);fmpq_clear(lc);}
Rat mul(const Rat&a,const Rat&b){Rat c;if(a.zero()||b.zero())return c;Rat ab,ba;fmpq_mpoly_t g;fmpq_mpoly_init(g,ctx);assert(fmpq_mpoly_gcd_cofactors(g,ab.n,ab.d,a.n,b.d,ctx));assert(fmpq_mpoly_gcd_cofactors(g,ba.n,ba.d,b.n,a.d,ctx));fmpq_mpoly_mul(c.n,ab.n,ba.n,ctx);fmpq_mpoly_mul(c.d,ab.d,ba.d,ctx);fmpq_mpoly_clear(g,ctx);monic(c);return c;}
Rat divi(const Rat&a,const Rat&b){assert(!b.zero());Rat bi;fmpq_mpoly_set(bi.n,b.d,ctx);fmpq_mpoly_set(bi.d,b.n,ctx);return mul(a,bi);}
Rat sub(const Rat&a,const Rat&b){if(b.zero())return a;Rat c;if(a.zero()){c=b;fmpq_mpoly_neg(c.n,c.n,ctx);return c;}Rat den;fmpq_mpoly_t g,h,t,ng,gg;fmpq_mpoly_init(g,ctx);fmpq_mpoly_init(h,ctx);fmpq_mpoly_init(t,ctx);fmpq_mpoly_init(ng,ctx);fmpq_mpoly_init(gg,ctx);assert(fmpq_mpoly_gcd_cofactors(g,den.n,den.d,a.d,b.d,ctx));fmpq_mpoly_mul(c.n,a.n,den.d,ctx);fmpq_mpoly_mul(t,b.n,den.n,ctx);fmpq_mpoly_sub(c.n,c.n,t,ctx);if(!c.zero()){assert(fmpq_mpoly_gcd_cofactors(h,ng,gg,c.n,g,ctx));fmpq_mpoly_swap(c.n,ng,ctx);fmpq_mpoly_mul(c.d,den.n,den.d,ctx);fmpq_mpoly_mul(c.d,c.d,gg,ctx);monic(c);}fmpq_mpoly_clear(g,ctx);fmpq_mpoly_clear(h,ctx);fmpq_mpoly_clear(t,ctx);fmpq_mpoly_clear(ng,ctx);fmpq_mpoly_clear(gg,ctx);return c;}
void readpoly(istream&in,fmpq_mpoly_t p){fmpq_mpoly_zero(p,ctx);int n;in>>n;fmpq_t c;fmpq_init(c);ulong exp[2];for(int i=0;i<n;i++){string v;in>>v>>exp[0]>>exp[1];assert(fmpq_set_str(c,v.c_str(),10)==0);fmpq_canonicalise(c);fmpq_mpoly_push_term_fmpq_ui(p,c,exp,ctx);}fmpq_mpoly_sort_terms(p,ctx);fmpq_mpoly_combine_like_terms(p,ctx);fmpq_clear(c);}
void writepoly(ostream&out,const fmpq_mpoly_t p){int n=fmpq_mpoly_length(p,ctx);out<<n<<' ';fmpq_t c;fmpq_init(c);ulong exp[2];for(int i=0;i<n;i++){fmpq_mpoly_get_term_coeff_fmpq(c,p,i,ctx);char*str=fmpq_get_str(nullptr,10,c);fmpq_mpoly_get_term_exp_ui(exp,p,i,ctx);out<<str<<' '<<exp[0]<<' '<<exp[1]<<' ';flint_free(str);}out<<'\n';fmpq_clear(c);}
using Row=map<int,Rat>;
