#pragma once
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <cassert>
#include <string>
#include <flint/fmpq_poly.h>
#include <flint/fmpq.h>
using namespace std;
struct URat{fmpq_poly_struct n[1],d[1];URat(){fmpq_poly_init(n);fmpq_poly_init(d);fmpq_poly_one(d);}URat(const URat&a):URat(){fmpq_poly_set(n,a.n);fmpq_poly_set(d,a.d);}URat(URat&&a):URat(){fmpq_poly_swap(n,a.n);fmpq_poly_swap(d,a.d);}URat&operator=(URat a){fmpq_poly_swap(n,a.n);fmpq_poly_swap(d,a.d);return *this;}~URat(){fmpq_poly_clear(n);fmpq_poly_clear(d);}bool zero()const{return fmpq_poly_is_zero(n);}};
void umonic(URat&a){if(a.zero()){fmpq_poly_one(a.d);return;}fmpq_t lc;fmpq_init(lc);fmpq_poly_get_coeff_fmpq(lc,a.d,fmpq_poly_degree(a.d));fmpq_poly_scalar_div_fmpq(a.n,a.n,lc);fmpq_poly_scalar_div_fmpq(a.d,a.d,lc);fmpq_clear(lc);}
URat umul(const URat&a,const URat&b){URat c;if(a.zero()||b.zero())return c;URat ab,ba;fmpq_poly_t g;fmpq_poly_init(g);fmpq_poly_gcd(g,a.n,b.d);assert(fmpq_poly_divides(ab.n,a.n,g));assert(fmpq_poly_divides(ab.d,b.d,g));fmpq_poly_gcd(g,b.n,a.d);assert(fmpq_poly_divides(ba.n,b.n,g));assert(fmpq_poly_divides(ba.d,a.d,g));fmpq_poly_mul(c.n,ab.n,ba.n);fmpq_poly_mul(c.d,ab.d,ba.d);fmpq_poly_clear(g);umonic(c);return c;}
URat udiv(const URat&a,const URat&b){assert(!b.zero());URat z;fmpq_poly_set(z.n,b.d);fmpq_poly_set(z.d,b.n);return umul(a,z);}
URat usub(const URat&a,const URat&b){if(b.zero())return a;URat c;if(a.zero()){c=b;fmpq_poly_neg(c.n,c.n);return c;}URat den;fmpq_poly_t g,t,h,q;for(auto p:{g,t,h,q})fmpq_poly_init(p);fmpq_poly_gcd(g,a.d,b.d);assert(fmpq_poly_divides(den.n,a.d,g));assert(fmpq_poly_divides(den.d,b.d,g));fmpq_poly_mul(c.n,a.n,den.d);fmpq_poly_mul(t,b.n,den.n);fmpq_poly_sub(c.n,c.n,t);if(!c.zero()){fmpq_poly_gcd(h,c.n,g);assert(fmpq_poly_divides(q,c.n,h));fmpq_poly_swap(c.n,q);assert(fmpq_poly_divides(q,g,h));fmpq_poly_mul(c.d,den.n,den.d);fmpq_poly_mul(c.d,c.d,q);umonic(c);}for(auto p:{g,t,h,q})fmpq_poly_clear(p);return c;}
void uread(istream&in,fmpq_poly_t p){int n;in>>n;fmpq_poly_zero(p);fmpq_t c;fmpq_init(c);for(int i=0;i<n;i++){string a;int e,d;in>>a>>e>>d;assert(d==0&&fmpq_set_str(c,a.c_str(),10)==0);fmpq_canonicalise(c);fmpq_poly_set_coeff_fmpq(p,e,c);}fmpq_clear(c);}
void uwrite(ostream&out,const fmpq_poly_t p){int count=0;fmpq_t c;fmpq_init(c);for(int i=0;i<fmpq_poly_length(p);i++){fmpq_poly_get_coeff_fmpq(c,p,i);if(!fmpq_is_zero(c))count++;}out<<count<<' ';for(int i=0;i<fmpq_poly_length(p);i++){fmpq_poly_get_coeff_fmpq(c,p,i);if(fmpq_is_zero(c))continue;char*s=fmpq_get_str(nullptr,10,c);out<<s<<' '<<i<<" 0 ";flint_free(s);}out<<'\n';fmpq_clear(c);}
using URow=map<int,URat>;
