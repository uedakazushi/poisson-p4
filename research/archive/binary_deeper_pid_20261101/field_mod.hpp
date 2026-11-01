
// Arithmetic in F_32003(e,d), for the explicitly labeled exploratory computations.
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cassert>
#include <flint/nmod_mpoly.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/ulong_extras.h>
using namespace std;nmod_mpoly_ctx_t ctx;ulong prime=32003;
struct Rat{nmod_mpoly_struct n[1],d[1];Rat(){nmod_mpoly_init(n,ctx);nmod_mpoly_init(d,ctx);nmod_mpoly_one(d,ctx);}Rat(const Rat&o):Rat(){nmod_mpoly_set(n,o.n,ctx);nmod_mpoly_set(d,o.d,ctx);}Rat(Rat&&o):Rat(){nmod_mpoly_swap(n,o.n,ctx);nmod_mpoly_swap(d,o.d,ctx);}~Rat(){nmod_mpoly_clear(n,ctx);nmod_mpoly_clear(d,ctx);}Rat&operator=(Rat o){nmod_mpoly_swap(n,o.n,ctx);nmod_mpoly_swap(d,o.d,ctx);return *this;}bool zero()const{return nmod_mpoly_is_zero(n,ctx);}void reduce(){if(zero()){nmod_mpoly_one(d,ctx);return;}nmod_mpoly_t g,a,b;nmod_mpoly_init(g,ctx);nmod_mpoly_init(a,ctx);nmod_mpoly_init(b,ctx);assert(nmod_mpoly_gcd_cofactors(g,a,b,n,d,ctx));nmod_mpoly_swap(n,a,ctx);nmod_mpoly_swap(d,b,ctx);ulong lc=n_invmod(nmod_mpoly_get_term_coeff_ui(d,0,ctx),prime);nmod_mpoly_scalar_mul_ui(n,n,lc,ctx);nmod_mpoly_scalar_mul_ui(d,d,lc,ctx);nmod_mpoly_clear(g,ctx);nmod_mpoly_clear(a,ctx);nmod_mpoly_clear(b,ctx);}};
void monic(Rat&c){if(c.zero()){nmod_mpoly_one(c.d,ctx);return;}ulong lc=n_invmod(nmod_mpoly_get_term_coeff_ui(c.d,0,ctx),prime);nmod_mpoly_scalar_mul_ui(c.n,c.n,lc,ctx);nmod_mpoly_scalar_mul_ui(c.d,c.d,lc,ctx);}
Rat mul(const Rat&a,const Rat&b){Rat c;if(a.zero()||b.zero())return c;Rat ab,ba;nmod_mpoly_t g;nmod_mpoly_init(g,ctx);assert(nmod_mpoly_gcd_cofactors(g,ab.n,ab.d,a.n,b.d,ctx));assert(nmod_mpoly_gcd_cofactors(g,ba.n,ba.d,b.n,a.d,ctx));nmod_mpoly_mul(c.n,ab.n,ba.n,ctx);nmod_mpoly_mul(c.d,ab.d,ba.d,ctx);nmod_mpoly_clear(g,ctx);monic(c);return c;}
Rat divi(const Rat&a,const Rat&b){assert(!b.zero());Rat bi;nmod_mpoly_set(bi.n,b.d,ctx);nmod_mpoly_set(bi.d,b.n,ctx);return mul(a,bi);}
Rat sub(const Rat&a,const Rat&b){if(b.zero())return a;Rat c;if(a.zero()){c=b;nmod_mpoly_neg(c.n,c.n,ctx);return c;}Rat den;nmod_mpoly_t g,h,t,ng,gg;nmod_mpoly_init(g,ctx);nmod_mpoly_init(h,ctx);nmod_mpoly_init(t,ctx);nmod_mpoly_init(ng,ctx);nmod_mpoly_init(gg,ctx);assert(nmod_mpoly_gcd_cofactors(g,den.n,den.d,a.d,b.d,ctx));nmod_mpoly_mul(c.n,a.n,den.d,ctx);nmod_mpoly_mul(t,b.n,den.n,ctx);nmod_mpoly_sub(c.n,c.n,t,ctx);if(!c.zero()){assert(nmod_mpoly_gcd_cofactors(h,ng,gg,c.n,g,ctx));nmod_mpoly_swap(c.n,ng,ctx);nmod_mpoly_mul(c.d,den.n,den.d,ctx);nmod_mpoly_mul(c.d,c.d,gg,ctx);monic(c);}nmod_mpoly_clear(g,ctx);nmod_mpoly_clear(h,ctx);nmod_mpoly_clear(t,ctx);nmod_mpoly_clear(ng,ctx);nmod_mpoly_clear(gg,ctx);return c;}
void readpoly(istream&in,nmod_mpoly_t p){nmod_mpoly_zero(p,ctx);int n;in>>n;fmpq_t c;fmpq_init(c);ulong exp[2];for(int i=0;i<n;i++){string v;in>>v>>exp[0]>>exp[1];assert(!fmpq_set_str(c,v.c_str(),10));fmpq_canonicalise(c);ulong a=fmpz_fdiv_ui(fmpq_numref(c),prime),b=fmpz_fdiv_ui(fmpq_denref(c),prime);assert(b);ulong q=n_mulmod2_preinv(a,n_invmod(b,prime),prime,n_preinvert_limb(prime));nmod_mpoly_push_term_ui_ui(p,q,exp,ctx);}nmod_mpoly_sort_terms(p,ctx);nmod_mpoly_combine_like_terms(p,ctx);fmpq_clear(c);}
void writepoly(ostream&out,const nmod_mpoly_t p){int n=nmod_mpoly_length(p,ctx);out<<n<<' ';ulong exp[2];for(int i=0;i<n;i++){nmod_mpoly_get_term_exp_ui(exp,p,i,ctx);out<<nmod_mpoly_get_term_coeff_ui(p,i,ctx)<<' '<<exp[0]<<' '<<exp[1]<<' ';}out<<'\n';}
using Row=map<int,Rat>;
