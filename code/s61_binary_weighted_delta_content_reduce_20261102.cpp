#include "s61_binary_deeper_quartic_field_20261101.hpp"
#include <chrono>
int main(int argc,char**argv){assert(argc==3);fmpq_mpoly_ctx_init(ctx,2,ORD_LEX);{ifstream in(argv[1]);int n,m;in>>n>>m;vector<Row>A(n);for(auto&r:A){int nn;in>>nn;for(int k=0;k<nn;k++){int j;Rat a;in>>j;readpoly(in,a.n);r.emplace(j,std::move(a));}}
 auto start=chrono::steady_clock::now();int valuation=0;ofstream meta(string(argv[2])+".meta");
 for(int pass=0;pass<2;pass++)for(int idx=0;idx<(pass?m:n);idx++){
  vector<fmpq_mpoly_struct*>v;
  if(!pass){for(auto&[j,a]:A[idx])v.push_back(a.n);}else for(auto&r:A){auto it=r.find(idx);if(it!=r.end())v.push_back(it->second.n);}
  assert(!v.empty());Rat g,tmp;
  for(auto p:v){assert(fmpq_mpoly_gcd(tmp.n,g.n,p,ctx));fmpq_mpoly_swap(g.n,tmp.n,ctx);if(fmpq_mpoly_is_one(g.n,ctx))break;}
  ulong e[2];long ord=100000;for(long k=0;k<fmpq_mpoly_length(g.n,ctx);k++){fmpq_mpoly_get_term_exp_ui(e,g.n,k,ctx);ord=min(ord,long(e[1]));}
  valuation+=ord;meta<<pass<<' '<<idx<<' '<<ord<<'\n';writepoly(meta,g.n);
  for(auto p:v){assert(fmpq_mpoly_divides(tmp.n,p,g.n,ctx));fmpq_mpoly_swap(p,tmp.n,ctx);}
  fmpz_t gn,dd;fmpz_init(gn);fmpz_init(dd);fmpz_zero(gn);fmpz_one(dd);
  for(auto p:v){fmpz_gcd(gn,gn,fmpq_numref(p->content));fmpz_lcm(dd,dd,fmpq_denref(p->content));}
  fmpq_t c;fmpq_init(c);fmpz_set(fmpq_numref(c),dd);fmpz_set(fmpq_denref(c),gn);fmpq_canonicalise(c);for(auto p:v)fmpq_mpoly_scalar_mul_fmpq(p,p,c,ctx);fmpq_clear(c);fmpz_clear(gn);fmpz_clear(dd);
 }
 ofstream out(argv[2]),lead(string(argv[2])+".leading");out<<n<<' '<<m<<'\n';lead<<n<<' '<<m<<'\n';fmpq_t c;fmpq_init(c);ulong exp[2];long terms=0;
 for(auto&r:A){out<<r.size()<<'\n';Row l;for(auto&[j,v]:r){out<<j<<'\n';writepoly(out,v.n);terms+=fmpq_mpoly_length(v.n,ctx);Rat a;for(long k=0;k<fmpq_mpoly_length(v.n,ctx);k++){fmpq_mpoly_get_term_exp_ui(exp,v.n,k,ctx);if(exp[1])continue;fmpq_mpoly_get_term_coeff_fmpq(c,v.n,k,ctx);fmpq_mpoly_push_term_fmpq_ui(a.n,c,exp,ctx);}fmpq_mpoly_sort_terms(a.n,ctx);fmpq_mpoly_combine_like_terms(a.n,ctx);if(!a.zero())l.emplace(j,std::move(a));}Rat g,tmp;for(auto&[j,v]:l){assert(fmpq_mpoly_gcd(tmp.n,g.n,v.n,ctx));fmpq_mpoly_swap(g.n,tmp.n,ctx);}lead<<l.size()<<'\n';for(auto&[j,v]:l){assert(fmpq_mpoly_divides(tmp.n,v.n,g.n,ctx));lead<<j<<'\n';writepoly(lead,tmp.n);}}
 fmpq_clear(c);meta<<"VALUATION "<<valuation<<'\n';cerr<<"CONTENT_EXACT_PASS EXTRA_VALUATION "<<valuation<<" TERMS "<<terms<<" SECONDS "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
