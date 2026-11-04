// Exact22-polynomial GCD and localizer-supported removal ONLY; no GB/factorization.
#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;fmpq_mpoly_ctx_t ctx;
struct P{fmpq_mpoly_t p;P(){fmpq_mpoly_init(p,ctx);}P(const P&v):P(){fmpq_mpoly_set(p,v.p,ctx);}P(P&&v):P(){fmpq_mpoly_swap(p,v.p,ctx);}P&operator=(P v){fmpq_mpoly_swap(p,v.p,ctx);return*this;}~P(){fmpq_mpoly_clear(p,ctx);}};
void readp(istream&f,P&p){long n;f>>n;assert(f&&n>0);fmpq_t c;fmpq_init(c);ulong ex[3];for(long i=0;i<n;i++){string s;f>>s>>ex[0]>>ex[1]>>ex[2];assert(f&&!fmpq_set_str(c,s.c_str(),10));fmpq_canonicalise(c);fmpq_mpoly_push_term_fmpq_ui(p.p,c,ex,ctx);}fmpq_mpoly_sort_terms(p.p,ctx);fmpq_mpoly_combine_like_terms(p.p,ctx);assert(fmpq_mpoly_length(p.p,ctx)==n);fmpq_clear(c);}
void monic(P&p){assert(!fmpq_mpoly_is_zero(p.p,ctx));fmpq_t c;fmpq_init(c);fmpq_mpoly_get_term_coeff_fmpq(c,p.p,0,ctx);fmpq_mpoly_scalar_div_fmpq(p.p,p.p,c,ctx);fmpq_clear(c);}
void jsonp(ostream&f,const P&p){f<<'[';fmpq_t c;fmpq_init(c);ulong ex[3];for(long i=0;i<fmpq_mpoly_length(p.p,ctx);i++){if(i)f<<',';fmpq_mpoly_get_term_coeff_fmpq(c,p.p,i,ctx);fmpq_mpoly_get_term_exp_ui(ex,p.p,i,ctx);char*s=fmpq_get_str(nullptr,10,c);f<<"[["<<ex[0]<<','<<ex[1]<<','<<ex[2]<<"],\""<<s<<"\"]";flint_free(s);}f<<']';fmpq_clear(c);assert(f.tellp()<30*1024*1024);}
P quotient(const P&a,const P&b){P q,c;assert(fmpq_mpoly_divides(q.p,a.p,b.p,ctx));fmpq_mpoly_mul(c.p,q.p,b.p,ctx);assert(fmpq_mpoly_equal(c.p,a.p,ctx));return q;}
int main(int argc,char**argv){assert(argc==3);auto start=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);
 {
  ifstream f(argv[1]);int n;f>>n;assert(n==22);vector<P>A(n);for(P&p:A)readp(f,p);P U;readp(f,U);string extra;assert(!(f>>extra));vector<P>chain;P g=A[0];monic(g);chain.push_back(g);
  for(int i=1;i<n;i++){P h;assert(fmpq_mpoly_gcd(h.p,g.p,A[i].p,ctx));monic(h);g=h;chain.push_back(g);assert(chrono::duration<double>(chrono::steady_clock::now()-start).count()<50);}
  vector<P>Q;for(const P&p:A)Q.push_back(quotient(p,g));P h=g;vector<P>removed,beforequot,Uquot;
  while(true){P t;assert(fmpq_mpoly_gcd(t.p,h.p,U.p,ctx));monic(t);if(fmpq_mpoly_is_fmpq(t.p,ctx))break;removed.push_back(t);h=quotient(h,t);beforequot.push_back(h);Uquot.push_back(quotient(U,t));assert(removed.size()<=15);}
  P coprime;assert(fmpq_mpoly_gcd(coprime.p,h.p,U.p,ctx));monic(coprime);assert(fmpq_mpoly_is_one(coprime.p,ctx));ofstream o(argv[2]);o<<"{\"gcd22\":";jsonp(o,g);o<<",\"gcd_chain22\":[";for(int i=0;i<n;i++){if(i)o<<',';jsonp(o,chain[i]);}o<<"],\"all22_literal_quotients\":[";for(int i=0;i<n;i++){if(i)o<<',';jsonp(o,Q[i]);}o<<"],\"localizer_removal_steps\":[";
  for(int i=0;i<(int)removed.size();i++){if(i)o<<',';o<<"{\"removed_gcd_with_U\":";jsonp(o,removed[i]);o<<",\"remaining_h_after\":";jsonp(o,beforequot[i]);o<<",\"U_quotient_by_removed\":";jsonp(o,Uquot[i]);o<<'}';}
  o<<"],\"remaining_h_coprime_U\":";jsonp(o,h);o<<",\"final_gcd_h_U\":";jsonp(o,coprime);o<<",\"all22_reconstructions_verified\":true,\"all_removed_factors_divide_U_exactly\":true,\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
 }
 fmpq_mpoly_ctx_clear(ctx);flint_cleanup();return 0;
}
