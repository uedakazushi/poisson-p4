// Atomic-denominator preflight ONLY: no scalar DAG or quartic coefficients.
#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <array>
using namespace std;fmpq_mpoly_ctx_t ctx;
struct P{fmpq_mpoly_t p;P(){fmpq_mpoly_init(p,ctx);}P(const P&v):P(){fmpq_mpoly_set(p,v.p,ctx);}P(P&&v):P(){fmpq_mpoly_swap(p,v.p,ctx);}P&operator=(P v){fmpq_mpoly_swap(p,v.p,ctx);return*this;}~P(){fmpq_mpoly_clear(p,ctx);}};
void readp(istream&f,P&q){
 long n;f>>n;assert(f&&n>=0);fmpq_t c;fmpq_init(c);ulong ex[3];
 for(long i=0;i<n;i++){string z;f>>z>>ex[0]>>ex[1]>>ex[2];assert(f&&!fmpq_set_str(c,z.c_str(),10));fmpq_canonicalise(c);assert(!fmpq_is_zero(c));fmpq_mpoly_push_term_fmpq_ui(q.p,c,ex,ctx);}
 fmpq_mpoly_sort_terms(q.p,ctx);fmpq_mpoly_combine_like_terms(q.p,ctx);assert(fmpq_mpoly_length(q.p,ctx)==n);fmpq_clear(c);
}
int main(int argc,char**argv){
 assert(argc==3);auto began=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);bool passed=true;
 {
  ifstream f(argv[1]);ofstream out(argv[2]);long nd,na,nk;f>>nd>>na>>nk;assert(nd==228&&na==2378&&nk==6);
  vector<P>known(nk),den(nd),reconstructed(nd);for(auto&p:known)readp(f,p);for(auto&p:den)readp(f,p);
  vector<array<long,6>>exponents(nd);for(auto&v:exponents)v.fill(0);vector<string>constants(nd);long divisions=0;fmpq_t c;fmpq_init(c);
  for(int i=0;i<nd;i++){
   P p=den[i],q;assert(!fmpq_mpoly_is_zero(p.p,ctx));
   for(int k=0;k<6;k++)while(fmpq_mpoly_divides(q.p,p.p,known[k].p,ctx)){p=q;exponents[i][k]++;divisions++;assert(exponents[i][k]<=1000);}
   if(!fmpq_mpoly_is_fmpq(p.p,ctx)){cerr<<"UNKNOWN_DENOMINATOR_FACTOR "<<i<<'\n';passed=false;break;}
   fmpq_mpoly_get_fmpq(c,p.p,ctx);assert(!fmpq_is_zero(c));char*str=fmpq_get_str(nullptr,10,c);constants[i]=str;flint_free(str);reconstructed[i]=p;
   for(int k=0;k<6;k++)for(int n=0;n<exponents[i][k];n++){P z;fmpq_mpoly_mul(z.p,reconstructed[i].p,known[k].p,ctx);reconstructed[i]=z;}
   assert(fmpq_mpoly_equal(reconstructed[i].p,den[i].p,ctx));
  }
  long crossproducts=0;long max_product_terms=0;
  if(passed)for(int j=0;j<na;j++){
   long di;f>>di;assert(0<=di&&di<nd);P n,l,r;readp(f,n);fmpq_mpoly_mul(l.p,n.p,den[di].p,ctx);fmpq_mpoly_mul(r.p,n.p,reconstructed[di].p,ctx);
   assert(fmpq_mpoly_equal(l.p,r.p,ctx));crossproducts++;max_product_terms=max(max_product_terms,fmpq_mpoly_length(l.p,ctx));
  }
  if(passed){string extra;assert(!(f>>extra));}
  out<<"{\"all_passed\":"<<(passed?"true":"false")<<",\"denominator_reconstructions\":"<<(passed?nd:0)<<",\"atomic_fraction_crossproducts\":"<<crossproducts<<",\"exact_named_divisions\":"<<divisions<<",\"max_fraction_crossproduct_terms\":"<<max_product_terms<<",\"factors\":[";
  if(passed)for(int i=0;i<nd;i++){if(i)out<<',';out<<"{\"constant\":\""<<constants[i]<<"\",\"a_d_F_B_C_J_powers\":[";for(int k=0;k<6;k++){if(k)out<<',';out<<exponents[i][k];}out<<"]}";}
  out<<"],\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";fmpq_clear(c);
 }
 fmpq_mpoly_ctx_clear(ctx);flint_cleanup();return passed?0:3;
}
