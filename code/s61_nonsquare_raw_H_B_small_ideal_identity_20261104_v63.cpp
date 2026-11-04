// Independent exact QQ4 polynomial matrix identity A*L=B; no Groebner computation.
#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;fmpq_mpoly_ctx_t ctx;
struct P{fmpq_mpoly_t p;P(){fmpq_mpoly_init(p,ctx);}P(const P&v):P(){fmpq_mpoly_set(p,v.p,ctx);}P(P&&v):P(){fmpq_mpoly_swap(p,v.p,ctx);}~P(){fmpq_mpoly_clear(p,ctx);}};
void readp(istream&f,P&p){long n;f>>n;assert(f&&n>=0);fmpq_t c;fmpq_init(c);ulong ex[4];for(long i=0;i<n;i++){string s;f>>s>>ex[0]>>ex[1]>>ex[2]>>ex[3];assert(f&&!fmpq_set_str(c,s.c_str(),10));fmpq_canonicalise(c);fmpq_mpoly_push_term_fmpq_ui(p.p,c,ex,ctx);}fmpq_mpoly_sort_terms(p.p,ctx);fmpq_mpoly_combine_like_terms(p.p,ctx);assert(fmpq_mpoly_length(p.p,ctx)==n);fmpq_clear(c);}
int main(int argc,char**argv){assert(argc==5);fmpq_mpoly_ctx_init(ctx,4,ORD_LEX);bool pass=true;long active=0,maximum=0;
 {
  ifstream af(argv[1]),lf(argv[2]),bf(argv[3]);long ar,ac,lr,lc,br,bc;af>>ar>>ac;lf>>lr>>lc;bf>>br>>bc;assert(ar==1&&ac==23&&lr==23&&br==1&&bc==lc);vector<P>A(ac),L(lr*lc),B(bc);for(P&p:A)readp(af,p);for(P&p:L)readp(lf,p);for(P&p:B)readp(bf,p);string extra;assert(!(af>>extra)&&!(lf>>extra)&&!(bf>>extra));
  for(long j=0;j<lc;j++){P sum,term;for(long i=0;i<23;i++)if(!fmpq_mpoly_is_zero(L[i*lc+j].p,ctx)){active++;fmpq_mpoly_mul(term.p,A[i].p,L[i*lc+j].p,ctx);fmpq_mpoly_add(sum.p,sum.p,term.p,ctx);maximum=max(maximum,fmpq_mpoly_length(sum.p,ctx));assert(maximum<=150000);}assert(fmpq_mpoly_equal(sum.p,B[j].p,ctx));}
  ofstream o(argv[4]);o<<"{\"all_passed\":true,\"source_generator_count\":23,\"verified_output_columns\":"<<lc<<",\"nonzero_lift_coefficients\":"<<active<<",\"max_intermediate_terms\":"<<maximum<<"}\n";
 }
 fmpq_mpoly_ctx_clear(ctx);flint_cleanup();return pass?0:3;
}
