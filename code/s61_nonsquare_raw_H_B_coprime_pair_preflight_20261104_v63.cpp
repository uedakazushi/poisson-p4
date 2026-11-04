// Fixed22-minor pair-GCD cost preflight. No resultant/GB/factorization/points.
#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <cassert>
#include <chrono>
#include <fstream>
#include <string>
#include <vector>
using namespace std;fmpq_mpoly_ctx_t ctx;
struct P{fmpq_mpoly_t p;P(){fmpq_mpoly_init(p,ctx);}P(const P&v):P(){fmpq_mpoly_set(p,v.p,ctx);}P(P&&v):P(){fmpq_mpoly_swap(p,v.p,ctx);}~P(){fmpq_mpoly_clear(p,ctx);}};
void readp(istream&f,P&p){long n;f>>n;assert(f&&n>0);fmpq_t c;fmpq_init(c);ulong ex[3];for(long i=0;i<n;i++){string s;f>>s>>ex[0]>>ex[1]>>ex[2];assert(f&&!fmpq_set_str(c,s.c_str(),10));fmpq_canonicalise(c);fmpq_mpoly_push_term_fmpq_ui(p.p,c,ex,ctx);}fmpq_mpoly_sort_terms(p.p,ctx);fmpq_mpoly_combine_like_terms(p.p,ctx);assert(fmpq_mpoly_length(p.p,ctx)==n);fmpq_clear(c);}
void monic(P&p){fmpq_t c;fmpq_init(c);assert(!fmpq_mpoly_is_zero(p.p,ctx));fmpq_mpoly_get_term_coeff_fmpq(c,p.p,0,ctx);fmpq_mpoly_scalar_div_fmpq(p.p,p.p,c,ctx);fmpq_clear(c);}
void jp(ostream&f,const P&p){f<<'[';fmpq_t c;fmpq_init(c);ulong ex[3];for(long i=0;i<fmpq_mpoly_length(p.p,ctx);i++){if(i)f<<',';fmpq_mpoly_get_term_coeff_fmpq(c,p.p,i,ctx);fmpq_mpoly_get_term_exp_ui(ex,p.p,i,ctx);char*s=fmpq_get_str(nullptr,10,c);f<<"[["<<ex[0]<<','<<ex[1]<<','<<ex[2]<<"],\""<<s<<"\"]";flint_free(s);}f<<']';fmpq_clear(c);}
int main(int argc,char**argv){assert(argc==3);auto start=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);{
 ifstream f(argv[1]);int n,k;f>>n;assert(n==22);vector<P>A(n);for(auto&p:A)readp(f,p);f>>k;assert(0<k&&k<=231);vector<pair<int,int>>pairs;for(int x=0;x<k;x++){int i,j;f>>i>>j;assert(0<=i&&i<j&&j<n);pairs.emplace_back(i,j);}string extra;assert(!(f>>extra));
 ofstream o(argv[2]);o<<"{\"tried_ordered_pairs\":[";bool found=false;int cnt=0;
 for(auto[i,j]:pairs){P g,q0,q1,t;assert(fmpq_mpoly_gcd(g.p,A[i].p,A[j].p,ctx));monic(g);assert(fmpq_mpoly_divides(q0.p,A[i].p,g.p,ctx));assert(fmpq_mpoly_divides(q1.p,A[j].p,g.p,ctx));fmpq_mpoly_mul(t.p,g.p,q0.p,ctx);assert(fmpq_mpoly_equal(t.p,A[i].p,ctx));fmpq_mpoly_mul(t.p,g.p,q1.p,ctx);assert(fmpq_mpoly_equal(t.p,A[j].p,ctx));if(cnt++)o<<',';o<<"{\"pair\":["<<i<<','<<j<<"],\"gcd\":";jp(o,g);o<<",\"literal_quotients\":[";jp(o,q0);o<<',';jp(o,q1);o<<"]}";assert(o.tellp()<30*1024*1024);if(fmpq_mpoly_is_one(g.p,ctx)){found=true;break;}assert(chrono::duration<double>(chrono::steady_clock::now()-start).count()<50);}
 o<<"],\"found_coprime_pair\":"<<(found?"true":"false")<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
}fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
