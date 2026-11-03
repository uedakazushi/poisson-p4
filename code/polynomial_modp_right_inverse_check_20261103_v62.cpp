// Independent FLINT exact polynomial check of M*L=g*identity over Fp[a,e,P].
// CLI: source-native right-inverse-native prime. Third variable d denotes P.
#include <flint/nmod_mpoly.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <chrono>
using namespace std;nmod_mpoly_ctx_t ctx;ulong prime;
struct Poly{nmod_mpoly_t p;Poly(){nmod_mpoly_init(p,ctx);}Poly(const Poly&q):Poly(){nmod_mpoly_set(p,q.p,ctx);}Poly(Poly&&q):Poly(){nmod_mpoly_swap(p,q.p,ctx);}Poly&operator=(Poly q){nmod_mpoly_swap(p,q.p,ctx);return *this;}~Poly(){nmod_mpoly_clear(p,ctx);}};
vector<vector<Poly>> read(const char*name){ifstream f(name);long n,m,k;f>>n>>m>>k;assert(f&&n>0&&m>0&&k==0);vector<vector<Poly>>A(n,vector<Poly>(m));ulong ex[3];for(auto&row:A)for(auto&p:row){long count;f>>count;assert(f&&count>=0);for(long j=0;j<count;j++){long long c;f>>c>>ex[0]>>ex[1]>>ex[2];assert(f);long long v=c%(long long)prime;if(v<0)v+=prime;assert(v);nmod_mpoly_push_term_ui_ui(p.p,(ulong)v,ex,ctx);}nmod_mpoly_sort_terms(p.p,ctx);nmod_mpoly_combine_like_terms(p.p,ctx);assert(nmod_mpoly_length(p.p,ctx)==count);}string extra;assert(!(f>>extra));return A;}
int main(int argc,char**argv){assert(argc==4);prime=stoul(argv[3]);nmod_mpoly_ctx_init(ctx,3,ORD_LEX,prime);{auto start=chrono::steady_clock::now();auto M=read(argv[1]);auto L=read(argv[2]);long n=M.size(),m=M[0].size();assert((long)L.size()==m&&L[0].size()==(size_t)n);Poly g;const char*vars[]={"a","e","d"};assert(!nmod_mpoly_set_str_pretty(g.p,"a*d^3*(100*a+18*e^2-9*e-2)*(5*d-6*e-4)^2",vars,ctx));long checks=0;for(long i=0;i<n;i++){for(long j=0;j<n;j++){Poly sum,tmp;for(long k=0;k<m;k++){nmod_mpoly_mul(tmp.p,M[i][k].p,L[k][j].p,ctx);nmod_mpoly_add(sum.p,sum.p,tmp.p,ctx);}assert(i==j?nmod_mpoly_equal(sum.p,g.p,ctx):nmod_mpoly_is_zero(sum.p,ctx));checks++;}cerr<<"EXACT_FP_ROW_CHECKED "<<i<<" SEC "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';}cout<<"EXACT_ALL_FP_RIGHT_INVERSE_IDENTITIES_PASS "<<checks<<'\n';}nmod_mpoly_ctx_clear(ctx);flint_cleanup();}
