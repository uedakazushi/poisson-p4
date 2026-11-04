// Independent fixed Laurent-polynomial identities, last-row Laplace and
// row-multilinearity norm pencil. No author arithmetic code is imported.
#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
fmpq_mpoly_ctx_t ctx;long peak=0;auto began=chrono::steady_clock::now();
struct P{fmpq_mpoly_t v;P(){fmpq_mpoly_init(v,ctx);}P(const P&p):P(){fmpq_mpoly_set(v,p.v,ctx);}P(P&&p):P(){fmpq_mpoly_swap(v,p.v,ctx);}P&operator=(P p){fmpq_mpoly_swap(v,p.v,ctx);return *this;}~P(){fmpq_mpoly_clear(v,ctx);}};
void guard(const P&p){assert(fmpq_mpoly_is_canonical(p.v,ctx));peak=max(peak,fmpq_mpoly_length(p.v,ctx));assert(peak<200000);assert(chrono::duration<double>(chrono::steady_clock::now()-began).count()<50);}
P shifted(const P&p,long n){assert(n>=0);P r;fmpq_t c;fmpq_init(c);ulong ex[3];for(long i=0;i<fmpq_mpoly_length(p.v,ctx);i++){fmpq_mpoly_get_term_coeff_fmpq(c,p.v,i,ctx);fmpq_mpoly_get_term_exp_ui(ex,p.v,i,ctx);ex[0]+=n;fmpq_mpoly_push_term_fmpq_ui(r.v,c,ex,ctx);}fmpq_mpoly_sort_terms(r.v,ctx);fmpq_mpoly_combine_like_terms(r.v,ctx);fmpq_clear(c);return r;}
struct L{P p;long a=0;L(){}L(P q,long z=0):p(std::move(q)),a(z){normalize();}void normalize(){if(fmpq_mpoly_is_zero(p.v,ctx)){a=0;return;}ulong ex[3];long m=10000000;for(long i=0;i<fmpq_mpoly_length(p.v,ctx);i++){fmpq_mpoly_get_term_exp_ui(ex,p.v,i,ctx);m=min(m,(long)ex[0]);}if(m){P r;fmpq_t c;fmpq_init(c);for(long i=0;i<fmpq_mpoly_length(p.v,ctx);i++){fmpq_mpoly_get_term_coeff_fmpq(c,p.v,i,ctx);fmpq_mpoly_get_term_exp_ui(ex,p.v,i,ctx);ex[0]-=m;fmpq_mpoly_push_term_fmpq_ui(r.v,c,ex,ctx);}fmpq_mpoly_sort_terms(r.v,ctx);fmpq_mpoly_combine_like_terms(r.v,ctx);fmpq_clear(c);p=std::move(r);a+=m;}guard(p);}};
L addLaurent(const L&x,const L&y,bool minus=false){long z=min(x.a,y.a);P p=shifted(x.p,x.a-z),q=shifted(y.p,y.a-z),r;if(minus)fmpq_mpoly_sub(r.v,p.v,q.v,ctx);else fmpq_mpoly_add(r.v,p.v,q.v,ctx);return L(std::move(r),z);}
L times(const L&x,const L&y){P p;fmpq_mpoly_mul(p.v,x.p.v,y.p.v,ctx);return L(std::move(p),x.a+y.a);}
L negative(const L&x){P p;fmpq_mpoly_neg(p.v,x.p.v,ctx);return L(std::move(p),x.a);}
bool equal(const L&x,const L&y){return x.a==y.a&&fmpq_mpoly_equal(x.p.v,y.p.v,ctx);}
L one(){P p;fmpq_mpoly_one(p.v,ctx);return L(std::move(p));}
L determinant(const vector<vector<L>>&A){int n=A.size();if(!n)return one();if(n==1)return A[0][0];L r;for(int j=0;j<n;j++){vector<vector<L>>B(n-1);for(int i=0;i<n-1;i++)for(int k=0;k<n;k++)if(k!=j)B[i].push_back(A[i][k]);L t=times(A[n-1][j],determinant(B));r=addLaurent(r,t,((n-1+j)&1));}return r;}
void monic(P&p){if(fmpq_mpoly_is_zero(p.v,ctx))return;fmpq_t c;fmpq_init(c);fmpq_mpoly_get_term_coeff_fmpq(c,p.v,0,ctx);fmpq_mpoly_scalar_div_fmpq(p.v,p.v,c,ctx);fmpq_clear(c);}
int main(int argc,char**argv){assert(argc==3);began=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,3,ORD_DEGREVLEX);long eqs=0,gcds=0,dets=0;{
 ifstream f(argv[1]);long count;f>>count;vector<L>v;v.reserve(count);for(long i=0;i<count;i++){char op;f>>op;if(op=='P'){long off,n;f>>off>>n;P p;fmpq_t c;fmpq_init(c);ulong ex[3];for(long j=0;j<n;j++){string s;f>>s>>ex[0]>>ex[1]>>ex[2];assert(!fmpq_set_str(c,s.c_str(),10));fmpq_canonicalise(c);fmpq_mpoly_push_term_fmpq_ui(p.v,c,ex,ctx);}fmpq_mpoly_sort_terms(p.v,ctx);fmpq_mpoly_combine_like_terms(p.v,ctx);assert(fmpq_mpoly_length(p.v,ctx)==n);fmpq_clear(c);v.emplace_back(std::move(p),off);
 }else if(op=='D'){int n;f>>n;vector<vector<L>>A(n);for(int j=0;j<n;j++)for(int k=0;k<n;k++){long t;f>>t;assert(t<i);A[j].push_back(v[t]);}v.push_back(determinant(A));dets++;
 }else{long x,y;f>>x>>y;assert(x<i&&y<i);v.push_back(op=='M'?times(v[x],v[y]):addLaurent(v[x],v[y],op=='S'));}}
 long assertions;f>>assertions;for(long i=0;i<assertions;i++){char op;long x,y;f>>op>>x>>y;assert(x<count&&y<count);if(op=='E'){if(!equal(v[x],v[y]))cerr<<"FAILED_EQ "<<i<<" NODES "<<x<<' '<<y<<" AOFFSET "<<v[x].a<<' '<<v[y].a<<" TERMS "<<fmpq_mpoly_length(v[x].p.v,ctx)<<' '<<fmpq_mpoly_length(v[y].p.v,ctx)<<'\n';assert(equal(v[x],v[y]));eqs++;}else{long z;f>>z;assert(op=='G'&&z<count&&v[x].a>=0&&v[y].a>=0&&v[z].a>=0);P p=shifted(v[x].p,v[x].a),q=shifted(v[y].p,v[y].a),g,expected=shifted(v[z].p,v[z].a);assert(fmpq_mpoly_gcd(g.v,p.v,q.v,ctx));monic(g);monic(expected);assert(fmpq_mpoly_equal(g.v,expected.v,ctx));gcds++;}guard(v[x].p);}string extra;assert(!(f>>extra));
 ofstream o(argv[2]);o<<"{\"all_passed\":true,\"exact_identities\":"<<eqs<<",\"independent_GCD_checks\":"<<gcds<<",\"last_row_Laplace_determinants\":"<<dets<<",\"peak_polynomial_terms\":"<<peak<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
