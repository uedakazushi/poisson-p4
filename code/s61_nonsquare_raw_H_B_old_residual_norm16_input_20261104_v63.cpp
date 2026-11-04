// Oldfinite residualmonic-e16 multiplicationmatrix INPUT ONLY. No determinant/resultant.
#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;fmpq_mpoly_ctx_t ctx;long peak=0;chrono::steady_clock::time_point start;
void tick(){assert(chrono::duration<double>(chrono::steady_clock::now()-start).count()<50);}
struct P{fmpq_mpoly_t p;P(){fmpq_mpoly_init(p,ctx);}P(const P&v):P(){fmpq_mpoly_set(p,v.p,ctx);}P(P&&v):P(){fmpq_mpoly_swap(p,v.p,ctx);}P&operator=(P v){fmpq_mpoly_swap(p,v.p,ctx);return*this;}~P(){fmpq_mpoly_clear(p,ctx);}};
void guard(const P&p){long n=fmpq_mpoly_length(p.p,ctx);peak=max(peak,n);assert(n<=150000);tick();}
P one(){P p;fmpq_mpoly_one(p.p,ctx);return p;}
P add(const P&a,const P&b){P r;fmpq_mpoly_add(r.p,a.p,b.p,ctx);guard(r);return r;}
P neg(const P&a){P r;fmpq_mpoly_neg(r.p,a.p,ctx);return r;}
P mul(const P&a,const P&b){P r;fmpq_mpoly_mul(r.p,a.p,b.p,ctx);guard(r);return r;}
bool eq(const P&a,const P&b){return fmpq_mpoly_equal(a.p,b.p,ctx);}
P scalar(const P&a,const string&s){P r;fmpq_t c;fmpq_init(c);assert(!fmpq_set_str(c,s.c_str(),10));fmpq_canonicalise(c);fmpq_mpoly_scalar_mul_fmpq(r.p,a.p,c,ctx);fmpq_clear(c);return r;}
P shift(const P&a,int var,long n){assert(n>=0);P r;ulong ex[3];fmpq_t c;fmpq_init(c);for(long i=0;i<fmpq_mpoly_length(a.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(ex,a.p,i,ctx);ex[var]+=n;fmpq_mpoly_get_term_coeff_fmpq(c,a.p,i,ctx);fmpq_mpoly_push_term_fmpq_ui(r.p,c,ex,ctx);}fmpq_mpoly_sort_terms(r.p,ctx);fmpq_clear(c);return r;}
P quotient(const P&a,const P&b){P q;assert(fmpq_mpoly_divides(q.p,a.p,b.p,ctx));assert(eq(mul(q,b),a));return q;}
void monic(P&p){assert(!fmpq_mpoly_is_zero(p.p,ctx));fmpq_t c;fmpq_init(c);fmpq_mpoly_get_term_coeff_fmpq(c,p.p,0,ctx);fmpq_mpoly_scalar_div_fmpq(p.p,p.p,c,ctx);fmpq_clear(c);}
void readp(istream&f,P&p){long n;f>>n;assert(f&&n>=0);fmpq_t c;fmpq_init(c);ulong ex[3];for(long i=0;i<n;i++){string s;f>>s>>ex[0]>>ex[1]>>ex[2];assert(f&&!fmpq_set_str(c,s.c_str(),10));fmpq_canonicalise(c);fmpq_mpoly_push_term_fmpq_ui(p.p,c,ex,ctx);}fmpq_mpoly_sort_terms(p.p,ctx);fmpq_mpoly_combine_like_terms(p.p,ctx);assert(fmpq_mpoly_length(p.p,ctx)==n);fmpq_clear(c);guard(p);}
struct L{P p;long a=0;L(){}L(P q,long z=0):p(std::move(q)),a(z){norm();}void norm(){if(fmpq_mpoly_is_zero(p.p,ctx)){a=0;return;}ulong ex[3];long m=1000000000;for(long i=0;i<fmpq_mpoly_length(p.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(ex,p.p,i,ctx);m=min(m,(long)ex[0]);}if(m){p=quotient(p,shift(one(),0,m));a+=m;}}};
L ladd(const L&a,const L&b){long z=min(a.a,b.a);return L(add(shift(a.p,0,a.a-z),shift(b.p,0,b.a-z)),z);}
L lneg(const L&a){return L(neg(a.p),a.a);}
L lmul(const L&a,const L&b){return L(mul(a.p,b.p),a.a+b.a);}
L lscalar(const L&a,const string&s){return L(scalar(a.p,s),a.a);}
bool leq(const L&a,const L&b){return a.a==b.a&&eq(a.p,b.p);}
L readd(istream&f){long sh;f>>sh;P p;readp(f,p);return L(p,sh);}
P det(const vector<vector<P>>&A){int n=A.size();vector<int>p(n);for(int i=0;i<n;i++)p[i]=i;P r;do{int inv=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)inv+=p[i]>p[j];P t=one();for(int i=0;i<n;i++)t=mul(t,A[i][p[i]]);r=add(r,inv%2?neg(t):t);}while(next_permutation(p.begin(),p.end()));return r;}
void jp(ostream&f,const P&p,long a=0){f<<'[';fmpq_t c;fmpq_init(c);ulong ex[3];for(long i=0;i<fmpq_mpoly_length(p.p,ctx);i++){if(i)f<<',';fmpq_mpoly_get_term_coeff_fmpq(c,p.p,i,ctx);fmpq_mpoly_get_term_exp_ui(ex,p.p,i,ctx);char*s=fmpq_get_str(nullptr,10,c);f<<"[["<<(long(ex[0])+a)<<','<<ex[1]<<','<<ex[2]<<"],\""<<s<<"\"]";flint_free(s);}f<<']';fmpq_clear(c);assert(f.tellp()<31*1024*1024);}
void jl(ostream&f,const L&a){jp(f,a.p,a.a);}
P ecoeff(const P&p,long k){P r;fmpq_t c;fmpq_init(c);ulong ex[3];for(long i=0;i<fmpq_mpoly_length(p.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(ex,p.p,i,ctx);if((long)ex[0]!=k)continue;ex[0]=0;fmpq_mpoly_get_term_coeff_fmpq(c,p.p,i,ctx);fmpq_mpoly_push_term_fmpq_ui(r.p,c,ex,ctx);}fmpq_mpoly_sort_terms(r.p,ctx);fmpq_clear(c);return r;}
int main(int argc,char**argv){assert(argc==3);start=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);{
 ifstream f(argv[1]);P H,K;readp(f,H);readp(f,K);string extra;assert(!(f>>extra));assert(fmpq_mpoly_degree_si(H.p,0,ctx)==16&&eq(ecoeff(H,16),one()));vector<P>quot(16),rem(16);vector<vector<P>>M(16,vector<P>(16));for(int j=0;j<16;j++){P src=shift(K,0,j);fmpq_mpoly_divrem(quot[j].p,rem[j].p,src.p,H.p,ctx);assert(eq(add(mul(quot[j],H),rem[j]),src));assert(fmpq_mpoly_degree_si(rem[j].p,0,ctx)<16);P rebuilt;for(int i=0;i<16;i++){M[i][j]=ecoeff(rem[j],i);rebuilt=add(rebuilt,shift(M[i][j],0,i));}assert(eq(rebuilt,rem[j]));guard(quot[j]);guard(rem[j]);}
 ofstream o(argv[2]);o<<"{\"native_variable_order\":[\"e\",\"a\",\"d_unused\"],\"column_quotients\":[";for(int j=0;j<16;j++){if(j)o<<',';jp(o,quot[j]);}o<<"],\"column_remainders\":[";for(int j=0;j<16;j++){if(j)o<<',';jp(o,rem[j]);}o<<"],\"matrix_rows\":[";for(int i=0;i<16;i++){if(i)o<<',';o<<'[';for(int j=0;j<16;j++){if(j)o<<',';jp(o,M[i][j]);}o<<']';}o<<"],\"all16_full_monic_division_identities_checked\":true,\"determinant_jobs\":0,\"peak_single_polynomial_support\":"<<peak<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
}fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
