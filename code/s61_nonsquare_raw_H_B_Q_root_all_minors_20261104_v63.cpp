// ALL22 fixedsource minors at rationalroot overirreducibleQ; no GB/points.
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

P coeff(const P&p,int var,long k){P r;ulong ex[3];fmpq_t c;fmpq_init(c);for(long i=0;i<fmpq_mpoly_length(p.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(ex,p.p,i,ctx);if((long)ex[var]!=k)continue;ex[var]=0;fmpq_mpoly_get_term_coeff_fmpq(c,p.p,i,ctx);fmpq_mpoly_push_term_fmpq_ui(r.p,c,ex,ctx);}fmpq_mpoly_sort_terms(r.p,ctx);fmpq_clear(c);return r;}
long degree(const P&p,int v){return fmpq_mpoly_degree_si(p.p,v,ctx);}
P evaluate(const P&p,const vector<P>&N,const vector<P>&D){long m=degree(p,2);assert(m>=0&&m<(long)N.size());P r;for(long j=0;j<=m;j++)r=add(r,mul(coeff(p,2,j),mul(N[j],D[m-j])));return r;}
struct Reduction{P value,quot,rem,norm,conj,normquot;};
Reduction reduce(const P&value,const P&Q){Reduction r;r.value=value;fmpq_mpoly_divrem(r.quot.p,r.rem.p,value.p,Q.p,ctx);assert(eq(add(mul(r.quot,Q),r.rem),value));assert(degree(r.rem,0)<=1);P u=coeff(r.rem,0,0),v=coeff(r.rem,0,1),q1=coeff(Q,0,1),q0=coeff(Q,0,0);r.norm=add(add(mul(u,u),neg(mul(q1,mul(u,v)))),mul(q0,mul(v,v)));r.conj=add(add(u,neg(mul(q1,v))),neg(shift(v,0,1)));P difference=add(mul(r.rem,r.conj),neg(r.norm));r.normquot=quotient(difference,Q);if(!fmpq_mpoly_is_zero(r.rem.p,ctx))assert(!fmpq_mpoly_is_zero(r.norm.p,ctx));return r;}
void jr(ostream&o,const Reduction&r){o<<"{\"cleared_evaluation\":";jp(o,r.value);o<<",\"Q_quotient\":";jp(o,r.quot);o<<",\"remainder\":";jp(o,r.rem);o<<",\"remainder_quadratic_norm_e\":";jp(o,r.norm);o<<",\"remainder_conjugate\":";jp(o,r.conj);o<<",\"norm_membership_Q_quotient\":";jp(o,r.normquot);o<<'}';}
int main(int argc,char**argv){assert(argc==3);start=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);{
 ifstream f(argv[1]);P Q,Numer,Den,U,Hclear;readp(f,Q);readp(f,Numer);readp(f,Den);readp(f,U);readp(f,Hclear);int n;f>>n;assert(n==22);vector<P>A(n);for(auto&p:A)readp(f,p);string extra;assert(!(f>>extra));assert(degree(Q,0)==2&&eq(coeff(Q,0,2),one()));vector<P>NP(8),DP(8);NP[0]=one();DP[0]=one();for(int i=1;i<8;i++){NP[i]=mul(NP[i-1],neg(Numer));DP[i]=mul(DP[i-1],Den);}
 vector<Reduction>R;vector<int>nz;for(int i=0;i<22;i++){R.push_back(reduce(evaluate(A[i],NP,DP),Q));if(!fmpq_mpoly_is_zero(R.back().rem.p,ctx))nz.push_back(i);tick();}
 Reduction RU=reduce(evaluate(U,NP,DP),Q),RH=reduce(evaluate(Hclear,NP,DP),Q);assert(!fmpq_mpoly_is_zero(RU.rem.p,ctx));assert(fmpq_mpoly_is_zero(RH.rem.p,ctx));
 ofstream o(argv[2]);o<<"{\"all22_evaluations\":[";for(int i=0;i<22;i++){if(i)o<<',';jr(o,R[i]);}o<<"],\"nonzero_remainder_indices\":[";for(int i=0;i<(int)nz.size();i++){if(i)o<<',';o<<nz[i];}o<<"],\"localizer_U_evaluation\":";jr(o,RU);o<<",\"Hclear_evaluation\":";jr(o,RH);o<<",\"all_reductions_and_norm_identities_checked\":true,\"peak_single_polynomial_support\":"<<peak<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
}fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
