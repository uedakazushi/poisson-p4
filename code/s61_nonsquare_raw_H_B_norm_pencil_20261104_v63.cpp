// Exact norm pencil mixed3 coefficients and commoncurve discrimination ONLY.
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
int main(int argc,char**argv){assert(argc==3);start=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);{
 ifstream f(argv[1]);vector<vector<P>>A(4,vector<P>(4)),B(4,vector<P>(4));for(auto&row:A)for(P&p:row)readp(f,p);for(auto&row:B)for(P&p:row)readp(f,p);P U,G,oldR0,oldR4;readp(f,U);readp(f,G);readp(f,oldR0);readp(f,oldR4);string extra;assert(!(f>>extra));
 vector<P>C(5);vector<int>perm={0,1,2,3};do{int inv=0;for(int i=0;i<4;i++)for(int j=i+1;j<4;j++)inv+=perm[i]>perm[j];vector<P>v(1);v[0]=one();for(int i=0;i<4;i++){vector<P>nv(v.size()+1);for(int j=0;j<(int)v.size();j++){nv[j]=add(nv[j],mul(v[j],A[i][perm[i]]));nv[j+1]=add(nv[j+1],mul(v[j],B[i][perm[i]]));}v=move(nv);}for(int j=0;j<5;j++)C[j]=add(C[j],inv%2?neg(v[j]):v[j]);}while(next_permutation(perm.begin(),perm.end()));
 P ratio4=scalar(shift(one(),0,4),"625000000000000"); // (5000*a)^4
 assert(eq(C[0],mul(ratio4,oldR0)));assert(eq(C[4],oldR4));cerr<<"ALL5_COEFFICIENTS_AND_ENDPOINTS_PASS\n";
 vector<P>prim(3),h(3);vector<vector<P>>removed(3),after(3),uq(3);vector<string>contents;
 for(int j=0;j<3;j++){if(fmpq_mpoly_is_zero(C[j+1].p,ctx)){contents.push_back("0");continue;}fmpq_t con;fmpq_init(con);fmpq_mpoly_content(con,C[j+1].p,ctx);char*cs=fmpq_get_str(nullptr,10,con);contents.push_back(cs);flint_free(cs);fmpq_mpoly_scalar_div_fmpq(prim[j].p,C[j+1].p,con,ctx);fmpq_clear(con);h[j]=prim[j];while(true){P t;assert(fmpq_mpoly_gcd(t.p,h[j].p,U.p,ctx));monic(t);if(fmpq_mpoly_is_fmpq(t.p,ctx))break;removed[j].push_back(t);h[j]=quotient(h[j],t);after[j].push_back(h[j]);uq[j].push_back(quotient(U,t));assert(removed[j].size()<=100);}P cop;assert(fmpq_mpoly_gcd(cop.p,h[j].p,U.p,ctx));monic(cop);assert(fmpq_mpoly_is_one(cop.p,ctx));}
 vector<P>chain,bq,sq;P common=G;monic(common);for(int j=0;j<3;j++){P cg;assert(fmpq_mpoly_gcd(cg.p,common.p,h[j].p,ctx));monic(cg);bq.push_back(quotient(common,cg));sq.push_back(quotient(h[j],cg));common=cg;chain.push_back(common);}cerr<<"UNIT_STRIP_AND_CURVE_GCD_PASS\n";
 ofstream o(argv[2]);o<<"{\"pencil_coefficients_lambda0_to4\":[";for(int j=0;j<5;j++){if(j)o<<',';jp(o,C[j]);}o<<"],\"mixed3_records\":[";for(int j=0;j<3;j++){if(j)o<<',';o<<"{\"lambda_degree\":"<<j+1<<",\"QQ_content\":\""<<contents[j]<<"\",\"primitive\":";jp(o,prim[j]);o<<",\"localizer_removal_steps\":[";for(int k=0;k<(int)removed[j].size();k++){if(k)o<<',';o<<"{\"removed_gcd_with_U\":";jp(o,removed[j][k]);o<<",\"remaining_after\":";jp(o,after[j][k]);o<<",\"U_quotient\":";jp(o,uq[j][k]);o<<'}';}o<<"],\"remaining_h_coprime_U\":";jp(o,h[j]);o<<",\"curve_gcd_after_this\":";jp(o,chain[j]);o<<",\"previous_curve_quotient_by_gcd\":";jp(o,bq[j]);o<<",\"mixed_quotient_by_gcd\":";jp(o,sq[j]);o<<'}';}o<<"],\"final_common_curve_gcd\":";jp(o,common);o<<",\"original_G_quotient_by_final\":";jp(o,quotient(G,common));o<<",\"mixed3_quotients_by_final\":[";for(int j=0;j<3;j++){if(j)o<<',';jp(o,quotient(h[j],common));}o<<"],\"endpoint_scaling_exact\":true,\"peak_single_polynomial_support\":"<<peak<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";cerr<<"EXPORT_PASS\n";
}fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
