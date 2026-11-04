// Exact fixed4 norm/adjugate/Bezout and localizer-content removal. No GB/points.
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
 ifstream f(argv[1]);vector<vector<P>>N(4,vector<P>(4));for(auto&row:N)for(P&p:row)readp(f,p);L H=readd(f),g=readd(f),fsource=readd(f),qparameter=readd(f);vector<L>Q(4);for(L&p:Q)p=readd(f);P U;readp(f,U);long aclear;string D;f>>aclear>>D;assert(aclear==10);string extra;assert(!(f>>extra));L s=L(scalar(one(),D),aclear);
 // Bind the submitted cleared matrix to the actual monic division source.
 for(int j=0;j<4;j++){L lhs=L(shift(g.p,2,j),g.a),rhs=lmul(H,Q[j]);for(int i=0;i<4;i++)rhs=ladd(rhs,L(scalar(shift(N[i][j],2,i),"1/"+D),-aclear));assert(leq(lhs,rhs));}
 cerr<<"INPUT_IDENTITIES_PASS\n";P R=det(N);assert(!fmpq_mpoly_is_zero(R.p,ctx));vector<P>v(4);
 for(int j=0;j<4;j++){vector<vector<P>>B(3,vector<P>(3));for(int i=1;i<4;i++){int k=0;for(int c=0;c<4;c++)if(c!=j)B[i-1][k++]=N[i][c];}v[j]=det(B);if(j%2)v[j]=neg(v[j]);}
 for(int i=0;i<4;i++){P z;for(int j=0;j<4;j++)z=add(z,mul(N[i][j],v[j]));assert(eq(z,i==0?R:P()));}cerr<<"DET_ADJ_PASS\n";
 L V,Qsum;for(int j=0;j<4;j++){V=ladd(V,L(shift(v[j],2,j)));Qsum=ladd(Qsum,lmul(L(v[j]),Q[j]));}
 assert(leq(ladd(lmul(s,lmul(g,V)),lneg(lmul(s,lmul(H,Qsum)))),L(R)));
 L alpha=L(scalar(one(),"-9/11"),-4),beta=lscalar(L(qparameter.p,qparameter.a-5),"9/11");assert(leq(H,ladd(lmul(alpha,fsource),lmul(beta,g))));
 L uf=lneg(lmul(s,lmul(Qsum,alpha))),vg=ladd(lmul(s,V),lneg(lmul(s,lmul(Qsum,beta))));assert(leq(ladd(lmul(uf,fsource),lmul(vg,g)),L(R)));cerr<<"LITERAL_BEZOUT_PASS\n";
 P primitive;fmpq_t content;fmpq_init(content);fmpq_mpoly_content(content,R.p,ctx);fmpq_mpoly_scalar_div_fmpq(primitive.p,R.p,content,ctx);P h=primitive;vector<P>removed,after,Uquot;
 while(true){P t;assert(fmpq_mpoly_gcd(t.p,h.p,U.p,ctx));monic(t);if(fmpq_mpoly_is_fmpq(t.p,ctx))break;removed.push_back(t);h=quotient(h,t);after.push_back(h);Uquot.push_back(quotient(U,t));assert(removed.size()<=100);}
 P cop;assert(fmpq_mpoly_gcd(cop.p,h.p,U.p,ctx));monic(cop);assert(fmpq_mpoly_is_one(cop.p,ctx));
 ofstream o(argv[2]);o<<"{\"det_cleared_matrix_R\":";jp(o,R);o<<",\"adjugate_first_column\":[";for(int j=0;j<4;j++){if(j)o<<',';jp(o,v[j]);}o<<"],\"monic_V\":";jl(o,V);o<<",\"monic_Qsum\":";jl(o,Qsum);o<<",\"literal_membership_f_cofactor\":";jl(o,uf);o<<",\"literal_membership_g_cofactor\":";jl(o,vg);char*cs=fmpq_get_str(nullptr,10,content);o<<",\"R_QQ_content\":\""<<cs<<"\",\"R_primitive\":";flint_free(cs);jp(o,primitive);o<<",\"localizer_removal_steps\":[";for(int j=0;j<(int)removed.size();j++){if(j)o<<',';o<<"{\"removed_gcd_with_U\":";jp(o,removed[j]);o<<",\"remaining_after\":";jp(o,after[j]);o<<",\"U_quotient\":";jp(o,Uquot[j]);o<<'}';}o<<"],\"remaining_h_coprime_U\":";jp(o,h);o<<",\"final_gcd_h_U\":";jp(o,cop);o<<",\"all4adj_and_literal_membership_identities_checked\":true,\"peak_single_polynomial_support\":"<<peak<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";fmpq_clear(content);cerr<<"EXPORT_PASS\n";
}fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
