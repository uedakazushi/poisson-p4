// Exactoldresidual16 polynomialmatrix det+onlyone RHS solve; no GB/points.
#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_mat.h>
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
P shift(const P&a,int var,long n){assert(n>=0);P r;ulong ex[3];fmpq_t c;fmpq_init(c);for(long i=0;i<fmpq_mpoly_length(a.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(ex,a.p,i,ctx);ex[var]+=n;fmpq_mpoly_get_term_coeff_fmpq(c,a.p,i,ctx);fmpq_mpoly_push_term_fmpq_ui(r.p,c,ex,ctx);}fmpq_mpoly_sort_terms(r.p,ctx);fmpq_mpoly_combine_like_terms(r.p,ctx);assert(fmpq_mpoly_is_canonical(r.p,ctx));fmpq_clear(c);return r;}
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
P fromZ(const fmpz_poly_t z){P p;fmpq_t c;fmpq_init(c);ulong ex[3]={0,0,0};fmpz_t v;fmpz_init(v);for(long j=0;j<=fmpz_poly_degree(z);j++){fmpz_poly_get_coeff_fmpz(v,z,j);if(fmpz_is_zero(v))continue;fmpq_set_fmpz(c,v);ex[0]=j;fmpq_mpoly_push_term_fmpq_ui(p.p,c,ex,ctx);}fmpq_mpoly_sort_terms(p.p,ctx);fmpq_mpoly_combine_like_terms(p.p,ctx);assert(fmpq_mpoly_is_canonical(p.p,ctx));fmpz_clear(v);fmpq_clear(c);return p;}
void toZ(fmpz_poly_t z,const P&p,const fmpz_t D){fmpq_t c;fmpq_init(c);fmpz_t v;fmpz_init(v);ulong ex[3];for(long i=0;i<fmpq_mpoly_length(p.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(ex,p.p,i,ctx);assert(ex[1]==0&&ex[2]==0);fmpq_mpoly_get_term_coeff_fmpq(c,p.p,i,ctx);fmpz_mul(v,fmpq_numref(c),D);assert(fmpz_divisible(v,fmpq_denref(c)));fmpz_divexact(v,v,fmpq_denref(c));fmpz_poly_set_coeff_fmpz(z,ex[0],v);}fmpz_clear(v);fmpq_clear(c);}
P divideScalar(const P&p,const string&v){fmpq_t c;fmpq_init(c);assert(!fmpq_set_str(c,v.c_str(),10));fmpq_canonicalise(c);P r;fmpq_mpoly_scalar_div_fmpq(r.p,p.p,c,ctx);fmpq_clear(c);return r;}
int main(int argc,char**argv){assert(argc==3);start=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);{
 ifstream f(argv[1]);vector<vector<P>>M(16,vector<P>(16));for(auto&row:M)for(auto&p:row)readp(f,p);P H,K;readp(f,H);readp(f,K);vector<P>quot(16);for(auto&p:quot)readp(f,p);string DS,HL,KL;f>>DS>>HL>>KL;string extra;assert(!(f>>extra));fmpz_t D;fmpz_init(D);assert(!fmpz_set_str(D,DS.c_str(),10));assert(fmpz_sgn(D)>0);
 for(int j=0;j<16;j++){P rhs=mul(H,quot[j]);for(int i=0;i<16;i++)rhs=add(rhs,shift(M[i][j],1,i));assert(eq(rhs,shift(K,1,j)));}cerr<<"ALL16_SOURCE_DIVISIONS_PASS\n";
 fmpz_poly_mat_t A,B,V,prod;fmpz_poly_mat_init(A,16,16);fmpz_poly_mat_init(B,16,1);fmpz_poly_mat_init(V,16,1);fmpz_poly_mat_init(prod,16,1);for(int i=0;i<16;i++)for(int j=0;j<16;j++)toZ(fmpz_poly_mat_entry(A,i,j),M[i][j],D);fmpz_poly_one(fmpz_poly_mat_entry(B,0,0));fmpz_poly_t detz,delta;fmpz_poly_init(detz);fmpz_poly_init(delta);tick();fmpz_poly_mat_det(detz,A);assert(!fmpz_poly_is_zero(detz));cerr<<"DET16_PASS degree="<<fmpz_poly_degree(detz)<<"\n";tick();assert(fmpz_poly_mat_solve_fflu(V,delta,A,B));assert(!fmpz_poly_is_zero(delta));fmpz_poly_mat_mul(prod,A,V);for(int i=0;i<16;i++)assert(i==0?fmpz_poly_equal(fmpz_poly_mat_entry(prod,i,0),delta):fmpz_poly_is_zero(fmpz_poly_mat_entry(prod,i,0)));cerr<<"ONE_RHS_FFLU_SOLVE_PASS\n";
 P determinant=fromZ(detz),ann=fromZ(delta);vector<P>v(16);for(int j=0;j<16;j++)v[j]=fromZ(fmpz_poly_mat_entry(V,j,0));{ofstream checkpoint(string(argv[2])+".completed-det-solve.json");checkpoint<<"{\"actual_det16\":";jp(checkpoint,determinant);checkpoint<<",\"fflu_solver_delta\":";jp(checkpoint,ann);checkpoint<<",\"only16_solution_cofactors\":[";for(int j=0;j<16;j++){if(j)checkpoint<<',';jp(checkpoint,v[j]);}checkpoint<<"]}\n";}P Ve,Qsum;for(int j=0;j<16;j++){Ve=add(Ve,shift(v[j],1,j));Qsum=add(Qsum,mul(v[j],quot[j]));}P hc=neg(scalar(Qsum,DS)),kc=scalar(Ve,DS);assert(eq(add(mul(hc,H),mul(kc,K)),ann));P Horig=scalar(H,HL),Korig=scalar(K,KL),hco=divideScalar(hc,HL),kco=divideScalar(kc,KL);assert(eq(add(mul(hco,Horig),mul(kco,Korig)),ann));cerr<<"ORIGINAL_H1_K1_BEZOUT_PASS\n";
 fmpq_t content;fmpq_init(content);fmpq_mpoly_content(content,ann.p,ctx);P primitive;fmpq_mpoly_scalar_div_fmpq(primitive.p,ann.p,content,ctx);ulong ex[3];long amin=1000000000;for(long i=0;i<fmpq_mpoly_length(primitive.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(ex,primitive.p,i,ctx);assert(ex[1]==0&&ex[2]==0);amin=min(amin,(long)ex[0]);}P apow=shift(one(),0,amin),res=quotient(primitive,apow);assert(fmpq_mpoly_length(res.p,ctx)>0);long total=fmpq_mpoly_length(Ve.p,ctx)+fmpq_mpoly_length(Qsum.p,ctx)+fmpq_mpoly_length(hc.p,ctx)+fmpq_mpoly_length(kc.p,ctx);for(auto&p:v)total+=fmpq_mpoly_length(p.p,ctx);assert(total<150000);
 ofstream o(argv[2]);o<<"{\"integer_clearing_D\":\""<<DS<<"\",\"actual_det16\":";jp(o,determinant);o<<",\"fflu_solver_delta\":";jp(o,ann);o<<",\"solver_delta_equals_det\":"<<(eq(ann,determinant)?"true":"false")<<",\"solver_delta_equals_negative_det\":"<<(eq(ann,neg(determinant))?"true":"false")<<",\"only16_solution_cofactors\":[";for(int j=0;j<16;j++){if(j)o<<',';jp(o,v[j]);}o<<"],\"V_e_polynomial\":";jp(o,Ve);o<<",\"Qsum\":";jp(o,Qsum);o<<",\"monic_h1_cofactor\":";jp(o,hc);o<<",\"monic_k1_cofactor\":";jp(o,kc);o<<",\"original_h1_cofactor\":";jp(o,hco);o<<",\"original_k1_cofactor\":";jp(o,kco);char*cs=fmpq_get_str(nullptr,10,content);o<<",\"delta_QQ_content\":\""<<cs<<"\",\"delta_primitive\":";flint_free(cs);jp(o,primitive);o<<",\"removed_known_a_power\":"<<amin<<",\"remaining_non_a_annihilator\":";jp(o,res);o<<",\"total_solution_polynomial_terms\":"<<total<<",\"peak_single_polynomial_support\":"<<peak<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";fmpq_clear(content);fmpz_poly_clear(detz);fmpz_poly_clear(delta);fmpz_poly_mat_clear(A);fmpz_poly_mat_clear(B);fmpz_poly_mat_clear(V);fmpz_poly_mat_clear(prod);fmpz_clear(D);
}fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
