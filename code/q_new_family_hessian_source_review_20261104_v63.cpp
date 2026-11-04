// Independent direct H16-column recurrence; no author source/helper include.
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
using namespace std;
fmpz_mpoly_ctx_t ctx;
long maxterms=0,products=0;
struct P{
 fmpz_mpoly_t v;P(){fmpz_mpoly_init(v,ctx);}P(const P&a):P(){fmpz_mpoly_set(v,a.v,ctx);}
 P(P&&a):P(){fmpz_mpoly_swap(v,a.v,ctx);}P&operator=(P a){fmpz_mpoly_swap(v,a.v,ctx);return *this;}
 ~P(){fmpz_mpoly_clear(v,ctx);}
};
void check(const P&a){long n=fmpz_mpoly_length(a.v,ctx);maxterms=max(maxterms,n);if(n>150000)throw runtime_error("term cap");}
P add(const P&a,const P&b){P c;fmpz_mpoly_add(c.v,a.v,b.v,ctx);check(c);return c;}
P scale(const P&a,long c){P b;fmpz_mpoly_scalar_mul_si(b.v,a.v,c,ctx);return b;}
P mul(const P&a,const P&b){P c;fmpz_mpoly_mul(c.v,a.v,b.v,ctx);products++;check(c);return c;}
void readpoly(istream&f,P&a){long n;f>>n;assert(n>=0);fmpz_t z;fmpz_init(z);ulong e[4];for(long k=0;k<n;k++){string c;f>>c;for(auto&x:e)f>>x;assert(f&&fmpz_set_str(z,c.c_str(),10)==0);fmpz_mpoly_push_term_fmpz_ui(a.v,z,e,ctx);}fmpz_mpoly_sort_terms(a.v,ctx);fmpz_mpoly_combine_like_terms(a.v,ctx);fmpz_clear(z);check(a);}
using M=vector<P>;
M mm(const M&a,const M&b,int n){M c(n*n);for(int i=0;i<n;i++)for(int j=0;j<n;j++)for(int k=0;k<n;k++)if(!fmpz_mpoly_is_zero(a[n*i+k].v,ctx)&&!fmpz_mpoly_is_zero(b[n*k+j].v,ctx))c[n*i+j]=add(c[n*i+j],mul(a[n*i+k],b[n*k+j]));return c;}
P pairtr(const M&a,const M&b,int n){P z;for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(!fmpz_mpoly_is_zero(a[n*i+j].v,ctx)&&!fmpz_mpoly_is_zero(b[n*j+i].v,ctx))z=add(z,mul(a[n*i+j],b[n*j+i]));return z;}
struct Target{string name;long degree,factor,bpow;fmpz_t den;P numerator;Target(){fmpz_init(den);}~Target(){fmpz_clear(den);}};
void verify(const P&actual,Target&t){
 // Both sides retain H^degree. Compare exact integer polynomials after
 // multiplying only a power of T, with B=T².
 P lhs(actual),rhs(t.numerator),shift;fmpz_t z;fmpz_init(z);
 fmpz_mpoly_scalar_mul_fmpz(lhs.v,lhs.v,t.den,ctx);
 fmpz_mpoly_gen(shift.v,3,ctx);P power;
 assert(fmpz_mpoly_pow_ui(power.v,shift.v,2*t.bpow,ctx));lhs=mul(lhs,power);
 assert(fmpz_mpoly_pow_ui(power.v,shift.v,t.degree,ctx));rhs=mul(rhs,power);
 fmpz_set_ui(z,1080);fmpz_pow_ui(z,z,t.degree);fmpz_mul_ui(z,z,t.factor);fmpz_mpoly_scalar_mul_fmpz(rhs.v,rhs.v,z,ctx);
 assert(fmpz_mpoly_equal(lhs.v,rhs.v,ctx));
 ulong e[4];for(long k=0;k<fmpz_mpoly_length(actual.v,ctx);k++){fmpz_mpoly_get_term_exp_ui(e,actual.v,k,ctx);assert(e[1]%2==0&&e[3]%2==0);}
 cerr<<"VERIFIED "<<t.name<<" RAW_TERMS "<<fmpz_mpoly_length(actual.v,ctx)<<" MAX_TERMS "<<maxterms<<endl;fmpz_clear(z);
}
int main(int argc,char**argv){
 assert(argc==3);fmpz_mpoly_ctx_init(ctx,4,ORD_LEX);{
 ifstream f(argv[1]);int n,width;f>>n>>width;assert(n==16&&width==136);vector<P>beta(16);for(auto&p:beta)readpoly(f,p);
 vector<pair<int,int>>mon(136);for(auto&v:mon)f>>v.first>>v.second;
 vector<P>quad;for(auto ij:mon)quad.push_back(mul(beta[ij.first],beta[ij.second]));
 M C(16),D(16),Theta(256),H(256);vector<M*>matrices={&C,&D,&Theta,&H};
 for(auto*M:matrices)for(auto&p:*M){for(int k=0;k<136;k++){long c;f>>c;if(c)p=add(p,scale(quad[k],c));}}
 Target target[4];for(auto&t:target){string den;f>>t.name>>t.degree>>t.factor>>t.bpow>>den;assert(fmpz_set_str(t.den,den.c_str(),10)==0);readpoly(f,t.numerator);}string extra;assert(f&&!(f>>extra));
 M C2=mm(C,C,4),C3=mm(C2,C,4),D2=mm(D,D,4),D3=mm(D2,D,4),CD=mm(C,D,4),C3D=mm(C3,D,4);
 P word0=pairtr(C3D,CD,4),word1=pairtr(C3,D3,4);verify(word0,target[0]);verify(word1,target[1]);
 P trace5,trace7;
 // Directly transport the16 genuine Hessian columns, one live column at a
 // time. This uses neither rank10 factorization nor a pure-character term.
 for(int column=0;column<16;column++){
  vector<P>level(16);for(int i=0;i<16;i++)level[i]=H[16*i+column];
  for(int degree=1;degree<=7;degree++){
   vector<P>next(16);for(int i=0;i<16;i++)for(int j=0;j<16;j++)if(!fmpz_mpoly_is_zero(Theta[16*i+j].v,ctx)&&!fmpz_mpoly_is_zero(level[j].v,ctx))next[i]=add(next[i],mul(Theta[16*i+j],level[j]));level=move(next);
   if(degree==5)trace5=add(trace5,level[column]);if(degree==7)trace7=add(trace7,level[column]);
  }
  cerr<<"DIRECT_HESSIAN_COLUMN "<<column<<" COMPLETE MAX_TERMS "<<maxterms<<endl;
 }
 verify(trace5,target[2]);verify(trace7,target[3]);
 ofstream out(argv[2]);out<<"{\"all_passed\":true,\"all4_exact_coefficient_maps\":true,\"direct_Hessian_columns\":16,\"max_terms\":"<<maxterms<<",\"products\":"<<products<<",\"raw_term_counts\":["<<fmpz_mpoly_length(word0.v,ctx)<<','<<fmpz_mpoly_length(word1.v,ctx)<<','<<fmpz_mpoly_length(trace5.v,ctx)<<','<<fmpz_mpoly_length(trace7.v,ctx)<<"]}\n";
 }fmpz_mpoly_ctx_clear(ctx);flint_cleanup();}
