// Fixed35 exact3x3 determinants, source-scaling reconstruction, named factors ONLY.
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
using namespace std;fmpq_mpoly_ctx_t ctx;
struct P{fmpq_mpoly_t p;P(){fmpq_mpoly_init(p,ctx);}P(const P&v):P(){fmpq_mpoly_set(p,v.p,ctx);}P(P&&v):P(){fmpq_mpoly_swap(p,v.p,ctx);}P&operator=(P v){fmpq_mpoly_swap(p,v.p,ctx);return*this;}~P(){fmpq_mpoly_clear(p,ctx);}};
long nt(const P&p){return fmpq_mpoly_length(p.p,ctx);}
void readp(istream&f,P&p){long n;f>>n;assert(f&&n>=0);fmpq_t c;fmpq_init(c);ulong ex[3];for(long i=0;i<n;i++){string s;f>>s>>ex[0]>>ex[1]>>ex[2];assert(f&&!fmpq_set_str(c,s.c_str(),10));fmpq_canonicalise(c);fmpq_mpoly_push_term_fmpq_ui(p.p,c,ex,ctx);}fmpq_mpoly_sort_terms(p.p,ctx);fmpq_mpoly_combine_like_terms(p.p,ctx);assert(nt(p)==n);fmpq_clear(c);}
P mul(const P&a,const P&b){P z;fmpq_mpoly_mul(z.p,a.p,b.p,ctx);assert(nt(z)<=150000);return z;}
P add(const P&a,const P&b,bool minus=false){P z;if(minus)fmpq_mpoly_sub(z.p,a.p,b.p,ctx);else fmpq_mpoly_add(z.p,a.p,b.p,ctx);assert(nt(z)<=150000);return z;}
P det_permutation(const vector<vector<P>>&A,const array<int,3>&I,const array<int,3>&J){P z;array<int,3>p={0,1,2};do{P term=mul(mul(A[I[0]][J[p[0]]],A[I[1]][J[p[1]]]),A[I[2]][J[p[2]]]);int inv=0;for(int i=0;i<3;i++)for(int j=i+1;j<3;j++)inv+=p[i]>p[j];z=add(z,term,inv%2);}while(next_permutation(p.begin(),p.end()));return z;}
P det_laplace(const vector<vector<P>>&A,const array<int,3>&I){P z;for(int j=0;j<3;j++){int k=(j+1)%3,l=(j+2)%3;P co=add(mul(A[I[1]][k],A[I[2]][l]),mul(A[I[1]][l],A[I[2]][k]),true);z=add(z,mul(A[I[0]][j],co));}return z;}
void jsonp(ostream&f,const P&p){f<<'[';fmpq_t c;fmpq_init(c);ulong e[3];for(long i=0;i<nt(p);i++){if(i)f<<',';fmpq_mpoly_get_term_coeff_fmpq(c,p.p,i,ctx);fmpq_mpoly_get_term_exp_ui(e,p.p,i,ctx);char*s=fmpq_get_str(nullptr,10,c);f<<"[["<<e[0]<<','<<e[1]<<','<<e[2]<<"],\""<<s<<"\"]";flint_free(s);}f<<']';fmpq_clear(c);assert(f.tellp()<30*1024*1024);}
string qstr(const fmpq_t q){char*s=fmpq_get_str(nullptr,10,q);string z=s;flint_free(s);return z;}
int main(int argc,char**argv){
 assert(argc==3);auto start=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);
 {
  ifstream f(argv[1]);ofstream out(argv[2]);int nr,nc,nk;f>>nr>>nc>>nk;assert(nr==7&&nc==3&&nk==6);vector<P>known(nk),ru(nr),cu(nc),D(nc);for(P&p:known)readp(f,p);
  vector<vector<P>>original(nr,vector<P>(nc)),A(nr,vector<P>(nc));for(auto&row:original)for(P&p:row)readp(f,p);for(auto&row:A)for(P&p:row)readp(f,p);for(P&p:ru)readp(f,p);for(P&p:cu)readp(f,p);for(P&p:D)readp(f,p);
  array<int,3>si,sj;for(int&i:si)f>>i;for(int&j:sj)f>>j;P oldDelta;readp(f,oldDelta);string extra;assert(!(f>>extra));
  for(int i=0;i<nr;i++)for(int j=0;j<nc;j++){P left=mul(mul(ru[i],A[i][j]),cu[j]),right=mul(D[j],original[i][j]);assert(fmpq_mpoly_equal(left.p,right.p,ctx));}
  P selected=det_permutation(original,si,sj);assert(fmpq_mpoly_equal(selected.p,oldDelta.p,ctx));P dp;fmpq_mpoly_one(dp.p,ctx);for(const P&p:D)dp=mul(dp,p);
  long count=0,zeros=0,units=0,totalterms=0,maxterms=0,primterms=0,maxprim=0;out<<"{\"all35_source_scalings_verified\":true,\"all21_entry_scalings_verified\":true,\"selected407Delta_exact_verified\":true,\"records\":[";
  for(int i=0;i<5;i++)for(int j=i+1;j<6;j++)for(int k=j+1;k<7;k++){
   array<int,3>I={i,j,k},J={0,1,2};P det=det_permutation(A,I,J),second=det_laplace(A,I);assert(fmpq_mpoly_equal(det.p,second.p,ctx));P source=det_laplace(original,I),scale;fmpq_mpoly_one(scale.p,ctx);for(int u:I)scale=mul(scale,ru[u]);for(const P&p:cu)scale=mul(scale,p);P left=mul(dp,source),right=mul(scale,det);assert(fmpq_mpoly_equal(left.p,right.p,ctx));
   if(count)out<<',';out<<"{\"rows7\":["<<i<<','<<j<<','<<k<<"],\"prepared_integer_minor\":";jsonp(out,det);out<<",\"actual_source_minor\":";jsonp(out,source);
   if(!nt(det)){zeros++;out<<",\"zero\":true}";count++;continue;}
   totalterms+=nt(det);maxterms=max(maxterms,nt(det));P q=det,next;fmpq_t c;fmpq_init(c);fmpq_mpoly_get_term_coeff_fmpq(c,det.p,0,ctx);string constant=qstr(c);fmpq_mpoly_scalar_div_fmpq(q.p,q.p,c,ctx);array<long,6>powers{};
   for(int u=0;u<6;u++)while(fmpq_mpoly_divides(next.p,q.p,known[u].p,ctx)){q=next;powers[u]++;}
   P primitive=q;fmpz_t den;fmpz_init(den);fmpq_mpoly_get_denominator(den,q.p,ctx);fmpq_t clear;fmpq_init(clear);fmpq_set_fmpz(clear,den);fmpq_mpoly_scalar_mul_fmpq(primitive.p,q.p,clear,ctx);
   P expected=primitive;for(int u=0;u<6;u++)for(int n=0;n<powers[u];n++)expected=mul(expected,known[u]);fmpq_div(c,c,clear);fmpq_mpoly_scalar_mul_fmpq(expected.p,expected.p,c,ctx);assert(fmpq_mpoly_equal(expected.p,det.p,ctx));
   bool unit=fmpq_mpoly_is_fmpq(q.p,ctx);units+=unit;primterms+=nt(primitive);maxprim=max(maxprim,nt(primitive));
   out<<",\"zero\":false,\"named_unit_constant\":\""<<qstr(c)<<"\",\"named_a_d_F_B_C_J_powers\":[";for(int u=0;u<6;u++){if(u)out<<',';out<<powers[u];}out<<"],\"primitive_integer_polynomial\":";jsonp(out,primitive);out<<",\"single_knownunit_minor\":"<<(unit?"true":"false")<<'}';
   fmpq_clear(c);fmpq_clear(clear);fmpz_clear(den);count++;assert(chrono::duration<double>(chrono::steady_clock::now()-start).count()<50);
  }
  P product;fmpq_mpoly_one(product.p,ctx);for(const P&p:known)product=mul(product,p);out<<"],\"FULLU_product\":";jsonp(out,product);
  out<<",\"all_passed\":true,\"minor_count\":"<<count<<",\"zero_minors\":"<<zeros<<",\"single_namedunit_minors\":"<<units<<",\"nonzero_prepared_minor_totalterms\":"<<totalterms<<",\"max_prepared_minor_terms\":"<<maxterms<<",\"primitive_totalterms_before_QQ_dedup\":"<<primterms<<",\"max_primitive_terms\":"<<maxprim<<",\"FULLU_product_terms\":"<<nt(product)<<",\"FULLU_product_degree\":"<<fmpq_mpoly_total_degree_si(product.p,ctx)<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";assert(count==35);
 }
 fmpq_mpoly_ctx_clear(ctx);flint_cleanup();return 0;
}
