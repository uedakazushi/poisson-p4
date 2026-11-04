// Independent P/a^k representation (no min-a normalization), DEGREVLEX storage.
// Fixed source DAG, determinant and supplied division witnesses only; no search.
#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;fmpq_mpoly_ctx_t ctx;struct V;set<const V*>alive;long peak=0,reserve=0,identities=0;auto began=chrono::steady_clock::now();
struct V{fmpq_mpoly_t p;long den=0;V(){fmpq_mpoly_init(p,ctx);alive.insert(this);}V(const V&b):V(){fmpq_mpoly_set(p,b.p,ctx);den=b.den;}V(V&&b):V(){fmpq_mpoly_swap(p,b.p,ctx);den=b.den;b.den=0;}V&operator=(V b){fmpq_mpoly_swap(p,b.p,ctx);swap(den,b.den);return*this;}~V(){alive.erase(this);fmpq_mpoly_clear(p,ctx);}};
void guard(){long n=reserve;for(auto*x:alive)n+=fmpq_mpoly_length(x->p,ctx);peak=max(peak,n);if(n>150000)throw runtime_error("150k_live_terms_guard");if(chrono::duration<double>(chrono::steady_clock::now()-began).count()>49)throw runtime_error("49s_guard");}
V shift(const V&a,long k){assert(k>=0);V b;fmpq_t q;fmpq_init(q);ulong ex[3];for(long i=0;i<fmpq_mpoly_length(a.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(ex,a.p,i,ctx);ex[0]+=k;fmpq_mpoly_get_term_coeff_fmpq(q,a.p,i,ctx);fmpq_mpoly_set_coeff_fmpq_ui(b.p,q,ex,ctx);}fmpq_clear(q);guard();return b;}
V one(){V r;fmpq_mpoly_one(r.p,ctx);return r;}V addv(const V&a,const V&b,int sign=1){V r;r.den=max(a.den,b.den);V x=shift(a,r.den-a.den),y=shift(b,r.den-b.den);if(sign>0)fmpq_mpoly_add(r.p,x.p,y.p,ctx);else fmpq_mpoly_sub(r.p,x.p,y.p,ctx);guard();return r;}
V times(const V&a,const V&b){V r;r.den=a.den+b.den;fmpq_mpoly_mul(r.p,a.p,b.p,ctx);guard();return r;}
void eq(const V&a,const V&b){identities++;V x=shift(a,b.den),y=shift(b,a.den);assert(fmpq_mpoly_equal(x.p,y.p,ctx));}
bool zero(const V&a){return fmpq_mpoly_is_zero(a.p,ctx);}
V readpoly(istream&f){long n;f>>n;assert(n>=0&&n<150000);vector<array<long,3>>ex(n);vector<string>coeff(n);long den=0;for(long i=0;i<n;i++){f>>ex[i][0]>>ex[i][1]>>ex[i][2]>>coeff[i];assert(ex[i][1]>=0&&ex[i][2]>=0);den=max(den,-ex[i][0]);}V r;r.den=den;fmpq_t q;fmpq_init(q);for(long i=0;i<n;i++){ulong e[3]={ulong(ex[i][0]+den),ulong(ex[i][1]),ulong(ex[i][2])};assert(!fmpq_set_str(q,coeff[i].c_str(),10));fmpq_canonicalise(q);fmpq_mpoly_set_coeff_fmpq_ui(r.p,q,e,ctx);}fmpq_clear(q);guard();return r;}
V divQQ(const V&a,const string&s){V r=a;fmpq_t q;fmpq_init(q);assert(!fmpq_set_str(q,s.c_str(),10));fmpq_canonicalise(q);assert(!fmpq_is_zero(q));fmpq_mpoly_scalar_div_fmpq(r.p,r.p,q,ctx);fmpq_clear(q);guard();return r;}
bool lex_reduced(const V&r,const V&f){assert(r.den==0&&f.den==0);array<ulong,3>lead{};ulong e[3];for(long i=0;i<fmpq_mpoly_length(f.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(e,f.p,i,ctx);array<ulong,3>x={e[0],e[1],e[2]};lead=max(lead,x);}for(long i=0;i<fmpq_mpoly_length(r.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(e,r.p,i,ctx);if(e[0]>=lead[0]&&e[1]>=lead[1]&&e[2]>=lead[2])return false;}return true;}
int main(int argc,char**argv){assert(argc==3);began=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,3,ORD_DEGREVLEX);try{{
 ifstream in(argv[1]);long n;in>>reserve>>n;assert(n==195);vector<V>v;v.reserve(n);
 for(long i=0;i<n;i++){char op;in>>op;V q;if(op=='z'){}else if(op=='o')q=one();else if(op=='p'){long aden;string c;in>>aden>>c;assert(aden>=0);q=divQQ(readpoly(in),c);q.den+=aden;}else{long a,b;in>>a>>b;assert(a<i&&b<i);q=op=='m'?times(v[a],v[b]):addv(v[a],v[b],op=='s'?-1:1);}v.push_back(move(q));guard();}
 long nc;in>>nc;assert(nc==108);for(long i=0;i<nc;i++){long j;in>>j;eq(v.at(j),readpoly(in));}
 vector<V>Q;for(long i=0;i<9;i++){long j;in>>j;V p=readpoly(in);eq(v.at(j),p);Q.push_back(move(p));}
 array<array<long,4>,3>roots;for(auto&row:roots)for(long&j:row)in>>j;int qidx,js[3];in>>qidx>>js[0]>>js[1]>>js[2];assert(qidx==7);vector<vector<V>>M(3);
 for(long i=0;i<3;i++)for(long j=0;j<3;j++){V p=addv(times(Q[qidx],v[roots[i][j]]),times(v[roots[i][3]],Q[js[j]]),-1);eq(p,readpoly(in));M[i].push_back(move(p));}
 vector<V>cof;cof.push_back(addv(times(M[1][1],M[2][2]),times(M[1][2],M[2][1]),-1));cof.push_back(addv(times(M[1][0],M[2][2]),times(M[1][2],M[2][0]),-1));cof.push_back(addv(times(M[1][0],M[2][1]),times(M[1][1],M[2][0]),-1));for(auto&c:cof)eq(c,readpoly(in));
 V laplace=addv(addv(times(M[0][0],cof[0]),times(M[0][1],cof[1]),-1),times(M[0][2],cof[2])),six;array<int,3>perm={0,1,2};do{V term=one();long parity=0;for(long i=0;i<3;i++){term=times(term,M[i][perm[i]]);for(long j=i+1;j<3;j++)parity+=perm[i]>perm[j];}six=addv(six,term,parity%2?-1:1);}while(next_permutation(perm.begin(),perm.end()));eq(laplace,six);V original=readpoly(in);eq(laplace,original);assert(!zero(original));
 V h=readpoly(in);eq(h,Q[qidx]);V h2=times(h,h),pure4=readpoly(in);eq(original,times(h2,pure4));assert(!zero(pure4));vector<V>factors;for(long i=0;i<6;i++)factors.push_back(readpoly(in));long steps;in>>steps;assert(steps==9);V current=pure4,product=one();vector<long>powers(6);long refusals=0;
 for(long i=0;i<steps;i++){long idx,power,exact;in>>idx>>power>>exact;V before=readpoly(in),quot=readpoly(in),rem=readpoly(in);eq(before,current);
  if(idx==0){assert(i==0&&power==4&&exact==1&&zero(rem));V ap=one();for(long k=0;k<power;k++)ap=times(ap,factors[0]);eq(before,times(ap,quot));product=times(product,ap);powers[0]=power;current=quot;}
  else{assert(1<=idx&&idx<6&&power==powers[idx]);eq(before,addv(times(factors[idx],quot),rem));assert(bool(exact)==zero(rem));assert(lex_reduced(rem,factors[idx]));if(exact){current=quot;product=times(product,factors[idx]);powers[idx]++;}else refusals++;}
 }
 V remaining=readpoly(in);eq(current,remaining);eq(pure4,times(product,remaining));eq(original,times(h2,times(product,remaining)));assert(powers==vector<long>({4,0,0,0,3,0})&&refusals==5);string extra;assert(!(in>>extra));guard();
 ofstream out(argv[2]);out<<"{\"all_passed\":true,\"independent_fraction_representation\":\"Numerator/a^k_without_min_a_normalization\",\"storage_order\":\"DEGREVLEX\",\"source_scalar_nodes\":195,\"original108_values_checked\":true,\"all9_source_Q9_tails_checked\":true,\"all9_actual_minor_entries_checked\":true,\"fixed3_cofactor_and_6permutation_checked\":true,\"h22_squared_times_actualpure4_equals_fixed3\":true,\"all9_supplied_division_steps_checked\":true,\"all5_refusal_remainders_lex_reduced\":true,\"named_powers\":[4,0,0,0,3,0],\"all_unit_product_reconstructions_checked\":true,\"literal_polynomial_identities\":"<<identities<<",\"peak_live_terms_including_inputreserve\":"<<peak<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";
 }}catch(const exception&e){cerr<<e.what()<<'\n';fmpq_mpoly_ctx_clear(ctx);return 2;}fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
