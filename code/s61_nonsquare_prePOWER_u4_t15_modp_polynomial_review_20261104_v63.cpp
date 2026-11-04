// Independent symbolic F32003[a,e,d,a^-1] replay. No parameter evaluation.
#include <flint/fmpq.h>
#include <flint/nmod_mpoly.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
constexpr ulong PRIME=32003;nmod_mpoly_ctx_t ctx;struct P;set<const P*>live;
long reserve_terms=0,peak=0,checks=0,denominator_checks=0;auto begun=chrono::steady_clock::now();ofstream progress;
struct P{nmod_mpoly_t p;long den=0;P(){nmod_mpoly_init(p,ctx);live.insert(this);}P(const P&a):P(){nmod_mpoly_set(p,a.p,ctx);den=a.den;}P(P&&a):P(){nmod_mpoly_swap(p,a.p,ctx);den=a.den;a.den=0;}P&operator=(P a){nmod_mpoly_swap(p,a.p,ctx);swap(den,a.den);return *this;}~P(){live.erase(this);nmod_mpoly_clear(p,ctx);}};
void guard(){long n=reserve_terms+15000;for(auto*p:live)n+=nmod_mpoly_length(p->p,ctx);peak=max(peak,n);if(n>150000)throw runtime_error("150000 live coefficient guard");if(chrono::duration<double>(chrono::steady_clock::now()-begun).count()>48)throw runtime_error("48 second guard");}
ulong pw(ulong a,ulong e){ulong x=1;for(;e;e>>=1,a=a*a%PRIME)if(e&1)x=x*a%PRIME;return x;}
ulong rational_coefficient(const string&s){fmpq_t q;fmpq_init(q);assert(!fmpq_set_str(q,s.c_str(),10));fmpq_canonicalise(q);ulong n=fmpz_fdiv_ui(fmpq_numref(q),PRIME),d=fmpz_fdiv_ui(fmpq_denref(q),PRIME);assert(d);denominator_checks++;fmpq_clear(q);return n*pw(d,PRIME-2)%PRIME;}
P one(){P p;nmod_mpoly_one(p.p,ctx);return p;}
P shift(const P&a,long k){assert(k>=0);P p;ulong ex[3];for(long i=0;i<nmod_mpoly_length(a.p,ctx);i++){nmod_mpoly_get_term_exp_ui(ex,a.p,i,ctx);ex[0]+=k;nmod_mpoly_set_coeff_ui_ui(p.p,nmod_mpoly_get_term_coeff_ui(a.p,i,ctx),ex,ctx);}guard();return p;}
P sum_poly(const P&a,const P&b,int s=1){P p;p.den=max(a.den,b.den);P x=shift(a,p.den-a.den),y=shift(b,p.den-b.den);if(s==1)nmod_mpoly_add(p.p,x.p,y.p,ctx);else nmod_mpoly_sub(p.p,x.p,y.p,ctx);guard();return p;}
P times(const P&a,const P&b){P p;p.den=a.den+b.den;nmod_mpoly_mul(p.p,a.p,b.p,ctx);guard();return p;}
P divide_scalar(const P&a,const string&s){P p=a;ulong q=rational_coefficient(s);assert(q);nmod_mpoly_scalar_mul_ui(p.p,p.p,pw(q,PRIME-2),ctx);guard();return p;}
P polynomial(istream&in){long n;in>>n;assert(n>=0&&n<150000);vector<array<long,3>>ex(n);vector<string>co(n);long den=0;for(long i=0;i<n;i++){in>>ex[i][0]>>ex[i][1]>>ex[i][2]>>co[i];assert(ex[i][1]>=0&&ex[i][2]>=0);den=max(den,-ex[i][0]);}P p;p.den=den;for(long i=0;i<n;i++){ulong e[3]={ulong(ex[i][0]+den),ulong(ex[i][1]),ulong(ex[i][2])};nmod_mpoly_set_coeff_ui_ui(p.p,rational_coefficient(co[i]),e,ctx);}guard();return p;}
void equal(const P&a,const P&b){P x=shift(a,b.den),y=shift(b,a.den);assert(nmod_mpoly_equal(x.p,y.p,ctx));checks++;}
P power(const P&a,long n){assert(n>=0);P p=one();for(long i=0;i<n;i++)p=times(p,a);return p;}
struct Node{char op;long a=-1,b=-1;unique_ptr<P>atom;};
int main(int argc,char**argv){assert(argc==4);nmod_mpoly_ctx_init(ctx,3,ORD_DEGREVLEX,PRIME);progress.open(argv[3]);try{{
 ifstream in(argv[1]);long n;in>>reserve_terms>>n;assert(n==768);vector<Node>nodes(n);vector<long>uses(n);vector<unique_ptr<P>>vals(n);
 for(long i=0;i<n;i++){in>>nodes[i].op;if(nodes[i].op=='p')nodes[i].atom=make_unique<P>(polynomial(in));else{in>>nodes[i].a>>nodes[i].b;assert(nodes[i].a<i&&nodes[i].b<i);uses[nodes[i].a]++;uses[nodes[i].b]++;}}
 long count;in>>count;for(long i=0,j;i<count;i++){in>>j;uses[j]++;}
 progress<<"all768_input_atoms_QQ_denominators_are_32003_units\n"<<flush;
 for(long i=0;i<n;i++){auto&x=nodes[i];if(x.op=='p')vals[i]=move(x.atom);else{P p=x.op=='m'?times(*vals[x.a],*vals[x.b]):sum_poly(*vals[x.a],*vals[x.b],x.op=='s'?-1:1);vals[i]=make_unique<P>(move(p));for(long j:{x.a,x.b})if(--uses[j]==0)vals[j].reset();}guard();if(!uses[i])vals[i].reset();}
 progress<<"all768_symbolic_modular_operations_completed\n"<<flush;
 long specials;in>>specials;assert(specials==3);for(long k=0;k<3;k++){long j;in>>j;equal(*vals[j],polynomial(in));}
 long output;in>>output;assert(output==767);P expected=polynomial(in),cleared=polynomial(in);equal(*vals[output],expected);equal(expected,cleared);assert(!nmod_mpoly_is_zero(expected.p,ctx));
 progress<<"all_public_QQ7022_coefficients_reduced_and_equal;nonzero\n"<<flush;
 array<array<long,4>,4>A;array<long,4>target;array<array<long,5>,3>source;array<long,3>srcids;array<long,5>bids;
 for(auto&r:A)for(auto&j:r)in>>j;for(auto&j:target)in>>j;for(auto&r:source)for(auto&j:r)in>>j;for(auto&j:srcids)in>>j;for(auto&j:bids)in>>j;
 P delta2=polynomial(in);assert(!nmod_mpoly_is_zero(delta2.p,ctx));array<array<P,2>,3>N;for(auto&r:N)for(auto&x:r)x=polynomial(in);vector<P>units;for(long i=0;i<6;i++){units.push_back(polynomial(in));assert(!nmod_mpoly_is_zero(units.back().p,ctx));}
 array<array<long,6>,3>powers;array<string,3>constants;array<long,6>lcm;for(auto&r:powers)for(auto&j:r)in>>j;for(auto&c:constants)in>>c;for(auto&j:lcm)in>>j;
 vector<P>co;co.push_back(delta2);co.push_back(sum_poly(times(N[0][0],N[2][1]),times(N[0][1],N[2][0]),-1));co.push_back(sum_poly(times(N[0][1],N[1][0]),times(N[0][0],N[1][1]),-1));
 for(long i=0;i<3;i++){P p=divide_scalar(co[i],constants[i]);if(i)p=sum_poly(P(),p,-1);for(long j=0;j<6;j++)p=times(p,power(units[j],lcm[j]-powers[i][j]));equal(p,*vals[srcids[i]]);}
 vector<P>bh;for(long j=0;j<5;j++){P p;for(long i=0;i<3;i++)p=sum_poly(p,times(*vals[srcids[i]],*vals[source[i][j]]));equal(p,*vals[bids[j]]);bh.push_back(move(p));}
 progress<<"actual3_source_cofactors_and_all5_Bhat_symbolic_identities_verified\n"<<flush;
 vector<vector<const P*>>border(5,vector<const P*>(5));for(long i=0;i<4;i++){for(long j=0;j<4;j++)border[i][j]=vals[A[i][j]].get();border[i][4]=vals[target[i]].get();}for(long j=0;j<5;j++)border[4][j]=&bh[j];
 array<int,5>perm={0,1,2,3,4};P det5;do{long s=1;for(long i=0;i<5;i++)for(long j=i+1;j<5;j++)if(perm[i]>perm[j])s=-s;P term=one();for(long i=0;i<5;i++)term=times(term,*border[i][perm[i]]);det5=sum_poly(det5,term,s);}while(next_permutation(perm.begin(),perm.end()));equal(det5,expected);
 progress<<"independent_actual_bordered5_determinant_equal_and_nonzero\n"<<flush;
 // New appended source chart polynomials: all must remain nonzero. The
 // fixed direct P4 determinant is also independently rebuilt over Fp.
 P savedD4=polynomial(in),remainingR=polynomial(in),pairp=polynomial(in),delta3=polynomial(in),Delta5=polynomial(in);
 for(auto*p:{&savedD4,&remainingR,&pairp,&delta3,&Delta5})assert(!nmod_mpoly_is_zero(p->p,ctx));
 array<int,4>p4={0,1,2,3};P det4;do{long s=1;for(long i=0;i<4;i++)for(long j=i+1;j<4;j++)if(p4[i]>p4[j])s=-s;P term=one();for(long i=0;i<4;i++)term=times(term,*vals[A[i][p4[i]]]);det4=sum_poly(det4,term,s);}while(next_permutation(p4.begin(),p4.end()));equal(det4,savedD4);
 equal(savedD4,times(times(power(units[0],4),power(units[4],3)),remainingR));equal(Delta5,times(delta2,delta3));
 equal(delta2,divide_scalar(sum_poly(P(),times(units[3],power(pairp,2)),-1),"5/6"));
 string extra;assert(!(in>>extra));guard();ulong ex[3];nmod_mpoly_get_term_exp_ui(ex,expected.p,0,ctx);ulong c=nmod_mpoly_get_term_coeff_ui(expected.p,0,ctx);
 ofstream out(argv[2]);out<<"{\"all_passed\":true,\"prime\":32003,\"method\":\"symbolic_nmod_mpoly_Laurent_P_over_a_power\",\"all768_source_operations_verified_mod_p\":true,\"all_public7022_QQ_coefficients_reduced_and_equal\":true,\"actual_fixed_bordered5_determinant_verified_nonzero_mod_p\":true,\"all3_source_cofactors_and5_Bhat_verified_mod_p\":true,\"P4_D4_R434_p_Delta5_and_all_named_denominator_polynomials_nonzero_mod_p\":true,\"QQ7022_exact_independent_verification_claimed\":false,\"modular_result_terms\":"<<nmod_mpoly_length(expected.p,ctx)<<",\"nonzero_term\":{\"Laurent_exponents\":["<<long(ex[0])-expected.den<<","<<ex[1]<<","<<ex[2]<<"],\"coefficient_mod32003\":"<<c<<"},\"QQ_denominator_unit_checks\":"<<denominator_checks<<",\"fixed_polynomial_identities\":"<<checks<<",\"peak_live_coefficients_with_reserves\":"<<peak<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-begun).count()<<"}\n";
 }}catch(const exception&e){progress<<"STOP "<<e.what()<<'\n'<<flush;cerr<<e.what()<<'\n';nmod_mpoly_ctx_clear(ctx);return 2;}nmod_mpoly_ctx_clear(ctx);flint_cleanup();}
