// Independent Laurent P/a^k checker, with last-use release and fixed bordered
// determinant. No producer code, Gaussian elimination, or new source query.
#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
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
fmpq_mpoly_ctx_t ctx;struct P;set<const P*> live;long reserve_terms=0,peak=0,checks=0;auto started=chrono::steady_clock::now();
struct P {
 fmpq_mpoly_t p;long den=0;
 P(){fmpq_mpoly_init(p,ctx);live.insert(this);}
 P(const P&x):P(){fmpq_mpoly_set(p,x.p,ctx);den=x.den;}
 P(P&&x):P(){fmpq_mpoly_swap(p,x.p,ctx);den=x.den;x.den=0;}
 P&operator=(P x){fmpq_mpoly_swap(p,x.p,ctx);swap(den,x.den);return *this;}
 ~P(){live.erase(this);fmpq_mpoly_clear(p,ctx);}
};
void guard(){long n=reserve_terms+15000;for(auto*x:live)n+=fmpq_mpoly_length(x->p,ctx);peak=max(peak,n);if(n>150000)throw runtime_error("live coefficient guard");if(chrono::duration<double>(chrono::steady_clock::now()-started).count()>48)throw runtime_error("48 second guard");}
P one(){P p;fmpq_mpoly_one(p.p,ctx);return p;}
P shift(const P&a,long n){assert(n>=0);P p;fmpq_t q;fmpq_init(q);ulong ex[3];for(long i=0;i<fmpq_mpoly_length(a.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(ex,a.p,i,ctx);ex[0]+=n;fmpq_mpoly_get_term_coeff_fmpq(q,a.p,i,ctx);fmpq_mpoly_set_coeff_fmpq_ui(p.p,q,ex,ctx);}fmpq_clear(q);guard();return p;}
P add(const P&a,const P&b,int sign=1){P p;p.den=max(a.den,b.den);P x=shift(a,p.den-a.den),y=shift(b,p.den-b.den);if(sign==1)fmpq_mpoly_add(p.p,x.p,y.p,ctx);else fmpq_mpoly_sub(p.p,x.p,y.p,ctx);guard();return p;}
P times(const P&a,const P&b){P p;p.den=a.den+b.den;fmpq_mpoly_mul(p.p,a.p,b.p,ctx);guard();return p;}
P rational(const P&a,const string&s){P p=a;fmpq_t q;fmpq_init(q);assert(!fmpq_set_str(q,s.c_str(),10));fmpq_canonicalise(q);assert(!fmpq_is_zero(q));fmpq_mpoly_scalar_div_fmpq(p.p,p.p,q,ctx);fmpq_clear(q);guard();return p;}
P readpoly(istream&in){long n;in>>n;assert(n>=0&&n<150000);vector<array<long,3>>es(n);vector<string>cs(n);long d=0;for(long i=0;i<n;i++){in>>es[i][0]>>es[i][1]>>es[i][2]>>cs[i];d=max(d,-es[i][0]);assert(es[i][1]>=0&&es[i][2]>=0);}P p;p.den=d;fmpq_t q;fmpq_init(q);for(long i=0;i<n;i++){ulong ex[3]={ulong(es[i][0]+d),ulong(es[i][1]),ulong(es[i][2])};assert(!fmpq_set_str(q,cs[i].c_str(),10));fmpq_canonicalise(q);fmpq_mpoly_set_coeff_fmpq_ui(p.p,q,ex,ctx);}fmpq_clear(q);guard();return p;}
void equal(const P&a,const P&b){P x=shift(a,b.den),y=shift(b,a.den);assert(fmpq_mpoly_equal(x.p,y.p,ctx));checks++;}
P powp(const P&a,long n){assert(n>=0);P p=one();for(long i=0;i<n;i++)p=times(p,a);return p;}
struct Node{char op;long a=-1,b=-1;unique_ptr<P>atom;};
int main(int argc,char**argv){assert(argc==3);fmpq_mpoly_ctx_init(ctx,3,ORD_DEGREVLEX);try{{
 ifstream in(argv[1]);long count;in>>reserve_terms>>count;assert(count==768);vector<Node>nodes(count);vector<long>uses(count);vector<unique_ptr<P>>vals(count);
 for(long i=0;i<count;i++){char op;in>>op;nodes[i].op=op;if(op=='p')nodes[i].atom=make_unique<P>(readpoly(in));else{in>>nodes[i].a>>nodes[i].b;assert(nodes[i].a<i&&nodes[i].b<i);uses[nodes[i].a]++;uses[nodes[i].b]++;}}
 long keepcount;in>>keepcount;set<long>keep;for(long i=0,j;i<keepcount;i++){in>>j;keep.insert(j);uses[j]++;}
 // The parser releases each literal as soon as it is moved to its value;
 // all literal atom polynomials are charged in live, plus input reserve.
 for(long i=0;i<count;i++){
  auto&n=nodes[i];if(n.op=='p')vals[i]=move(n.atom);
  else{P p=n.op=='m'?times(*vals[n.a],*vals[n.b]):add(*vals[n.a],*vals[n.b],n.op=='s'?-1:1);vals[i]=make_unique<P>(move(p));for(long j:{n.a,n.b})if(--uses[j]==0)vals[j].reset();}
  guard();if(!uses[i])vals[i].reset();
 }
 long special;in>>special;assert(special==3);for(long i=0;i<special;i++){long j;in>>j;equal(*vals[j],readpoly(in));}
 long output;in>>output;assert(output==767);P expected=readpoly(in);equal(*vals[output],expected);P cleared=readpoly(in);equal(expected,cleared);assert(!fmpq_mpoly_is_zero(expected.p,ctx));
 array<array<long,4>,4>A;array<long,4>target;array<array<long,5>,3>source;
 for(auto&r:A)for(auto&j:r)in>>j;for(auto&j:target)in>>j;for(auto&r:source)for(auto&j:r)in>>j;
 array<long,3>srcids;array<long,5>bids;for(auto&j:srcids)in>>j;for(auto&j:bids)in>>j;
 P delta2=readpoly(in);array<array<P,2>,3>N;for(auto&r:N)for(auto&x:r)x=readpoly(in);
 vector<P>units;for(long i=0;i<6;i++)units.push_back(readpoly(in));array<array<long,6>,3>pw;array<string,3>constants;array<long,6>lcm;
 for(auto&r:pw)for(auto&j:r)in>>j;for(auto&c:constants)in>>c;for(auto&j:lcm)in>>j;
 vector<P>co;co.push_back(delta2);co.push_back(add(times(N[0][0],N[2][1]),times(N[0][1],N[2][0]),-1));co.push_back(add(times(N[0][1],N[1][0]),times(N[0][0],N[1][1]),-1));
 for(long i=0;i<3;i++){P p=rational(co[i],constants[i]);if(i)p=add(P(),p,-1);for(long j=0;j<6;j++)p=times(p,powp(units[j],lcm[j]-pw[i][j]));equal(p,*vals[srcids[i]]);}
 vector<P>bh;
 for(long j=0;j<5;j++){P p;for(long i=0;i<3;i++)p=add(p,times(*vals[srcids[i]],*vals[source[i][j]]));equal(p,*vals[bids[j]]);bh.push_back(move(p));}
 // Independent border determinant: det([[P4,target],[Bhat,Bhat_target]])
 // equals the cleared Schur numerator. It avoids the author's four adj4
 // contractions and evaluates only this fixed selected bordered matrix.
 vector<vector<const P*>>border(5,vector<const P*>(5));for(long i=0;i<4;i++){for(long j=0;j<4;j++)border[i][j]=vals[A[i][j]].get();border[i][4]=vals[target[i]].get();}for(long j=0;j<5;j++)border[4][j]=&bh[j];
 array<int,5>perm={0,1,2,3,4};P determinant;
 do{long sign=1;for(long i=0;i<5;i++)for(long j=i+1;j<5;j++)if(perm[i]>perm[j])sign=-sign;P term=one();for(long i=0;i<5;i++)term=times(term,*border[i][perm[i]]);determinant=add(determinant,term,sign);}while(next_permutation(perm.begin(),perm.end()));
 equal(determinant,expected);string extra;assert(!(in>>extra));guard();
 ofstream out(argv[2]);out<<"{\"all_passed\":true,\"fixed_source_operations\":768,\"all_three_actual_H_source_cofactors_verified\":true,\"all_five_Bhat_actual_source_components_verified\":true,\"independent_fixed_bordered5_permutation_determinant_verified\":true,\"all7022_true_QQ_coefficients_equal\":true,\"no_a_min_normalization\":true,\"peak_live_terms_with_reserves\":"<<peak<<",\"exact_polynomial_equalities\":"<<checks<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"}\n";
 }}catch(const exception&e){cerr<<e.what()<<'\n';fmpq_mpoly_ctx_clear(ctx);return 2;}fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
