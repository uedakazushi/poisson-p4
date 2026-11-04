// Exact univariate source arithmetic for the whole P80 curve.
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
using namespace std;auto began=chrono::steady_clock::now();long peak=0,stored=0;
struct P{fmpq_poly_t v;P(){fmpq_poly_init(v);}P(const P&a):P(){fmpq_poly_set(v,a.v);}P(P&&a):P(){fmpq_poly_swap(v,a.v);}P&operator=(P a){fmpq_poly_swap(v,a.v);return*this;}~P(){fmpq_poly_clear(v);}};
void guard(const P&p){peak=max(peak,fmpq_poly_length(p.v));assert(peak<=5000);assert(chrono::duration<double>(chrono::steady_clock::now()-began).count()<50);}
P pc(const string&s){P p;fmpq_t q;fmpq_init(q);assert(!fmpq_set_str(q,s.c_str(),10));fmpq_canonicalise(q);fmpq_poly_set_fmpq(p.v,q);fmpq_clear(q);return p;}
P one(){return pc("1");}P plusp(const P&a,const P&b){P p;fmpq_poly_add(p.v,a.v,b.v);guard(p);return p;}P negp(const P&a){P p;fmpq_poly_neg(p.v,a.v);return p;}P mulp(const P&a,const P&b){P p;fmpq_poly_mul(p.v,a.v,b.v);guard(p);return p;}P scalar(const P&a,const string&s){P p;fmpq_t q;fmpq_init(q);assert(!fmpq_set_str(q,s.c_str(),10));fmpq_canonicalise(q);fmpq_poly_scalar_mul_fmpq(p.v,a.v,q);fmpq_clear(q);guard(p);return p;}
P shiftp(const P&a,long k){assert(k>=0);P p;fmpq_poly_shift_left(p.v,a.v,k);guard(p);return p;}P quotient(const P&a,const P&b){P q,r;fmpq_poly_divrem(q.v,r.v,a.v,b.v);assert(fmpq_poly_is_zero(r.v));assert(fmpq_poly_equal(mulp(q,b).v,a.v));return q;}
struct R{P n,d;R():d(one()){}R(P x):n(x),d(one()){}R(P x,P y):n(x),d(y){norm();}void norm(){assert(!fmpq_poly_is_zero(d.v));if(fmpq_poly_is_zero(n.v)){d=one();return;}P g;fmpq_poly_gcd(g.v,n.v,d.v);n=quotient(n,g);d=quotient(d,g);fmpq_t c;fmpq_init(c);fmpq_poly_get_coeff_fmpq(c,d.v,fmpq_poly_degree(d.v));fmpq_poly_scalar_div_fmpq(n.v,n.v,c);fmpq_poly_scalar_div_fmpq(d.v,d.v,c);fmpq_clear(c);guard(n);guard(d);}};
R rc(const string&s){return R(pc(s));}R add(const R&a,const R&b){return R(plusp(mulp(a.n,b.d),mulp(b.n,a.d)),mulp(a.d,b.d));}R neg(const R&a){return R(negp(a.n),a.d);}R sub(const R&a,const R&b){return add(a,neg(b));}R mul(const R&a,const R&b){return R(mulp(a.n,b.n),mulp(a.d,b.d));}bool zero(const R&a){return fmpq_poly_is_zero(a.n.v);}bool eq(const R&a,const R&b){return fmpq_poly_equal(a.n.v,b.n.v)&&fmpq_poly_equal(a.d.v,b.d.v);}map<string,long>inversions;
R divunit(const R&a,const R&b,const string&why){assert(!zero(b));inversions[why]++;return R(mulp(a.n,b.d),mulp(a.d,b.n));}
R powr(R p,long k){if(k<0){p=divunit(rc("1"),p,"certified_FULLU_unit_negative_power");k=-k;}R r=rc("1");for(long j=0;j<k;j++)r=mul(r,p);return r;}
struct T{string q;long a,e,d;};using C=vector<T>;
C readc(istream&f){long n;f>>n;assert(n>=0);C p;for(long i=0;i<n;i++){T t;f>>t.q>>t.a>>t.e>>t.d;assert(t.a>=0&&t.e>=0&&t.d>=0);p.push_back(t);}return p;}
vector<P>ap,ep,dp;P z;
R eval(const C&p){long clear=0;for(auto&t:p)clear=max(clear,2*t.a+t.e);P r;for(auto&t:p){while(long(ap.size())<=t.a)ap.push_back(mulp(ap.back(),ap[1]));while(long(ep.size())<=t.e)ep.push_back(mulp(ep.back(),ep[1]));while(long(dp.size())<=t.d)dp.push_back(mulp(dp.back(),dp[1]));P q=scalar(shiftp(mulp(mulp(ap[t.a],ep[t.e]),dp[t.d]),clear-2*t.a-t.e),t.q);r=plusp(r,q);}return R(r,shiftp(one(),clear));}
void jp(ostream&o,const P&p){stored+=fmpq_poly_length(p.v);assert(stored<=150000);o<<'[';fmpq_t q;fmpq_init(q);bool first=true;for(long i=0;i<fmpq_poly_length(p.v);i++){fmpq_poly_get_coeff_fmpq(q,p.v,i);if(fmpq_is_zero(q))continue;if(!first)o<<',';first=false;char*s=fmpq_get_str(nullptr,10,q);o<<"[["<<i<<"],\""<<s<<"\"]";flint_free(s);}fmpq_clear(q);o<<']';}
void jr(ostream&o,const R&r){o<<"{\"numerator\":";jp(o,r.n);o<<",\"denominator\":";jp(o,r.d);o<<'}';}
P unitproduct;
bool unit_numerator(const R&r){if(zero(r))return false;P p=r.n;for(int step=0;step<100;step++){P g;fmpq_poly_gcd(g.v,p.v,unitproduct.v);if(fmpq_poly_degree(g.v)==0)break;p=quotient(p,g);}return fmpq_poly_degree(p.v)==0;}
int main(int argc,char**argv){assert(argc==3);began=chrono::steady_clock::now();{
 ifstream f(argv[1]);z=shiftp(one(),1);P zp20=plusp(z,pc("20"));ap={one(),scalar(plusp(scalar(z,"7"),pc("80")),"-8/45")};ep={one(),scalar(plusp(z,pc("40")),"2/9")};dp={one(),scalar(mulp(z,zp20),"-3/80")};R a(ap[1],shiftp(one(),2)),e(ep[1],z),d(dp[1]);
 R Cc=sub(sub(add(mul(rc("100"),a),mul(rc("18"),powr(e,2))),mul(rc("9"),e)),rc("2"));
 R F=sub(add(add(mul(rc("500"),a),mul(rc("90"),powr(e,2))),mul(rc("3"),e)),rc("24"));
 R Bb=add(mul(sub(sub(add(mul(rc("50"),a),mul(rc("18"),powr(e,2))),mul(rc("9"),e)),rc("2")),d),sub(mul(rc("60"),e),rc("40")));
 R J=sub(sub(mul(rc("5"),mul(a,d)),mul(rc("6"),e)),rc("4"));R g1=add(mul(sub(mul(rc("3"),e),rc("2")),d),rc("20"));R g2=add(mul(add(mul(rc("6"),e),rc("1")),d),rc("40"));
 assert(eq(powr(Bb,2),mul(Cc,mul(g1,g2))));array<R,6>units={a,d,F,Bb,Cc,J};unitproduct=mulp(mulp(mulp(mulp(mulp(z,zp20),plusp(scalar(z,"7"),pc("300"))),plusp(scalar(z,"7"),pc("80"))),plusp(scalar(z,"17"),pc("216"))),plusp(scalar(z,"7"),pc("60")));
 vector<vector<R>>raw(16,vector<R>(7));long n,m;f>>n>>m;assert(n==16&&m==7);for(auto&row:raw)for(auto&v:row)v=eval(readc(f));
 vector<int>ids(5);for(auto&i:ids)f>>i;assert(ids==vector<int>({1,7,8,10,13}));map<int,vector<R>>M;for(int i:ids){M[i]=vector<R>(6);for(auto&v:M[i])v=eval(readc(f));}
 R Delta=eval(readc(f));assert(unit_numerator(Delta));vector<vector<long>>H(7,vector<long>(7));for(auto&row:H)for(auto&v:row)f>>v;
 vector<array<long,8>>mons(29);for(auto&v:mons)for(auto&x:v)f>>x;vector<vector<C>>quad(6,vector<C>(29));for(auto&row:quad)for(auto&p:row)p=readc(f);
 string extra;assert(!(f>>extra));

 // The explicit source-unit bridge to the missing functional.
 R B4=powr(Bb,4),N=mul(rc("37500"),mul(powr(Cc,2),g1));assert(eq(mul(M[8][0],divunit(N,B4,"row8_inverse_B4")),rc("1")));
 vector<int>pr={7,13,10},pcv={2,4,3};vector<vector<R>>A(3,vector<R>(3)),adj(3,vector<R>(3));vector<R>S(3);
 for(int i=0;i<3;i++){S[i]=sub(mul(B4,M[pr[i]][5]),mul(mul(M[pr[i]][0],N),M[8][5]));for(int j=0;j<3;j++)A[i][j]=M[pr[i]][pcv[j]];}
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){vector<int>rr,cc;for(int k=0;k<3;k++){if(k!=j)rr.push_back(k);if(k!=i)cc.push_back(k);}R v=sub(mul(A[rr[0]][cc[0]],A[rr[1]][cc[1]]),mul(A[rr[0]][cc[1]],A[rr[1]][cc[0]]));adj[i][j]=(i+j)%2?neg(v):v;}
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){R v,w;for(int k=0;k<3;k++){v=add(v,mul(A[i][k],adj[k][j]));w=add(w,mul(adj[i][k],A[k][j]));}assert(eq(v,i==j?Delta:R())&&eq(w,i==j?Delta:R()));}
 vector<R>w(3);for(int i=0;i<3;i++){for(int j=0;j<3;j++)w[i]=add(w[i],mul(adj[i][j],S[j]));w[i]=divunit(w[i],Delta,"proven_global_Theta_Delta_unit");}
 vector<R>eta6(6);eta6[0]=neg(mul(N,M[8][5]));eta6[2]=neg(w[0]);eta6[4]=neg(w[1]);eta6[3]=neg(w[2]);eta6[5]=B4;for(int j:{0,2,3,4,5})eta6[1]=sub(eta6[1],mul(M[1][j],eta6[j]));
 assert(eq(M[1][1],mul(rc("-1/14062500"),mul(Bb,Cc))));eta6[1]=divunit(eta6[1],M[1][1],"saved_unnormalized_checkpoint_pivot_minus_BC_over14062500");
 vector<R>etaX(7),eta(7);for(int j=0;j<5;j++)etaX[j]=eta6[j];R hh=sub(mul(Cc,g1),Bb);etaX[5]=mul(hh,B4);etaX[6]=neg(mul(mul(a,mul(g1,g2)),B4));
 for(int p=0;p<7;p++)for(int j=0;j<7;j++)eta[p]=add(eta[p],mul(rc(to_string(H[j][p])),etaX[j]));assert(eq(eta[3],etaX[6])&&unit_numerator(eta[3]));R norm=eta[3];for(auto&v:eta)v=divunit(v,norm,"known_unit_physical_v3_normalization");for(auto&v:etaX)v=divunit(v,norm,"known_unit_physical_v3_normalization");assert(eq(eta[3],rc("1")));
 for(int i=0;i<16;i++){R v;for(int j=0;j<7;j++)v=add(v,mul(raw[i][j],eta[j]));assert(zero(v));}
 vector<R>v8(8);vector<int>physical={0,1,2,3,4,6,7};for(int j=0;j<7;j++)v8[physical[j]]=eta[j];vector<R>qe(6);
 for(int i=0;i<6;i++)for(int k=0;k<29;k++){if(mons[k][5]){assert(quad[i][k].empty());continue;}R t=eval(quad[i][k]);for(int j=0;j<8;j++)if(mons[k][j])t=mul(t,powr(v8[j],mons[k][j]));qe[i]=add(qe[i],t);}
 bool descent=true;for(int i=0;i<5;i++)descent&=zero(qe[i]);
 ofstream o(argv[2]);o<<"{\"curve\":\"P80\",\"eta_physical7\":[";for(int j=0;j<7;j++){if(j)o<<',';jr(o,eta[j]);}o<<"],\"eta_oldX7\":[";for(int j=0;j<7;j++){if(j)o<<',';jr(o,etaX[j]);}o<<"],\"v3_is_one\":true,\"all16_actual_rawH_annihilations_exact\":true,\"quadratic_evaluations\":[";for(int i=0;i<6;i++){if(i)o<<',';jr(o,qe[i]);}o<<"],\"quadratic_zero_flags\":[";for(int i=0;i<6;i++){if(i)o<<',';o<<(zero(qe[i])?"true":"false");}o<<"],\"quadratic_unit_flags\":[";for(int i=0;i<6;i++){if(i)o<<',';o<<(unit_numerator(qe[i])?"true":"false");}o<<"],\"all5_eliminated_Hquad_directions_annihilated\":"<<(descent?"true":"false")<<",\"no_cache_row_query\":true,\"peak_polynomial_length\":"<<peak<<",\"stored_output_terms\":"<<stored<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";
 }flint_cleanup();}
