// Independent fixed subset-Laplace determinant and literal certificate checker.
// No author helper, elimination, gcd, factorization, point or rank search.
#include <flint/fmpq_poly.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/nmod_poly.h>
#include <flint/ulong_extras.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
static long live=0, peak=0, reserve_slots=51034, denominator_checks=0;
static auto began=chrono::steady_clock::now();
static constexpr ulong prime=32003;
struct P {
 fmpq_poly_t v; long slots=0;
 P(){fmpq_poly_init(v);}
 P(const P&a):P(){fmpq_poly_set(v,a.v);update();}
 P(P&&a):P(){fmpq_poly_swap(v,a.v);slots=a.slots;a.slots=0;}
 P&operator=(P a){fmpq_poly_swap(v,a.v);swap(slots,a.slots);return *this;}
 ~P(){live-=slots;fmpq_poly_clear(v);}
 void update(){live+=fmpq_poly_length(v)-slots;slots=fmpq_poly_length(v);peak=max(peak,live);
  if(peak+reserve_slots>150000)throw runtime_error("150000 live coefficient slots exceeded");
  if(chrono::duration<double>(chrono::steady_clock::now()-began).count()>53)throw runtime_error("53 second arithmetic guard");}
};
P readpoly(istream&f){long n;f>>n;assert(n>=0);P p;fmpq_t q;fmpq_init(q);long last=-1;
 for(long i=0;i<n;i++){long e;string s;f>>e>>s;assert(e>last);last=e;
  assert(!fmpq_set_str(q,s.c_str(),10));fmpq_canonicalise(q);assert(!fmpq_is_zero(q));fmpq_poly_set_coeff_fmpq(p.v,e,q);}
 fmpq_clear(q);p.update();return p;}
P mul(const P&a,const P&b){P c;fmpq_poly_mul(c.v,a.v,b.v);c.update();return c;}
P add(const P&a,const P&b,int sign=1){P c;if(sign==1)fmpq_poly_add(c.v,a.v,b.v);else fmpq_poly_sub(c.v,a.v,b.v);c.update();return c;}
P one(){P a;fmpq_poly_one(a.v);a.update();return a;}
P linear(long c0,long c1=1){P a;fmpq_poly_set_coeff_si(a.v,0,c0);fmpq_poly_set_coeff_si(a.v,1,c1);a.update();return a;}
bool eq(const P&a,const P&b){return fmpq_poly_equal(a.v,b.v);}
void p_integral(const P&a){fmpq_t q;fmpq_init(q);for(long i=0;i<fmpq_poly_length(a.v);i++){
 fmpq_poly_get_coeff_fmpq(q,a.v,i);assert(fmpz_fdiv_ui(fmpq_denref(q),prime));denominator_checks++;}fmpq_clear(q);}
void lc_p_unit(const P&a){fmpq_t q;fmpq_init(q);fmpq_poly_get_coeff_fmpq(q,a.v,fmpq_poly_degree(a.v));
 assert(fmpz_fdiv_ui(fmpq_denref(q),prime)&&fmpz_fdiv_ui(fmpq_numref(q),prime));fmpq_clear(q);}
P exact_divide(const P&a,const P&b){P q,r;assert(!fmpq_poly_is_zero(b.v));fmpq_poly_divrem(q.v,r.v,a.v,b.v);q.update();r.update();
 assert(fmpq_poly_is_zero(r.v));assert(eq(mul(q,b),a));return q;}
struct Mod {nmod_poly_t v;Mod(){nmod_poly_init(v,prime);}~Mod(){nmod_poly_clear(v);}};
void reduce(const P&a,Mod&b){fmpq_t q;fmpq_init(q);for(long i=0;i<fmpq_poly_length(a.v);i++){
 fmpq_poly_get_coeff_fmpq(q,a.v,i);ulong d=fmpz_fdiv_ui(fmpq_denref(q),prime);assert(d);
 // Both factors are <32003, so this ordinary product fits in ulong.
 ulong c=(fmpz_fdiv_ui(fmpq_numref(q),prime)*n_invmod(d,prime))%prime;
 nmod_poly_set_coeff_ui(b.v,i,c);}fmpq_clear(q);}
int main(int argc,char**argv){assert(argc==2);try{
 ifstream f(argv[1]);size_t count;f>>count;vector<P> pool;pool.reserve(count);for(size_t i=0;i<count;i++)pool.push_back(readpoly(f));
 vector<vector<int>> ex(9,vector<int>(5));for(auto&r:ex)for(int&x:r){f>>x;assert(x>=0);}
 vector<vector<pair<int,int>>> refs(9,vector<pair<int,int>>(9));for(auto&r:refs)for(auto&x:r)f>>x.first>>x.second;
 vector<int> id(13);for(int&x:id)f>>x;
 vector<vector<int>> old(5,vector<int>(3));for(auto&r:old)for(int&x:r)f>>x;
 string extra;assert(!(f>>extra));
 const P&N=pool[id[0]],&H=pool[id[1]],&Q=pool[id[2]],&R=pool[id[3]],&PH=pool[id[4]],&PR=pool[id[5]],
  &CH=pool[id[6]],&CR=pool[id[7]],&G=pool[id[8]],&S=pool[id[9]],&T=pool[id[10]],&QH=pool[id[11]],&QR=pool[id[12]];
 vector<P> factors;factors.push_back(linear(0));factors.push_back(linear(20));factors.push_back(linear(60));
 factors.push_back(linear(2700,203));factors.push_back(linear(12));
 vector<vector<P>> A(9,vector<P>(9));long lcd_checks=0;
 for(int i=0;i<9;i++){P lcd=one();for(int k=0;k<5;k++)for(int j=0;j<ex[i][k];j++)lcd=mul(lcd,factors[k]);
  for(int j=0;j<9;j++){const P&num=pool[refs[i][j].first],&den=pool[refs[i][j].second];p_integral(num);p_integral(den);
   A[i][j]=mul(num,exact_divide(lcd,den));assert(eq(mul(A[i][j],den),mul(num,lcd)));p_integral(A[i][j]);lcd_checks++;}}
 p_integral(H);p_integral(N);p_integral(Q);p_integral(R);
 // Determinant by fixed column-set dynamic Laplace expansion. All subsets of
 // the fixed nine columns, no pivoting or division. Only adjacent layers live.
 vector<P> dp(512);dp[0]=one();long products=0;
 for(int row=0;row<9;row++){vector<P> next(512);
  for(int mask=0;mask<512;mask++)if(__builtin_popcount((unsigned)mask)==row){
   if(fmpq_poly_is_zero(dp[mask].v))continue;
   for(int col=0;col<9;col++)if(!(mask&(1<<col))){int to=mask|(1<<col);
    int sign=(__builtin_popcount((unsigned)mask>>(col+1))&1)?-1:1;
    next[to]=add(next[to],mul(dp[mask],A[row][col]),sign);products++;}}
  dp=move(next);}
 assert(eq(dp[511],N));assert(fmpq_poly_degree(N.v)==202);dp.clear();
 assert(eq(add(mul(H,Q),R),N));assert(fmpq_poly_degree(H.v)==15&&fmpq_poly_degree(Q.v)==187&&fmpq_poly_degree(R.v)==14);
 assert(eq(mul(CH,PH),H)&&eq(mul(CR,PR),R));assert(!fmpq_poly_is_zero(CH.v)&&!fmpq_poly_is_zero(CR.v));
 assert(eq(add(mul(S,PH),mul(T,PR)),G));assert(eq(mul(G,QH),H)&&eq(mul(G,QR),R));
 P expected;fmpq_t q;fmpq_init(q);fmpq_set_si(q,22800,13);fmpq_poly_set_coeff_fmpq(expected.v,0,q);
 fmpq_set_si(q,200,13);fmpq_poly_set_coeff_fmpq(expected.v,1,q);fmpq_poly_set_coeff_si(expected.v,2,1);expected.update();fmpq_clear(q);assert(eq(G,expected));
 // Existing five-minor common H certificate: literal divisions and saved
 // old-prime cofactor Bezout, with p-integrality and degree preservation.
 Mod sum,a,b,prod;
 for(int i=0;i<5;i++){const P&primitive=pool[old[i][0]],&quot=pool[old[i][1]],&bezout=pool[old[i][2]];
  assert(eq(mul(H,quot),primitive));p_integral(quot);lc_p_unit(quot);reduce(quot,a);reduce(bezout,b);
  nmod_poly_mul(prod.v,a.v,b.v);nmod_poly_add(sum.v,sum.v,prod.v);}
 assert(nmod_poly_is_one(sum.v));
 cout<<"{\"all_passed\":true,\"all81_actual_source_LCD_checks\":"<<lcd_checks
  <<",\"fixed_subset_Laplace_products\":"<<products<<",\"determinant_degree\":202,\"exact_N202_H15_Q187_R14\":true"
  <<",\"QQ_gcd_from_literal_Bezout_and_both_divisions\":true,\"exact_gcd\":\"z^2+200z/13+22800/13\""
  <<",\"old_five_H_divisions_and_saved_modular_cofactor_Bezout\":true,\"p_integral_coefficient_checks\":"<<denominator_checks
  <<",\"no_elimination_gcd_factor_or_rank_search\":true,\"literal_reserve\":"<<reserve_slots
  <<",\"peak_live_slots_including_reserve\":"<<peak+reserve_slots
  <<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";
 }catch(const exception&e){cerr<<e.what()<<'\n';return 2;}flint_cleanup();}
