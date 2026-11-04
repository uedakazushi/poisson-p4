// Independent fixed polynomial identities only: no gcd, solve or rank search.
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

long live_slots=0,peak_slots=0;constexpr long reserve_slots=44870,cap_slots=150000;
void need(bool b,const char*message){if(!b)throw std::runtime_error(message);}
struct P{
 fmpq_poly_t x;long counted=0;
 void sync(){long n=fmpq_poly_length(x);live_slots+=n-counted;counted=n;peak_slots=std::max(peak_slots,live_slots);need(live_slots+reserve_slots<=cap_slots,"150000 live polynomial/reserve slots");}
 P(){fmpq_poly_init(x);}P(const P&p):P(){fmpq_poly_set(x,p.x);sync();}
 P(P&&p)noexcept:P(){fmpq_poly_swap(x,p.x);counted=p.counted;p.counted=0;}
 P&operator=(const P&p){fmpq_poly_set(x,p.x);sync();return *this;}
 P&operator=(P&&p)noexcept{fmpq_poly_swap(x,p.x);std::swap(counted,p.counted);return *this;}
 ~P(){live_slots-=counted;fmpq_poly_clear(x);}
};
bool same(const P&a,const P&b){return fmpq_poly_equal(a.x,b.x);}
bool zero(const P&a){return fmpq_poly_is_zero(a.x);}
P constant(const std::string&s){P p;fmpq_t q;fmpq_init(q);need(!fmpq_set_str(q,s.c_str(),10),"rational input");fmpq_canonicalise(q);fmpq_poly_set_fmpq(p.x,q);fmpq_clear(q);p.sync();return p;}
P plus(const P&a,const P&b){P p;fmpq_poly_add(p.x,a.x,b.x);p.sync();return p;}
P times(const P&a,const P&b){P p;fmpq_poly_mul(p.x,a.x,b.x);p.sync();return p;}
P scaled(const P&a,const std::string&s){P p;fmpq_t q;fmpq_init(q);need(!fmpq_set_str(q,s.c_str(),10),"rational multiplier");fmpq_canonicalise(q);need(!fmpq_is_zero(q),"nonzero QQ unit");fmpq_poly_scalar_mul_fmpq(p.x,a.x,q);fmpq_clear(q);p.sync();return p;}
P power(const P&a,long e){need(e>=0&&e<1000,"power bound");P p;fmpq_poly_pow(p.x,a.x,e);p.sync();return p;}
P quotient(const P&a,const P&b){need(!zero(b),"nonzero exact divisor");P p,r;fmpq_poly_divrem(p.x,r.x,a.x,b.x);p.sync();r.sync();need(zero(r)&&same(times(p,b),a),"exact division/remultiplication");return p;}
P readpoly(){long n;need(bool(std::cin>>n)&&n>=0&&n<2000,"polynomial input count");P p;fmpq_t q;fmpq_init(q);for(long k=0;k<n;k++){long e;std::string c;need(bool(std::cin>>e>>c)&&e>=0&&e<2000,"polynomial term input");need(!fmpq_set_str(q,c.c_str(),10),"polynomial QQ input");fmpq_canonicalise(q);fmpq_poly_set_coeff_fmpq(p.x,e,q);}fmpq_clear(q);p.sync();return p;}
std::array<P,5> factors;
P unit_product(const std::array<int,5>&e){P p=constant("1");for(int i=0;i<5;i++)p=times(p,power(factors[i],e[i]));return p;}
struct Split{P full,primitive;};
Split readsplit(const P&actual){Split s;P full=readpoly();s.primitive=readpoly();P unit=readpoly();std::string qq;std::cin>>qq;std::array<int,5>e;for(int&i:e)std::cin>>i;
 need(same(full,actual),"computed minor equals frozen full polynomial");need(same(unit_product(e),unit),"five named factor multiplicities");
 need(same(scaled(times(unit,s.primitive),qq),actual),"full=QQunit*knownunits*primitive");
 if(!zero(actual)){
  fmpq_t lc;fmpq_init(lc);fmpq_poly_get_coeff_fmpq(lc,s.primitive.x,fmpq_poly_degree(s.primitive.x));need(fmpq_is_one(lc),"primitive is monic");fmpq_clear(lc);
  for(int i=0;i<5;i++){P q,r;fmpq_poly_divrem(q.x,r.x,s.primitive.x,factors[i].x);q.sync();r.sync();need(!zero(r),"maximal named-unit stripping");}
 }else need(zero(s.primitive),"zero minor primitive");
 s.full=std::move(full);return s;
}
int main(){try{
 auto start=std::chrono::steady_clock::now();
 for(P&f:factors)f=readpoly();P D=readpoly();std::vector<std::vector<P>>S(8,std::vector<P>(8)),H(8,std::vector<P>(8));
 for(auto&r:S)for(P&p:r)p=readpoly();for(auto&r:H)for(P&p:r)p=readpoly();
 for(int i=0;i<8;i++)for(int j=0;j<8;j++){
  P a,b;for(int k=0;k<8;k++){a=plus(a,times(S[i][k],H[k][j]));b=plus(b,times(H[i][k],S[k][j]));}
  need(same(a,i==j?D:P())&&same(b,i==j?D:P()),"all128 old adjugate products");
 }
 std::vector<std::array<int,5>>exponents(7);for(auto&e:exponents)for(int&v:e)std::cin>>v;
 std::vector<std::vector<P>>V(7,std::vector<P>(8));
 for(int i=0;i<7;i++){
  P lcd=unit_product(exponents[i]);
  for(int j=0;j<8;j++){P n=readpoly(),d=readpoly();P q=quotient(lcd,d);V[i][j]=times(n,q);need(same(times(V[i][j],d),times(n,lcd)),"all56 actual LCD RF clearings");}
 }
 std::vector<Split>splits;splits.push_back(readsplit(D));
 int count;std::cin>>count;need(count==15,"exactly15 fixed replacements");
 for(int k=0;k<count;k++){
  int light,row;std::cin>>light>>row;need(light==(k<8?0:1)&&row==(k<8?k:k-8),"fixed prefix replacement labels");
  P p,q;for(int j=0;j<8;j++)p=plus(p,times(V[light][j],H[j][row]));for(int j=7;j>=0;j--)q=plus(q,times(V[light][j],H[j][row]));
  need(same(p,q),"fixed forward/reverse8 dot");splits.push_back(readsplit(p));
 }
 P gcd=readpoly();need(fmpq_poly_is_one(gcd.x),"saved final gcd is1 (no gcd computed)");
 std::array<int,5>e;for(int&k:e)std::cin>>k;need(e==std::array<int,5>{21,93,6,0,1},"known E powers");P E=readpoly();need(same(E,unit_product(e)),"E is literal whole-P100 known unit");
 int n;std::cin>>n;need(n==5,"five primitive Bezout contributors");P sum;
 std::vector<int>indices;
 for(int j=0;j<n;j++){int i;std::cin>>i;need(i>=0&&i<16,"primitive contributor index");P a=readpoly();sum=plus(sum,times(a,splits[i].primitive));indices.push_back(i);}
 need(indices==std::vector<int>{0,7,9,10,15}&&fmpq_poly_is_one(sum.x),"literal sum primitive cofactors=1");
 std::cin>>n;need(n==5,"five actual-minor contributors");P actualsum;indices.clear();
 for(int j=0;j<n;j++){int i;std::cin>>i;need(i>=0&&i<16,"actual contributor index");P b=readpoly();actualsum=plus(actualsum,times(b,splits[i].full));indices.push_back(i);}
 need(indices==std::vector<int>{0,7,9,10,15}&&same(actualsum,E),"literal sum actual minors=E");
 std::string extra;need(!(std::cin>>extra),"trailing input");
 std::cout<<"{\"all_passed\":true,\"old_adjugate_products_checked\":128,\"actual_LCD_clearings_checked\":56,\"fixed_replacement_dots_checked\":15,\"minor_factorizations_checked\":16,\"primitive_Bezout_equals1\":true,\"actual_minor_Bezout_equals_E\":true,\"no_gcd_rank_search_or_unknown_factor_inverse\":true,\"peak_live_polynomial_slots\":"<<peak_slots<<",\"literal_reserve\":"<<reserve_slots<<",\"peak_with_reserve\":"<<peak_slots+reserve_slots<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
