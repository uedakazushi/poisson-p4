#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <cassert>
#include <fstream>
#include <sstream>
#include <map>
#include <array>
#include <vector>
#include <string>
#include <iostream>
using namespace std;
struct P{fmpq_poly_t p;P(){fmpq_poly_init(p);}P(const P&o):P(){fmpq_poly_set(p,o.p);}P(P&&o):P(){fmpq_poly_swap(p,o.p);}P&operator=(P o){fmpq_poly_swap(p,o.p);return *this;}~P(){fmpq_poly_clear(p);}};
struct Term{array<int,7> mon;P num,den;};
void readpoly(istream&in,P&f){int n;in>>n;fmpq_t c;fmpq_init(c);for(int i=0;i<n;i++){string s;in>>s;assert(!fmpq_set_str(c,s.c_str(),10));fmpq_canonicalise(c);fmpq_poly_set_coeff_fmpq(f.p,i,c);}fmpq_clear(c);}
void exact(P&o,const P&a,const P&b){P rem;fmpq_poly_divrem(o.p,rem.p,a.p,b.p);assert(fmpq_poly_is_zero(rem.p));}
string pretty(const P&p){ostringstream o;bool first=true;fmpq_t c;fmpq_init(c);for(slong i=0;i<fmpq_poly_length(p.p);i++){fmpq_poly_get_coeff_fmpq(c,p.p,i);if(fmpq_is_zero(c))continue;char*s=fmpq_get_str(nullptr,10,c);if(!first)o<<"+";first=false;o<<"("<<s<<")";flint_free(s);if(i)o<<"*rr^"<<i;}fmpq_clear(c);return first?"0":o.str();}
int main(int argc,char**argv){assert(argc==4);ifstream raw(argv[1]),meta(argv[2]);ofstream out(argv[3]);assert(raw&&meta&&out);P ud,td;readpoly(meta,ud);readpoly(meta,td);int F,np,n,D;raw>>F>>np>>n>>D;assert(F==135&&np==1&&n==7&&D==10);vector<vector<Term>>rows(136);long count=0;
for(int d=2;d<=D;d++){P den,a,b;fmpq_poly_pow(a.p,ud.p,d-1);fmpq_poly_pow(b.p,td.p,d);fmpq_poly_mul(den.p,a.p,b.p);
 for(int kind=0;kind<2;kind++)for(int i=0;i<(kind?135:45);i++){long N;raw>>N;map<array<int,7>,P>terms;fmpq_t c;fmpq_init(c);
  for(long j=0;j<N;j++){string s;int e;array<int,7>key;raw>>s>>e;for(auto&z:key)raw>>z;assert(raw);if(!kind&&i!=44)continue;assert(!fmpq_set_str(c,s.c_str(),10));fmpq_canonicalise(c);fmpq_poly_set_coeff_fmpq(terms[key].p,e,c);}
  if(!kind&&d>=7)assert(N==0);
  if(kind||i==44)for(auto&[key,cp]:terms){P g;fmpq_poly_gcd(g.p,cp.p,den.p);Term t;t.mon=key;exact(t.num,cp,g);exact(t.den,den,g);fmpq_t lead;fmpq_init(lead);fmpq_poly_get_coeff_fmpq(lead,t.den.p,fmpq_poly_degree(t.den.p));fmpq_poly_scalar_div_fmpq(t.num.p,t.num.p,lead);fmpq_poly_scalar_div_fmpq(t.den.p,t.den.p,lead);fmpq_clear(lead);rows[kind?i:135].push_back(std::move(t));count++;}fmpq_clear(c);
 }
}
auto writepoly=[&](const P&p){out<<fmpq_poly_length(p.p);fmpq_t c;fmpq_init(c);for(slong j=0;j<fmpq_poly_length(p.p);j++){fmpq_poly_get_coeff_fmpq(c,p.p,j);char*v=fmpq_get_str(nullptr,10,c);out<<' '<<v;flint_free(v);}fmpq_clear(c);out<<'\n';};
out<<136<<'\n';int weights[7]={1,1,2,2,1,2,2};
for(int i=0;i<136;i++){
 P common;fmpq_poly_one(common.p);for(auto&t:rows[i]){P g,quot,product;fmpq_poly_gcd(g.p,common.p,t.den.p);exact(quot,t.den,g);fmpq_poly_mul(product.p,common.p,quot.p);common=std::move(product);}
 P content;for(auto&t:rows[i]){P mult,num,g;exact(mult,common,t.den);fmpq_poly_mul(num.p,t.num.p,mult.p);t.num=std::move(num);fmpq_poly_gcd(g.p,content.p,t.num.p);content=std::move(g);}
 int weight=-1;
 for(auto&t:rows[i]){P num;exact(num,t.num,content);t.num=std::move(num);int w=0;for(int k=0;k<7;k++)w+=weights[k]*t.mon[k];if(weight<0)weight=w;assert(weight==w);}
 out<<i<<' '<<weight<<' '<<rows[i].size()<<'\n';writepoly(common);writepoly(content);
 for(auto&t:rows[i]){for(auto z:t.mon)out<<z<<' ';out<<'\n';writepoly(t.num);}
 cerr<<"ROW "<<i<<" weight "<<weight<<" terms "<<rows[i].size()<<" lcd "<<fmpq_poly_degree(common.p)<<" content "<<fmpq_poly_degree(content.p)<<'\n';
}
cerr<<"NORMALIZED_COEFFICIENTS "<<count<<"\n";}
