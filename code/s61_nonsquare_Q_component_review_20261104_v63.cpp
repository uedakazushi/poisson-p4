#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <gmp.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>

// No field division, inverse, gcd or Euclid routine is called. All identities
// use saved coefficients and unreduced QQ[e] numerator/denominator pairs.
struct P {
 fmpq_poly_t x;
 P(){fmpq_poly_init(x);} explicit P(long n):P(){fmpq_poly_set_si(x,n);}
 P(const P& p):P(){fmpq_poly_set(x,p.x);} P(P&&p):P(){fmpq_poly_swap(x,p.x);}
 ~P(){fmpq_poly_clear(x);} P&operator=(const P&p){fmpq_poly_set(x,p.x);return *this;}
 bool zero()const{return fmpq_poly_is_zero(x);} bool operator==(const P&p)const{return fmpq_poly_equal(x,p.x);}
};
P operator+(const P&a,const P&b){P r;fmpq_poly_add(r.x,a.x,b.x);return r;}
P operator-(const P&a,const P&b){P r;fmpq_poly_sub(r.x,a.x,b.x);return r;}
P operator*(const P&a,const P&b){P r;fmpq_poly_mul(r.x,a.x,b.x);return r;}
P readP(std::istream&f){int n;f>>n;assert(f.good());P r;fmpq_t c;fmpq_init(c);for(int i=0;i<n;i++){int e;std::string v;f>>e>>v;assert(e>=0);assert(fmpq_set_str(c,v.c_str(),10)==0);fmpq_canonicalise(c);fmpq_poly_set_coeff_fmpq(r.x,e,c);}fmpq_clear(c);return r;}
P q0,q1;long peak=0;
struct E {
 P c0,c1,den;
 E():c0(0),c1(0),den(1){}
 E(P a,P b,P d):c0(std::move(a)),c1(std::move(b)),den(std::move(d)){
  assert(!den.zero());if(c0.zero()&&c1.zero())den=P(1);
  peak=std::max(peak,std::max({fmpq_poly_degree(c0.x),fmpq_poly_degree(c1.x),fmpq_poly_degree(den.x)}));
 }
 bool zero()const{return c0.zero()&&c1.zero();}
};
E operator+(const E&a,const E&b){if(a.zero())return b;if(b.zero())return a;if(a.den==b.den)return E(a.c0+b.c0,a.c1+b.c1,a.den);return E(a.c0*b.den+b.c0*a.den,a.c1*b.den+b.c1*a.den,a.den*b.den);}
E operator-(const E&a){return E(P(0)-a.c0,P(0)-a.c1,a.den);}
E operator*(const E&a,const E&b){if(a.zero()||b.zero())return E();P c=a.c1*b.c1;return E(a.c0*b.c0-q0*c,a.c0*b.c1+a.c1*b.c0-q1*c,a.den*b.den);}
bool eq(const E&a,const E&b){return a.c0*b.den==b.c0*a.den && a.c1*b.den==b.c1*a.den;}
E conj(const E&a){return E(a.c0-q1*a.c1,P(0)-a.c1,a.den);}
E readE(std::istream&f){int n;f>>n;E r,A(P(0),P(1),P(1));for(int j=0;j<n;j++){int a;f>>a;assert(a>=0&&a<=16);P nu=readP(f),de=readP(f);E t(nu,P(0),de);for(int i=0;i<a;i++)t=t*A;r=r+t;}return r;}
using D=std::vector<E>;
void trim(D&a){while(!a.empty()&&a.back().zero())a.pop_back();}
D readD(std::istream&f){int n;f>>n;D a;for(int i=0;i<n;i++){int k;f>>k;assert(k>=0&&k<40);if((int)a.size()<=k)a.resize(k+1);a[k]=a[k]+readE(f);}trim(a);return a;}
D add(const D&a,const D&b){D r(std::max(a.size(),b.size()));for(size_t i=0;i<r.size();i++)r[i]=(i<a.size()?a[i]:E())+(i<b.size()?b[i]:E());trim(r);return r;}
D neg(const D&a){D r;for(auto&x:a)r.push_back(-x);return r;}
D mul(const D&a,const D&b){if(a.empty()||b.empty())return {};D r(a.size()+b.size()-1);for(size_t i=0;i<a.size();i++)for(size_t j=0;j<b.size();j++)r[i+j]=r[i+j]+a[i]*b[j];trim(r);return r;}
bool eq(const D&a,const D&b){for(size_t i=0;i<std::max(a.size(),b.size());i++)if(!eq(i<a.size()?a[i]:E(),i<b.size()?b[i]:E()))return false;return true;}
int main(int argc,char**argv){assert(argc==2);std::ifstream f(argv[1]);q0=readP(f);q1=readP(f);int n;f>>n;std::vector<D>v;long checks=0;for(int z=0;z<n;z++){
 std::string op;f>>op;assert(f.good());int i,j;
 if(op=="P")v.push_back(readD(f));
 else if(op=="A"||op=="M"){f>>i>>j;assert(i<(int)v.size()&&j<(int)v.size());v.push_back(op=="A"?add(v[i],v[j]):mul(v[i],v[j]));}
 else if(op=="N"){f>>i;v.push_back(neg(v.at(i)));}
 else if(op=="LC"){f>>i;assert(!v.at(i).empty());v.push_back(D{v[i].back()});}
 else if(op=="C"){f>>i;assert(v.at(i).size()<=1);v.push_back(v[i].empty()?D{}:D{conj(v[i][0])});}
 else if(op=="EV"){f>>i>>j;assert(v.at(j).size()<=1);D r;for(int k=(int)v.at(i).size()-1;k>=0;k--)r=add(mul(r,v[j]),D{v[i][k]});v.push_back(r);}
 else if(op=="E"){f>>i>>j;assert(eq(v.at(i),v.at(j)));checks++;}
 else if(op=="NZ"){f>>i;assert(!v.at(i).empty());checks++;}
 else if(op=="DEG"){f>>i>>j;assert((int)v.at(i).size()-1==j);checks++;}
 else{std::cerr<<"unknown operation "<<op<<"\n";return 3;}
 }
 std::cout<<"{\"all_passed\":true,\"exact_checks\":"<<checks<<",\"registers\":"<<v.size()<<",\"peak_unreduced_e_degree\":"<<peak<<",\"field_inverse_or_GCD_calls\":0}\n";
}
