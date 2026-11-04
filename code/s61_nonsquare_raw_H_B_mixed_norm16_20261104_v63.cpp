// ONE mixed-m1 determinant + primitive univariate GCD, no matrix-vector solve.
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
struct P{fmpz_poly_t v;P(){fmpz_poly_init(v);}P(const P&x):P(){fmpz_poly_set(v,x.v);}P(P&&x):P(){fmpz_poly_swap(v,x.v);}~P(){fmpz_poly_clear(v);}P&operator=(P x){fmpz_poly_swap(v,x.v);return*this;}};
struct Z{fmpz_t v;Z(){fmpz_init(v);}~Z(){fmpz_clear(v);}};
void readp(istream&f,P&p){long n;f>>n;assert(f&&n>=0&&n<=150000);Z z;vector<long>s;for(long i=0;i<n;i++){long d;string c;f>>d>>c;assert(f&&d>=0&&find(s.begin(),s.end(),d)==s.end());s.push_back(d);assert(!fmpz_set_str(z.v,c.c_str(),10)&&!fmpz_is_zero(z.v));fmpz_poly_set_coeff_fmpz(p.v,d,z.v);}}
long jsonbytes(const P&p){long n=2;Z z;for(long i=0;i<=fmpz_poly_degree(p.v);i++){fmpz_poly_get_coeff_fmpz(z.v,p.v,i);if(!fmpz_is_zero(z.v))n+=40+fmpz_sizeinbase(z.v,10);}return n;}
void jp(ostream&o,const P&p){o<<'[';bool first=true;Z z;for(long i=fmpz_poly_degree(p.v);i>=0;i--){fmpz_poly_get_coeff_fmpz(z.v,p.v,i);if(fmpz_is_zero(z.v))continue;if(!first)o<<',';first=false;char*s=fmpz_get_str(nullptr,10,z.v);o<<'['<<i<<",\""<<s<<"\"]";flint_free(s);}o<<']';}
void jz(ostream&o,const fmpz_t z){char*s=fmpz_get_str(nullptr,10,z);o<<'"'<<s<<'"';flint_free(s);}
int main(int argc,char**argv){
 assert(argc==3);auto start=chrono::steady_clock::now();ifstream in(argv[1]);vector<vector<P>>A(16,vector<P>(16));for(auto&r:A)for(auto&p:r)readp(in,p);P old;readp(in,old);assert(!fmpz_poly_is_zero(old.v));string extra;assert(!(in>>extra));
 P previous;fmpz_poly_one(previous.v);long divisions=0,peakdegree=0,peakbits=0;int sign=1;bool zero=false;
 for(int k=0;k<15;k++){
  assert(chrono::duration<double>(chrono::steady_clock::now()-start).count()<52);
  int row=15;while(row>=k&&fmpz_poly_is_zero(A[row][k].v))row--;if(row<k){zero=true;break;}
  if(row!=k){swap(A[row],A[k]);sign=-sign;}P pivot=A[k][k];
  for(int i=k+1;i<16;i++)for(int j=k+1;j<16;j++){
   P n,t,q,c;fmpz_poly_mul(n.v,pivot.v,A[i][j].v);fmpz_poly_mul(t.v,A[i][k].v,A[k][j].v);fmpz_poly_sub(n.v,n.v,t.v);assert(fmpz_poly_divides(q.v,n.v,previous.v));fmpz_poly_mul(c.v,q.v,previous.v);assert(fmpz_poly_equal(c.v,n.v));divisions++;
   peakdegree=max(peakdegree,(long)fmpz_poly_degree(n.v));peakbits=max(peakbits,labs((long)fmpz_poly_max_bits(n.v)));A[i][j]=std::move(q);
  }
  for(int i=k+1;i<16;i++)fmpz_poly_zero(A[i][k].v);previous=std::move(pivot);
 }
 P delta;if(!zero){delta=A[15][15];if(sign<0)fmpz_poly_neg(delta.v,delta.v);}zero=fmpz_poly_is_zero(delta.v);assert(fmpz_poly_degree(delta.v)<=152&&labs((long)fmpz_poly_max_bits(delta.v))<=5190);
 {ofstream o(string(argv[2])+".det-checkpoint.json");o<<"{\"determinant_zero\":"<<(zero?"true":"false")<<",\"exact_divisions\":"<<divisions<<",\"determinant\":";jp(o,delta);o<<"}\n";}
 if(zero){ofstream o(argv[2]);o<<"{\"determinant_zero\":true,\"stopped_before_gcd\":true,\"exact_divisions\":"<<divisions<<"}\n";return 0;}
 Z content;fmpz_poly_content(content.v,delta.v);assert(fmpz_sgn(content.v)>0);P primitive,res,c;fmpz_poly_scalar_divexact_fmpz(primitive.v,delta.v,content.v);fmpz_poly_scalar_mul_fmpz(c.v,primitive.v,content.v);assert(fmpz_poly_equal(c.v,delta.v));
 long amin=0;Z v;while(true){fmpz_poly_get_coeff_fmpz(v.v,primitive.v,amin);if(!fmpz_is_zero(v.v))break;amin++;}fmpz_poly_shift_right(res.v,primitive.v,amin);fmpz_poly_shift_left(c.v,res.v,amin);assert(fmpz_poly_equal(c.v,primitive.v));
 P g,qn,qo;fmpz_poly_gcd_subresultant(g.v,res.v,old.v);if(fmpz_sgn(fmpz_poly_lead(g.v))<0)fmpz_poly_neg(g.v,g.v);assert(fmpz_poly_divides(qn.v,res.v,g.v)&&fmpz_poly_divides(qo.v,old.v,g.v));fmpz_poly_mul(c.v,g.v,qn.v);assert(fmpz_poly_equal(c.v,res.v));fmpz_poly_mul(c.v,g.v,qo.v);assert(fmpz_poly_equal(c.v,old.v));
 {ofstream o(string(argv[2])+".gcd-checkpoint.json");o<<"{\"gcd\":";jp(o,g);o<<",\"new_primitive_without_a\":";jp(o,res);o<<",\"old_a117\":";jp(o,old);o<<",\"quotient_new\":";jp(o,qn);o<<",\"quotient_old\":";jp(o,qo);o<<"}\n";}
 long td=max(0L,(long)fmpz_poly_degree(qn.v))+max(0L,(long)fmpz_poly_degree(qo.v));long mb=max(labs((long)fmpz_poly_max_bits(qn.v)),labs((long)fmpz_poly_max_bits(qo.v)));long bb=(td+2)*(mb+20);long bytebound=(td+3)*(80+(bb*30103LL+99999)/100000);
 bool dobez=bytebound<24*1024*1024&&chrono::duration<double>(chrono::steady_clock::now()-start).count()<40;P S,T;Z R;bool bez=false;
 if(dobez){
  if(fmpz_poly_degree(qn.v)==0){fmpz_poly_get_coeff_fmpz(R.v,qn.v,0);fmpz_poly_one(S.v);}
  else if(fmpz_poly_degree(qo.v)==0){fmpz_poly_get_coeff_fmpz(R.v,qo.v,0);fmpz_poly_one(T.v);}
  else fmpz_poly_xgcd_modular(R.v,S.v,T.v,qn.v,qo.v);
  assert(!fmpz_is_zero(R.v));P x,y,r;fmpz_poly_mul(x.v,S.v,qn.v);fmpz_poly_mul(y.v,T.v,qo.v);fmpz_poly_add(x.v,x.v,y.v);fmpz_poly_set_fmpz(r.v,R.v);assert(fmpz_poly_equal(x.v,r.v));
  fmpz_poly_mul(x.v,S.v,res.v);fmpz_poly_mul(y.v,T.v,old.v);fmpz_poly_add(x.v,x.v,y.v);fmpz_poly_scalar_mul_fmpz(r.v,g.v,R.v);assert(fmpz_poly_equal(x.v,r.v));bez=true;
 }
 long base=jsonbytes(delta)+jsonbytes(primitive)+jsonbytes(res)+jsonbytes(old)+jsonbytes(g)+jsonbytes(qn)+jsonbytes(qo)+10000;bool savebez=bez&&base+jsonbytes(S)+jsonbytes(T)+fmpz_sizeinbase(R.v,10)<28*1024*1024;assert(base<28*1024*1024);
 ofstream o(argv[2]);o<<"{\"determinant_zero\":false,\"determinant\":";jp(o,delta);o<<",\"determinant_QQ_content\":";jz(o,content.v);o<<",\"primitive\":";jp(o,primitive);o<<",\"removed_known_a_power\":"<<amin<<",\"new_primitive_without_a\":";jp(o,res);o<<",\"old_primitive_a117\":";jp(o,old);o<<",\"gcd\":";jp(o,g);o<<",\"quotient_new\":";jp(o,qn);o<<",\"quotient_old\":";jp(o,qo);
 o<<",\"Bezout_computed\":"<<(bez?"true":"false")<<",\"Bezout_saved\":"<<(savebez?"true":"false")<<",\"Bezout_size_bound_bytes\":"<<bytebound;
 if(savebez){o<<",\"Bezout_S\":";jp(o,S);o<<",\"Bezout_T\":";jp(o,T);o<<",\"Bezout_constant\":";jz(o,R.v);}
 o<<",\"exact_Bareiss_divisions\":"<<divisions<<",\"maximum_numerator_degree\":"<<peakdegree<<",\"maximum_numerator_bits\":"<<peakbits<<",\"gcd_degree\":"<<fmpz_poly_degree(g.v)<<",\"gcd_is_QQunit\":"<<(fmpz_poly_degree(g.v)==0?"true":"false")<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";assert(o.tellp()<28*1024*1024);cerr<<"DETERMINANT_GCD_PRODUCT_IDENTITIES_PASS\n";
}
