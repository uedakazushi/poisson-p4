// Independent Brown PRS: five-variable reversed LEX coefficient context.
// Deliberately does not call resultant, gcd, interpolation, or a finite field.
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_mpoly.h>
#include <openssl/evp.h>
using namespace std;
fmpz_mpoly_ctx_t ctx;
struct P {
 fmpz_mpoly_t a;
 P(){fmpz_mpoly_init(a,ctx);} P(long n):P(){fmpz_mpoly_set_si(a,n,ctx);}
 P(const P&b):P(){fmpz_mpoly_set(a,b.a,ctx);} P(P&&b):P(){fmpz_mpoly_swap(a,b.a,ctx);}
 P&operator=(const P&b){if(this!=&b)fmpz_mpoly_set(a,b.a,ctx);return *this;}
 P&operator=(P&&b){if(this!=&b)fmpz_mpoly_swap(a,b.a,ctx);return *this;}
 ~P(){fmpz_mpoly_clear(a,ctx);}
};
P operator+(const P&a,const P&b){P c;fmpz_mpoly_add(c.a,a.a,b.a,ctx);return c;}
P operator-(const P&a,const P&b){P c;fmpz_mpoly_sub(c.a,a.a,b.a,ctx);return c;}
P operator*(const P&a,const P&b){P c;fmpz_mpoly_mul(c.a,a.a,b.a,ctx);return c;}
P neg(const P&a){P c;fmpz_mpoly_neg(c.a,a.a,ctx);return c;}
P pw(const P&a,unsigned n){P c;assert(fmpz_mpoly_pow_ui(c.a,a.a,n,ctx));return c;}
bool zero(const P&a){return fmpz_mpoly_is_zero(a.a,ctx);}
bool equal(const P&a,const P&b){return fmpz_mpoly_equal(a.a,b.a,ctx);}
long terms(const P&a){return fmpz_mpoly_length(a.a,ctx);}
string termline(const P&a,long i){fmpz_t c;fmpz_init(c);ulong e[5];fmpz_mpoly_get_term_coeff_fmpz(c,a.a,i,ctx);fmpz_mpoly_get_term_exp_ui(e,a.a,i,ctx);char*p=fmpz_get_str(nullptr,10,c);ostringstream o;o<<p;flint_free(p);for(auto x:e)o<<' '<<x;o<<'\n';fmpz_clear(c);return o.str();}
string digest(const P&a){EVP_MD_CTX*d=EVP_MD_CTX_new();assert(d&&EVP_DigestInit_ex(d,EVP_sha256(),nullptr));for(long i=0;i<terms(a);i++){string l=termline(a,i);assert(EVP_DigestUpdate(d,l.data(),l.size()));}unsigned char h[32];unsigned n;assert(EVP_DigestFinal_ex(d,h,&n)&&n==32);EVP_MD_CTX_free(d);const char*x="0123456789abcdef";string s;for(auto c:h){s+=x[c>>4];s+=x[c&15];}return s;}
using U=vector<P>;
void trim(U&a){while(!a.empty()&&zero(a.back()))a.pop_back();}
int deg(const U&a){return int(a.size())-1;}
long allterms(const U&a){long n=0;for(auto&p:a)n+=terms(p);return n;}
struct Division {U remainder;vector<int> degree_history;unsigned exponent;};
Division prem(const U&A,const U&B){
 assert(deg(A)>=deg(B)&&!B.empty());U r=A;P l=B.back();unsigned exponent=deg(A)-deg(B)+1;vector<int> hist;
 for(unsigned i=0;i<exponent;i++){
  int before=deg(r);hist.push_back(before);
  if(before<deg(B)){for(auto&c:r)c=c*l;trim(r);continue;}
  int shift=before-deg(B);P c=r.back();for(auto&x:r)x=x*l;
  for(unsigned j=0;j<B.size();j++)r[shift+j]=r[shift+j]-c*B[j];
  trim(r);assert(deg(r)<before);
 }
 assert(deg(r)<deg(B));return{move(r),hist,exponent};
}
U divexact(const U&a,const P&den){U b(a.size());for(unsigned i=0;i<a.size();i++){assert(fmpz_mpoly_divides(b[i].a,a[i].a,den.a,ctx));assert(equal(b[i]*den,a[i]));}trim(b);return b;}
string output_dir;uint64_t bytes_used=0;const uint64_t cap=31*1024*1024;
void savepoly(ostream&o,const P&p){o<<terms(p)<<'\n';for(long i=0;i<terms(p);i++)o<<termline(p,i);}
void saveu(const string&name,const U&a){
 uint64_t size=128;for(auto&p:a)for(long i=0;i<terms(p);i++)size+=termline(p,i).size();
 if(size+bytes_used>cap)throw runtime_error("32MiB evidence cap before "+name);
 ofstream o(output_dir+"/"+name+".txt");o<<"5 "<<a.size()<<"\n";for(auto&p:a)savepoly(o,p);o.close();bytes_used+=size;
}
void front(const string&stage,long n=-1){ofstream f(output_dir+"/native-frontier.json");f<<"{\"stage\":\""<<stage<<"\",\"latest_total_coefficient_terms\":"<<n<<"}\n";cerr<<stage<<" TERMS "<<n<<endl;}
void division_record(const string&name,const Division&d,const string&description){ofstream o(output_dir+"/"+name+"-recipe.json");o<<"{\"pseudo_division_exponent\":"<<d.exponent<<",\"successive_remainder_degrees\":[";for(unsigned i=0;i<d.degree_history.size();i++){if(i)o<<',';o<<d.degree_history[i];}o<<"],\"quotient_DAG\":\"Q_next=LC(B)*Q+LC(R)*s^(degR-degB); use zero addition if degR<degB\",\"remainder_DAG\":\"R_next=LC(B)*R-LC(R)*s^(degR-degB)*B, repeated exponent times; padded by LC(B) if degree already low\",\"literal_identity\":\"LC(B)^exponent*A=Q_DAG*B+remainder\",\"normalization\":\""<<description<<"\"}\n";}
U constants(initializer_list<long> v){U a;for(auto x:v)a.emplace_back(x);trim(a);return a;}
void fixture(const U&A,const U&B,const string&name){
 assert(deg(A)==5&&deg(B)==3);auto d2=prem(A,B);U c2=d2.remainder;for(auto&x:c2)x=neg(x);assert(deg(c2)==2);
 auto d1=prem(B,c2);U c1=divexact(d1.remainder,pw(B.back(),3));assert(deg(c1)==1);
 auto d0=prem(c2,c1);U c0=divexact(d0.remainder,pw(c2.back(),2));assert(deg(c0)==0);
 fmpz_mat_t m;fmpz_mat_init(m,8,8);fmpz_t z;fmpz_init(z);ulong e[5]={0,0,0,0,0};
 for(int i=0;i<5;i++)for(int j=0;j<=3;j++){fmpz_mpoly_get_coeff_fmpz_ui(z,B[3-j].a,e,ctx);fmpz_set(fmpz_mat_entry(m,i,i+j),z);}
 for(int i=0;i<3;i++)for(int j=0;j<=5;j++){fmpz_mpoly_get_coeff_fmpz_ui(z,A[5-j].a,e,ctx);fmpz_set(fmpz_mat_entry(m,5+i,i+j),z);}
 fmpz_t det,got;fmpz_init(det);fmpz_init(got);fmpz_mat_det(det,m);fmpz_mpoly_get_coeff_fmpz_ui(got,c0[0].a,e,ctx);fmpz_neg(got,got);assert(fmpz_equal(det,got)&&!fmpz_is_zero(det));
 char*t=fmpz_get_str(nullptr,10,det);ofstream out(output_dir+"/fixture-"+name+".json");out<<"{\"passed\":true,\"independent8x8_Bareiss_determinant\":\""<<t<<"\",\"PRS_degrees\":[5,3,2,1,0],\"author_input_order_resultant_equals_negative_C0\":true}\n";flint_free(t);fmpz_clear(det);fmpz_clear(got);fmpz_clear(z);fmpz_mat_clear(m);
}
U load_input(istream&f){long n;f>>n;assert(n>0);U a;fmpz_t c;fmpz_init(c);ulong e[6],co[5];
 for(long i=0;i<n;i++){string t;f>>t;for(auto&x:e)f>>x;assert(f&&fmpz_set_str(c,t.c_str(),10)==0);while(a.size()<=e[0])a.emplace_back();for(int j=0;j<5;j++)co[j]=e[5-j];fmpz_mpoly_push_term_fmpz_ui(a[e[0]].a,c,co,ctx);}
 for(auto&p:a){fmpz_mpoly_sort_terms(p.a,ctx);fmpz_mpoly_combine_like_terms(p.a,ctx);}fmpz_clear(c);trim(a);return a;
}
P load_result(const string&path){ifstream f(path);int nv;long n;f>>nv>>n;assert(nv==6&&n==47114);P p;fmpz_t c;fmpz_init(c);ulong e[6],co[5];for(long i=0;i<n;i++){string t;f>>t;for(auto&x:e)f>>x;assert(f&&e[0]==0&&fmpz_set_str(c,t.c_str(),10)==0);for(int j=0;j<5;j++)co[j]=e[5-j];fmpz_mpoly_push_term_fmpz_ui(p.a,c,co,ctx);}string extra;assert(!(f>>extra));fmpz_mpoly_sort_terms(p.a,ctx);fmpz_mpoly_combine_like_terms(p.a,ctx);fmpz_clear(c);assert(terms(p)==47114);return p;}
int main(int argc,char**argv){assert(argc==4);output_dir=argv[3];fmpz_mpoly_ctx_init(ctx,5,ORD_LEX);int status=0;auto start=chrono::steady_clock::now();{
 try{
  fixture(constants({1,3,-1,0,2,1}),constants({-1,2,1,1}),"monic");
  fixture(constants({1,-5,2,-1,3,2}),constants({2,4,-2,3}),"nonmonic");front("tiny_monic_and_nonmonic_fixture_passed");
  ifstream in(argv[1]);int nv,elim;in>>nv>>elim;assert(nv==6&&elim==0);U B=load_input(in),A=load_input(in);string extra;assert(!(in>>extra));assert(deg(A)==5&&deg(B)==3);saveu("input-A5",A);saveu("input-B3",B);front("starting_C2_pseudo_remainder");
  auto d2=prem(A,B);U c2=move(d2.remainder);for(auto&x:c2)x=neg(x);assert(deg(c2)==2);saveu("C2",c2);division_record("C2",d2,"C2=-prem(A5,B3)");front("C2_saved",allterms(c2));
  auto d1=prem(B,c2);front("C1_pseudo_remainder_ready",allterms(d1.remainder));P den1=pw(B.back(),3);U c1=divexact(d1.remainder,den1);assert(deg(c1)==1);saveu("C1",c1);division_record("C1",d1,"C1=prem(B3,C2)/LC(B3)^3; exact coefficient divisions all checked by remultiplication");saveu("denominator-C1",U{den1});front("C1_saved",allterms(c1));
  auto d0=prem(c2,c1);front("C0_pseudo_remainder_ready",allterms(d0.remainder));P den0=pw(c2.back(),2);U c0=divexact(d0.remainder,den0);assert(deg(c0)==0);division_record("C0",d0,"C0=prem(C2,C1)/LC(C2)^2; exact coefficient division checked by remultiplication");saveu("C0",c0);saveu("denominator-C0",U{den0});front("C0_saved",allterms(c0));
  P expected=load_result(argv[2]);assert(equal(neg(c0[0]),expected));front("all47114_raw_resultant_coefficients_match_negative_C0",terms(expected));
  ofstream report(output_dir+"/native-report.json");report<<"{\"all_passed\":true,\"PRS_degrees\":[5,3,2,1,0],\"coefficient_variable_order\":[\"Z2\",\"Z1\",\"Z0\",\"R\",\"r\"],\"coefficient_term_order\":\"LEX\",\"raw_all47114_integer_coefficients_match\":true,\"C2_terms\":"<<allterms(c2)<<",\"C1_terms\":"<<allterms(c1)<<",\"C0_terms\":"<<allterms(c0)<<",\"C0_reversedLEX_sha256\":\""<<digest(c0[0])<<"\",\"no_resultant_gcd_point_prime_or_source26_call\":true,\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
 }catch(const exception&e){ofstream o(output_dir+"/native-failure.json");o<<"{\"all_passed\":false,\"error\":\""<<e.what()<<"\"}\n";cerr<<e.what()<<endl;status=2;}
 }fmpz_mpoly_ctx_clear(ctx);flint_cleanup();return status;}
