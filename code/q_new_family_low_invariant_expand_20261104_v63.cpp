// One bounded integer FLINT expansion, no beta point and no Inv24 evaluation.
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
using namespace std;
fmpz_mpoly_ctx_t ctx;
const long TERM_CAP=150000,OUTPUT_TERMS=100000;const long BYTE_CAP=32L*1024*1024;
long max_terms=0,multiplications=0,output_bytes=0;
struct P{
 fmpz_mpoly_t v;P(){fmpz_mpoly_init(v,ctx);}P(const P&p):P(){fmpz_mpoly_set(v,p.v,ctx);}
 P(P&&p):P(){fmpz_mpoly_swap(v,p.v,ctx);}P&operator=(P p){fmpz_mpoly_swap(v,p.v,ctx);return *this;}
 ~P(){fmpz_mpoly_clear(v,ctx);}
};
void guard(const P&p){long n=fmpz_mpoly_length(p.v,ctx);max_terms=max(max_terms,n);if(n>TERM_CAP)throw runtime_error("intermediate term cap");}
P add(const P&a,const P&b){P c;fmpz_mpoly_add(c.v,a.v,b.v,ctx);guard(c);return c;}
P scale(long x,const P&a){P c;fmpz_mpoly_scalar_mul_si(c.v,a.v,x,ctx);guard(c);return c;}
P mul(const P&a,const P&b){P c;fmpz_mpoly_mul(c.v,a.v,b.v,ctx);multiplications++;guard(c);return c;}
using M=vector<P>;
M mm(const M&a,const M&b){M c(16);for(int i=0;i<4;i++)for(int j=0;j<4;j++)for(int k=0;k<4;k++)c[4*i+j]=add(c[4*i+j],mul(a[4*i+k],b[4*k+j]));return c;}
P tr(const M&a){P q;for(int i=0;i<4;i++)q=add(q,a[5*i]);return q;}
P trpair(const M&a,const M&b){P q;for(int i=0;i<4;i++)for(int j=0;j<4;j++)q=add(q,mul(a[4*i+j],b[4*j+i]));return q;}
void readpoly(istream&f,P&p){long n;f>>n;assert(f&&n>=0);fmpz_t z;fmpz_init(z);ulong e[4];for(long k=0;k<n;k++){string c;f>>c;for(int i=0;i<4;i++)f>>e[i];assert(f&&fmpz_set_str(z,c.c_str(),10)==0);fmpz_mpoly_push_term_fmpz_ui(p.v,z,e,ctx);}fmpz_mpoly_sort_terms(p.v,ctx);fmpz_mpoly_combine_like_terms(p.v,ctx);fmpz_clear(z);guard(p);}
void writebytes(ofstream&f,const string&s){if(output_bytes+(long)s.size()>BYTE_CAP)throw runtime_error("32MiB coefficient serialization cap");f<<s;output_bytes+=s.size();}
string zstr(const fmpz_t z){char*s=fmpz_get_str(nullptr,10,z);string a(s);flint_free(s);return a;}
string save(const string&dir,const string&name,const P&raw,long beta_degree,long scale_factor,long scalar_degree,long L){
 P p(raw);long n=fmpz_mpoly_length(p.v,ctx);if(n>OUTPUT_TERMS)throw runtime_error("100k output term cap");
 fmpz_t den,g,z;fmpz_init(den);fmpz_init(g);fmpz_init(z);fmpz_set_si(den,L);fmpz_pow_ui(den,den,beta_degree);fmpz_mul_si(den,den,scale_factor);fmpz_set(g,den);
 for(long k=0;k<n&&!fmpz_is_one(g);k++){fmpz_mpoly_get_term_coeff_fmpz(z,p.v,k,ctx);fmpz_gcd(g,g,z);}
 if(!fmpz_is_one(g)){fmpz_mpoly_scalar_divexact_fmpz(p.v,p.v,g,ctx);fmpz_divexact(den,den,g);}
 ulong e[4];long cancelB=beta_degree;
 for(long k=0;k<n;k++){fmpz_mpoly_get_term_exp_ui(e,p.v,k,ctx);cancelB=min(cancelB,(long)e[3]);}
 if(cancelB){P B,power,quot;fmpz_mpoly_gen(B.v,3,ctx);assert(fmpz_mpoly_pow_ui(power.v,B.v,cancelB,ctx));assert(fmpz_mpoly_divides(quot.v,p.v,power.v,ctx));p=quot;}
 ofstream file(dir+"/"+name+".json");assert(file);ostringstream header;
 header<<"{\"name\":\""<<name<<"\",\"variables\":[\"a=u/X\",\"Y=X^2\",\"W\",\"B\"],\"source_c_power\":"<<scalar_degree
 <<",\"denominator_integer\":\""<<zstr(den)<<"\",\"denominator_B_power\":"<<beta_degree-cancelB
 <<",\"denominator_D_power\":"<<beta_degree<<",\"D\":\"147456*Y-(3072*W-1)^2\",\"terms\":[";
 writebytes(file,header.str());long bits=0;
 for(long k=0;k<n;k++){fmpz_mpoly_get_term_coeff_fmpz(z,p.v,k,ctx);fmpz_mpoly_get_term_exp_ui(e,p.v,k,ctx);assert((e[0]+e[1])%2==0);bits=max(bits,(long)fmpz_bits(z));ostringstream term;if(k)term<<',';term<<"[\""<<zstr(z)<<"\","<<e[0]<<','<<(e[0]+e[1])/2<<','<<e[2]<<','<<e[3]<<']';writebytes(file,term.str());}
 writebytes(file,"]}\n");file.close();ostringstream meta;
 meta<<"{\"name\":\""<<name<<"\",\"terms\":"<<n<<",\"max_coefficient_bits\":"<<bits<<",\"B_factors_cancelled\":"<<cancelB<<",\"D_factors_cancelled\":0,\"denominator_integer\":\""<<zstr(den)<<"\"}";
 cerr<<"SAVED "<<name<<" TERMS "<<n<<" BITS "<<bits<<" TOTAL_BYTES "<<output_bytes<<endl;fmpz_clear(den);fmpz_clear(g);fmpz_clear(z);return meta.str();
}
int main(int argc,char**argv){assert(argc==3);fmpz_mpoly_ctx_init(ctx,4,ORD_DEGLEX);auto start=chrono::steady_clock::now();int result=0;{
 try{
  ifstream in(argv[1]);int vars,count,width;long L;in>>vars>>count>>width>>L;assert(vars==4&&count==16&&width==136&&L==1080);
  vector<P>b(16);for(auto&p:b)readpoly(in,p);vector<pair<int,int>>mon(width);for(auto&ij:mon){in>>ij.first>>ij.second;assert(0<=ij.first&&ij.first<=ij.second&&ij.second<16);}
  vector<vector<long>>coeff(32,vector<long>(width));for(auto&row:coeff)for(auto&c:row)in>>c;assert(in);string extra;assert(!(in>>extra));
  vector<P>products;for(auto ij:mon)products.push_back(mul(b[ij.first],b[ij.second]));M C(16),D(16);
  for(int row=0;row<32;row++){P sum;for(int k=0;k<width;k++)if(coeff[row][k])sum=add(sum,scale(coeff[row][k],products[k]));(row<16?C[row]:D[row-16])=sum;}
  products.clear();b.clear();M CC=mm(C,C),CD=mm(C,D),DD=mm(D,D);P qC=tr(CC),qD=tr(DD),qCD=tr(CD);
  assert(fmpz_mpoly_is_zero(add(qC,scale(-2,qD)).v,ctx));assert(fmpz_mpoly_is_zero(add(qCD,qD).v,ctx));
  P F=add(add(qC,scale(-4,qCD)),scale(4,qD));vector<string>stats;
  stats.push_back(save(argv[2],"F4",F,4,16,2,L));
  vector<P>pencil(5);vector<int>perm={0,1,2,3};
  do{int sign=1;for(int i=0;i<4;i++)for(int j=i+1;j<4;j++)if(perm[i]>perm[j])sign=-sign;vector<P>cur(1);fmpz_mpoly_set_si(cur[0].v,sign,ctx);
   for(int i=0;i<4;i++){vector<P>next(cur.size()+1);for(int k=0;k<(int)cur.size();k++){next[k]=add(next[k],mul(cur[k],C[4*i+perm[i]]));next[k+1]=add(next[k+1],mul(cur[k],D[4*i+perm[i]]));}cur=move(next);}
   for(int k=0;k<5;k++)pencil[k]=add(pencil[k],cur[k]);
  }while(next_permutation(perm.begin(),perm.end()));
  assert(fmpz_mpoly_is_zero(add(add(mul(F,F),scale(-400,pencil[0])),scale(-1600,pencil[4])).v,ctx));
  for(int k=0;k<5;k++)stats.push_back(save(argv[2],"p"+to_string(k),pencil[k],8,256,4,L));pencil.clear();
  P T4=trpair(CC,DD),D4=trpair(DD,DD);M E=mm(C,DD);P T6=trpair(E,E);E.clear();
  P N=add(add(scale(8,T6),scale(-4,mul(qD,T4))),add(mul(mul(qD,qD),qC),scale(-2,mul(D4,qC))));
  stats.push_back(save(argv[2],"N12",N,12,32768,6,L));
  ofstream report(string(argv[2])+"/native-report.json");report<<"{\"all_passed\":true,\"outputs\":[";for(int i=0;i<(int)stats.size();i++){if(i)report<<',';report<<stats[i];}
  report<<"],\"max_intermediate_terms\":"<<max_terms<<",\"polynomial_multiplications\":"<<multiplications<<",\"coefficient_output_bytes\":"<<output_bytes
  <<",\"quartic_relations_and_F4_square_pencil_identity_verified\":true,\"integer_word_N12\":\"8tr(CD2CD2)-4tr(D2)tr(C2D2)+tr(D2)^2tr(C2)-2tr(D4)tr(C2)\",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
  cerr<<"ALL_SEVEN_EXACT_LOW_INVARIANT_POLYNOMIALS_PASS\n";
 }catch(const exception&e){cerr<<"BOUNDED_EXPANSION_STOP "<<e.what()<<endl;result=2;}
 }fmpz_mpoly_ctx_clear(ctx);flint_cleanup();return result;}
