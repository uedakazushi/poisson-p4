// ONE fixed-prime Kronecker gcd, no multivariate factorization or parameter sample.
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <stdexcept>
#include <chrono>
#include <sys/resource.h>
using P=std::vector<unsigned>;const unsigned MOD=32003;unsigned long long products=0;int phase=0,kcursor=-1,pass=0;auto start=std::chrono::steady_clock::now();
unsigned mul(unsigned a,unsigned b){++products;if(products>150000000)throw std::runtime_error("product_guard");return (a*b)%MOD;}
unsigned inv(unsigned a){if(!a)throw std::runtime_error("zero_inverse");unsigned v=1;for(unsigned n=MOD-2;n;n>>=1,a=mul(a,a))if(n&1)v=mul(v,a);return v;}
void trim(P&f){while(!f.empty()&&!f.back())f.pop_back();}
void clockcheck(){if(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()>45)throw std::runtime_error("45s_soft_guard");}
P read(const char*path,const std::string&magic,int n,int degree){
 std::ifstream s(path);std::string m,D,c;int n0;s>>m>>D>>n0;if(m!=magic||n!=n0||D!="164144011096127878172550148707469312481250000000000000000000")throw std::runtime_error("input_header");P f(degree+1);long long old=-1;
 for(int i=0;i<n;++i){long long key;s>>key>>c;if(!s||key<=old||c.empty())throw std::runtime_error("input_body");old=key;int a=key/1048576,e=(key%1048576)/1024,d=key%1024;if(a<0||a>16||e<0||e>33||d<0||d>18)throw std::runtime_error("input_box");int w=646*a+19*e+d;if(w>degree||f[w])throw std::runtime_error("Kronecker_collision_or_degree");unsigned v=0;bool neg=c[0]=='-';for(size_t j=neg;j<c.size();++j){if(c[j]<'0'||c[j]>'9')throw std::runtime_error("integer_digit");v=(v*10+c[j]-'0')%MOD;}f[w]=neg&&v?MOD-v:v;}
 if(s>>c)throw std::runtime_error("extra_input");if(f.empty()||!f.back())throw std::runtime_error("bad_prime_leading_degree_drop");return f;
}
int main(int argc,char**argv){if(argc==2&&std::string(argv[1])=="--help"){std::cout<<"--run N0 N2 OUTDIR; prime32003, weights646/19/1; 25000 owned,31bit scalar,150M products,45soft60hard,2GB,RAW4M; no auto retry.\n";return 0;}if(argc!=5||std::string(argv[1])!="--run")return 2;rlimit ram{2000000000,2000000000},cpu{60,60};setrlimit(RLIMIT_AS,&ram);setrlimit(RLIMIT_CPU,&cpu);start=std::chrono::steady_clock::now();P f,g;
 try{f=read(argv[2],"T0C4QV001",4500,10544);g=read(argv[3],"Q1B4QV001",2317,7975);phase=1;if(f.capacity()+g.capacity()+4096>25000)throw std::runtime_error("owned_guard");
 while(!g.empty()){
  unsigned il=inv(g.back());
  while(f.size()>=g.size()&&!f.empty()){
   kcursor=(int)(f.size()-g.size());unsigned c=mul(f.back(),il);
   for(size_t j=0;j<g.size();++j){unsigned v=mul(c,g[j]);f[kcursor+j]=(f[kcursor+j]+MOD-v)%MOD;}
   trim(f);if((products&8191)<g.size())clockcheck();
  }
  f.swap(g);++pass;clockcheck();
 }
 phase=2;unsigned il=inv(f.back());for(auto&x:f)x=mul(x,il);std::ofstream o(std::string(argv[4])+"/monic-gcd-mod32003.txt");o<<f.size()-1<<'\n';for(size_t i=0;i<f.size();++i)if(f[i])o<<i<<' '<<f[i]<<'\n';
 std::ofstream r(std::string(argv[4])+"/native-report.json");r<<"{\"all_passed\":true,\"prime\":32003,\"input_degrees_preserved\":[10544,7975],\"gcd_degree\":"<<f.size()-1<<",\"gcd_is_one\":"<<(f.size()==1?"true":"false")<<",\"euclidean_passes\":"<<pass<<",\"products\":"<<products<<",\"owned_upper\":25000,\"no_char0_factor_claim_if_nonconstant\":true}\n";return 0;
 }catch(const std::exception&e){std::ofstream o(std::string(argv[4])+"/prefix.txt");o<<"KRONPREFIX1\n"<<phase<<' '<<pass<<' '<<kcursor<<' '<<products<<'\n';for(P*p:{&f,&g}){o<<p->size()<<'\n';for(auto x:*p)o<<x<<' ';o<<'\n';}std::ofstream r(std::string(argv[4])+"/native-report.json");r<<"{\"all_passed\":false,\"prefix_saved\":true,\"reason\":\""<<e.what()<<"\"}\n";return 10;}
}
