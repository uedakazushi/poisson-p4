// Q1 ONLY, fixed known-unit J Hasse valuation. --help never evaluates anything.
#include <flint/fmpz.h>
#include <sys/resource.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
struct Z{fmpz_t v;Z(){fmpz_init(v);}Z(Z&&o)noexcept{fmpz_init(v);fmpz_swap(v,o.v);}Z(const Z&o){fmpz_init(v);fmpz_set(v,o.v);}~Z(){fmpz_clear(v);}};
struct Rec{long long key;Z c;};
struct Stop:std::runtime_error{using std::runtime_error::runtime_error;};
using Poly=std::map<int,Z>;
static constexpr size_t OWN=29336,META=16000;
static constexpr long long MA=1048576,ME=1024;
static size_t ninput=0,powslots=0,accsize=0,peak=0;
static unsigned long long products=0,scalings=0,updates=0;
static int phase=0,kcur=-1,icursor=-1,rcursor=-1,ph=-1,pr=-1;
static std::chrono::steady_clock::time_point started;
static std::string zz(const Z&z){char*s=fmpz_get_str(nullptr,10,z.v);std::string r(s);flint_free(s);return r;}
static void own(size_t ni,size_t np,size_t na){size_t n=2*(ni+np+na)+META;peak=std::max(peak,n);if(n>OWN)throw Stop("owned29336_guard");}
static void bits(const Z&z){if(fmpz_bits(z.v)>518)throw Stop("518bit_guard");}
static void clockcheck(){if(std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()>45)throw Stop("45s_soft_guard");}
static void update(Poly&a,int key,const Z&v){
 int degreea=key/1024,e=key%1024;if(degreea<0||degreea>34||e<0||e>51-kcur||degreea+e>68-kcur)throw Stop("bivariate_support_guard");
 auto it=a.find(key);Z x;if(it==a.end())fmpz_set(x.v,v.v);else fmpz_add(x.v,it->second.v,v.v);bits(x);
 if(fmpz_is_zero(x.v)){if(it!=a.end()){a.erase(it);--accsize;}return;}
 if(it==a.end()){own(ninput,powslots,accsize+1);a.emplace(key,std::move(x));++accsize;}else fmpz_swap(it->second.v,x.v);
 ++updates;
}
static void save_prefix(const std::string&p,const std::string&why,const std::vector<Rec>&in,const std::vector<std::vector<Z>>&P,const Poly&a,const Z&D,const Z&input_l1){
 std::ofstream o(p);o<<"Q1JHP001\n"<<phase<<' '<<kcur<<' '<<icursor<<' '<<rcursor<<' '<<ph<<' '<<pr<<' '<<products<<' '<<scalings<<' '<<updates<<' '<<peak<<'\n'<<why<<'\n'<<zz(D)<<' '<<zz(input_l1)<<'\n';
 o<<in.size()<<'\n';for(auto&r:in)o<<r.key<<' '<<zz(r.c)<<'\n';
 o<<P.size()<<'\n';for(auto&v:P){o<<v.size()<<'\n';for(auto&c:v)o<<zz(c)<<'\n';}
 o<<a.size()<<'\n';for(auto&[key,c]:a)o<<key<<' '<<zz(c)<<'\n';
 if(!o)throw std::runtime_error("prefix_device_failure");
}
int main(int argc,char**argv){
 if(argc==2&&std::string(argv[1])=="--help"){
  std::cout<<"--run INPUT OUTDIR BASE_RAW. ONE Q1/J Hasse valuation. d-root=(6e+4)/(5a), externalden=a^18; 0..18 firstnonzero, no quotient/gcd; own29336,bits518,45soft60hard,2GB,RAW8000000.\n";return 0;
 }
 if(argc!=5||std::string(argv[1])!="--run")return 2;
 if(std::stoull(argv[4])+4000000+200000>8000000){std::cerr<<"RAW closure fails before arithmetic\n";return 3;}
 rlimit ram{2000000000,2000000000},cpu{60,60};setrlimit(RLIMIT_AS,&ram);setrlimit(RLIMIT_CPU,&cpu);started=std::chrono::steady_clock::now();
 std::ifstream f(argv[2]);std::string magic,s;size_t n;int dmax,shift;f>>magic>>n>>dmax>>shift;if(magic!="Q1JH001"||n!=4811||dmax!=18||shift!=0)return 4;
 Z D,input_l1,abs_c;f>>s;fmpz_set_str(D.v,s.c_str(),10);std::vector<Rec>in;in.reserve(n);std::vector<std::vector<Z>>P;Poly acc;
 try{
  // The copied input is the literal ROOT-normalized Q1 integer/common-D representation.
  for(size_t i=0;i<n;++i){icursor=(int)i;Rec r;if(!(f>>r.key>>s))throw std::runtime_error("literal_input_IO_failure");fmpz_set_str(r.c.v,s.c_str(),10);bits(r.c);
   Z nextL;fmpz_abs(abs_c.v,r.c.v);fmpz_add(nextL.v,input_l1.v,abs_c.v);if(fmpz_bits(nextL.v)>301)throw Stop("input_l1_301bit_guard");
   long long a=r.key/MA,rem=r.key%MA;int e=(int)(rem/ME),d=(int)(rem%ME);
   if(a<0||a>16||e<0||e>33||d<0||d>18||a+e+d>50)throw Stop("literal_Q1_support");
   own(ninput+1,powslots,accsize);in.push_back(std::move(r));++ninput;fmpz_swap(input_l1.v,nextL.v);
  }
  if(!f)throw Stop("literal_input_truncated");
  // Fixed metadata is counted in META16000. No powers of source Q1.
  unsigned long long choose[19][19]{};for(int i=0;i<=18;++i){choose[i][0]=choose[i][i]=1;for(int j=1;j<i;++j)choose[i][j]=choose[i-1][j-1]+choose[i-1][j];}
  std::vector<Z>scale5(19);fmpz_one(scale5[0].v);for(int i=1;i<=18;++i){fmpz_mul_ui(scale5[i].v,scale5[i-1].v,5);++scalings;}
  phase=1;P.resize(19);for(int h=0;h<=18;++h){own(ninput,powslots+h+1,accsize);P[h].resize(h+1);powslots+=h+1;}
  fmpz_one(P[0][0].v);int root[2]={4,6};
  for(int h=1;h<=18;++h)for(int r=0;r<(int)P[h-1].size();++r)for(int j=0;j<2;++j){
   ph=h;pr=2*r+j;Z t,next;fmpz_mul_si(t.v,P[h-1][r].v,root[j]);++scalings;fmpz_add(next.v,P[h][r+j].v,t.v);bits(next);fmpz_swap(next.v,P[h][r+j].v);
  }
  std::ofstream zeros(std::string(argv[3])+"/completed-zero-orders.txt");
  auto convolution=[&](int k,int pass){
   phase=pass;Z x,y,z;int begin=pass==3?(int)in.size()-1:0,end=pass==3?-1:(int)in.size(),step=pass==3?-1:1;
   for(int ii=begin;ii!=end;ii+=step){
    icursor=ii;long long aa=in[ii].key/MA,rem=in[ii].key%MA;int e=(int)(rem/ME),d=(int)(rem%ME),i=(int)aa;
    if(d<k)continue;int h=d-k;
    int rb=pass==3?(int)P[h].size()-1:0,re=pass==3?-1:(int)P[h].size(),rs=pass==3?-1:1;rcursor=rb;
    fmpz_mul_ui(x.v,scale5[18-h].v,choose[d][k]);fmpz_mul(y.v,in[ii].c.v,x.v);scalings+=2;bits(y);
    for(int r=rb;r!=re;r+=rs){
     rcursor=r;if((products&4095)==0)clockcheck();fmpz_mul(z.v,y.v,P[h][r].v);if(pass==3)fmpz_neg(z.v,z.v);bits(z);update(acc,(i+18-h)*1024+e+r,z);++products;
    }
   }
  };
  for(int k=0;k<=18;++k){
   kcur=k;phase=2;icursor=rcursor=0;convolution(k,2);
   if(acc.empty()){zeros<<k<<" ZERO\n";zeros.flush();continue;}
   // Verify only the first nonzero layer, in-place, then reproduce it.
   convolution(k,3);if(!acc.empty())throw Stop("reverse_Hasse_nonzero");
   convolution(k,4);if(acc.empty())throw Stop("reproduction_unexpected_zero");
   phase=5;icursor=rcursor=-1;clockcheck();
   Z Dclear;fmpz_mul(Dclear.v,D.v,scale5[18].v);++scalings;bits(Dclear);
   std::ofstream o(std::string(argv[3])+"/first-nonzero-J-Hasse.txt");o<<"Q1JHV001\n"<<k<<' '<<zz(Dclear)<<' '<<18<<'\n'<<acc.size()<<'\n';for(auto&[key,c]:acc)o<<key/1024<<' '<<key%1024<<' '<<zz(c)<<'\n';
   if(!o)throw std::runtime_error("output_device_failure_no_second_full_prefix");
   std::ofstream r(std::string(argv[3])+"/native-report.json");r<<"{\"all_passed\":true,\"J_valuation\":"<<k<<",\"first_nonzero_terms\":"<<acc.size()<<",\"peak_owned\":"<<peak<<",\"integer_products\":"<<products<<",\"scalar_root_products\":"<<scalings<<",\"updates\":"<<updates<<",\"reverse_first_layer_zero\":true,\"reproduced_nonzero\":true}\n";return 0;
  }
  throw Stop("all19_layers_zero_contradicts_nonzero_input_no_retry");
 }catch(const Stop&e){
  save_prefix(std::string(argv[3])+"/Q1JHP001-prefix.txt",e.what(),in,P,acc,D,input_l1);
  std::ofstream r(std::string(argv[3])+"/native-report.json");r<<"{\"all_passed\":false,\"prefix_only\":true,\"reason\":\""<<e.what()<<"\",\"peak_owned\":"<<peak<<",\"integer_products\":"<<products<<"}\n";return 10;
 }
}
