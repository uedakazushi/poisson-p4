// ONE Q1 named-B case. Preparation --help never builds roots or Hasse values.
#include <flint/fmpz.h>
#include <sys/resource.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
struct Z{fmpz_t v;Z(){fmpz_init(v);}Z(const Z&o){fmpz_init(v);fmpz_set(v,o.v);}Z(Z&&o)noexcept{fmpz_init(v);fmpz_swap(v,o.v);}~Z(){fmpz_clear(v);}};
struct Rec{long long key;Z c;};using Poly=std::map<int,Z>;
struct Stop:std::runtime_error{using std::runtime_error::runtime_error;};
static constexpr size_t OWN=37436,META=16000;static constexpr int EBASE=128;
static constexpr long long MA=1048576,ME=1024;
static size_t ni=0,nr=0,np=0,nn=0,ns=0,na=0,peak=0;
static unsigned long long products=0,scalings=0,updates=0;
static int phase=0,tc=-1,hc=-1,rc=-1,jc=-1,bc=-1,ic=-1,kc=-1,pk=0,sk=0,base_ready=0;
static std::chrono::steady_clock::time_point started;
static std::string str(const Z&z){char*s=fmpz_get_str(nullptr,10,z.v);std::string v(s);flint_free(s);return v;}
static void bits(const Z&z){if(fmpz_bits(z.v)>518)throw Stop("518bit_before_commit");}
static void own(size_t extra=0){size_t z=2*(ni+nr+np+nn+ns+na+extra)+META;if(z>OWN)throw Stop("37436_owned_before_commit");peak=std::max(peak,z);}
static void clockcheck(){if(std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()>45)throw Stop("45s_soft_guard");}
static void update(Poly&p,int key,const Z&v,size_t&n,int kind){
 int a=key/EBASE,e=key%EBASE;bool ok=a>=0&&e>=0;
 if(kind==0)ok=ok&&a<=18-hc&&e<=36-hc-2*a;
 if(kind==1)ok=ok&&a<=tc+1&&e<=2*(tc+1-a);
 if(kind==2)ok=ok&&a<=34&&e<=69&&a+e<=86-kc;
 if(!ok)throw Stop("exact_source_support_guard");
 auto it=p.find(key);Z z;if(it==p.end())fmpz_set(z.v,v.v);else fmpz_add(z.v,it->second.v,v.v);bits(z);
 if(fmpz_is_zero(z.v)){if(it!=p.end()){p.erase(it);--n;}return;}
 if(it==p.end()){own(1);p.emplace(key,std::move(z));++n;}else fmpz_swap(it->second.v,z.v);++updates;
}
static void polyout(std::ostream&o,const Poly&p){o<<p.size()<<'\n';for(auto&[key,c]:p)o<<key<<' '<<str(c)<<'\n';}
static void prefix(const std::string&path,const std::string&why,const Z&D,const Z&L,const Z&base,const std::vector<Rec>&in,const std::vector<std::vector<Z>>&R,const Poly&P,const Poly&N,const std::vector<Poly>&S,const Poly&A){
 std::ofstream o(path);o<<"Q1BHP001\n"<<phase<<' '<<tc<<' '<<hc<<' '<<rc<<' '<<jc<<' '<<bc<<' '<<ic<<' '<<kc<<' '<<pk<<' '<<sk<<' '<<base_ready<<' '<<products<<' '<<scalings<<' '<<updates<<' '<<peak<<'\n'<<why<<'\n'<<str(D)<<' '<<str(L)<<' '<<str(base)<<'\n';
 o<<in.size()<<'\n';for(auto&r:in)o<<r.key<<' '<<str(r.c)<<'\n';o<<R.size()<<'\n';for(auto&r:R){o<<r.size()<<'\n';for(auto&c:r)o<<str(c)<<'\n';}
 polyout(o,P);polyout(o,N);o<<S.size()<<'\n';for(auto&p:S)polyout(o,p);polyout(o,A);
 if(!o)throw std::runtime_error("prefix_IO_failure_no_duplicate_state");
}
int main(int argc,char**argv){
 if(argc==2&&std::string(argv[1])=="--help"){
  std::cout<<"--run INPUT OUTDIR BASE_RAW. ONE canonicalQ1/B valuation 0..18 firstnonzero. Cache S_h=R^h*b1^(18-h); denominator D*b1^18 SYMBOLIC, b1 NOT assumed FULLUunit. owned37436/bits518/45soft60hard/2GB/RAW8M.\n";return 0;
 }
 if(argc!=5||std::string(argv[1])!="--run")return 2;
 if(std::stoull(argv[4])+5000000+200000>8000000){std::cerr<<"RAW closure fails BEFORE arithmetic\n";return 3;}
 rlimit ram{2000000000,2000000000},cpu{60,60};setrlimit(RLIMIT_AS,&ram);setrlimit(RLIMIT_CPU,&cpu);started=std::chrono::steady_clock::now();
 std::ifstream f(argv[2]);std::string magic,cstr;size_t n;int dm,shift;f>>magic>>n>>dm>>shift;
 if(magic!="Q1BH001"||n!=4811||dm!=18||shift!=0)return 4;
 Z D,L,activebase;f>>cstr;fmpz_set_str(D.v,cstr.c_str(),10);
 std::vector<Rec>in;in.reserve(n);std::vector<std::vector<Z>>R;Poly P,N,A;std::vector<Poly>S;
 try{
  for(size_t i=0;i<n;++i){ic=(int)i;Rec q;if(!(f>>q.key>>cstr))throw std::runtime_error("literal_input_IO_failure");fmpz_set_str(q.c.v,cstr.c_str(),10);bits(q.c);
   int a=(int)(q.key/MA),e=(int)((q.key%MA)/ME),d=(int)(q.key%ME);if(a<0||a>16||e<0||e>33||d<0||d>18||a+e+d>50||fmpz_is_zero(q.c.v))throw Stop("literal_Q1_support");
   Z ac,next;fmpz_abs(ac.v,q.c.v);fmpz_add(next.v,L.v,ac.v);if(fmpz_bits(next.v)>301)throw Stop("input_L1_301bit_before_commit");own(1);in.push_back(std::move(q));++ni;fmpz_swap(L.v,next.v);
  }
  unsigned long long choose[19][19]{};for(int l=0;l<=18;++l){choose[l][0]=choose[l][l]=1;for(int k=1;k<l;++k)choose[l][k]=choose[l-1][k-1]+choose[l-1][k];}
  phase=1;R.resize(19);for(int h=0;h<=18;++h){own(h+1);R[h].resize(h+1);nr+=h+1;}fmpz_one(R[0][0].v);int rr[2]={40,-60};
  for(int h=1;h<=18;++h)for(int r=0;r<h;++r)for(int j=0;j<2;++j){hc=h;rc=r;jc=j;Z v,z;fmpz_mul_si(v.v,R[h-1][r].v,rr[j]);bits(v);fmpz_add(z.v,R[h][r+j].v,v.v);bits(z);fmpz_swap(z.v,R[h][r+j].v);++scalings;}
  S.resize(19);Z one;fmpz_one(one.v);own(1);P.emplace(0,std::move(one));np=1;int bkey[4]={128,2,1,0},bcoef[4]={50,18,-9,-2};
  for(int t=0;t<=18;++t){tc=t;hc=18-t;phase=2;
   for(auto&[pkey,pv]:P)for(int r=0;r<=hc;++r){pk=pkey;rc=r;if(((products+scalings)&4095)==0)clockcheck();Z v;fmpz_mul(v.v,pv.v,R[hc][r].v);bits(v);update(S[hc],pkey+r,v,ns,0);++products;}
   if(t==18)break;phase=3;
   for(auto&[pkey,pv]:P)for(int j=0;j<4;++j){pk=pkey;bc=j;if(((products+scalings)&4095)==0)clockcheck();Z v;fmpz_mul_si(v.v,pv.v,bcoef[j]);bits(v);update(N,pkey+bkey[j],v,nn,1);++scalings;}
   P.clear();np=0;P.swap(N);std::swap(np,nn);
  }
  // Cache construction is DONE. Release ALL R and current/next b1 powers before A.
  phase=4;R.clear();R.shrink_to_fit();nr=0;P.clear();np=0;N.clear();nn=0;own();
  std::ofstream zeros(std::string(argv[3])+"/completed-zero-orders.txt");
  auto convolution=[&](int k,int pass){
   phase=pass;int begin=pass==6?(int)in.size()-1:0,end=pass==6?-1:(int)in.size(),step=pass==6?-1:1;
   for(int i=begin;i!=end;i+=step){ic=i;base_ready=0;int d=(int)(in[i].key%ME);if(d<k)continue;hc=d-k;
    Z z;fmpz_mul_ui(z.v,in[i].c.v,choose[d][k]);bits(z);fmpz_swap(activebase.v,z.v);++scalings;base_ready=1;
    int ia=(int)(in[i].key/MA),ie=(int)((in[i].key%MA)/ME);int offset=ia*EBASE+ie;
    auto term=[&](int key,const Z&v){sk=key;if(((products+scalings)&4095)==0)clockcheck();Z out;fmpz_mul(out.v,activebase.v,v.v);if(pass==6)fmpz_neg(out.v,out.v);bits(out);update(A,offset+key,out,na,2);++products;};
    if(pass==6){for(auto it=S[hc].rbegin();it!=S[hc].rend();++it)term(it->first,it->second);}else for(auto&[key,v]:S[hc])term(key,v);
   }base_ready=0;
  };
  for(int k=0;k<=18;++k){kc=k;convolution(k,5);
   if(A.empty()){zeros<<k<<" ZERO\n";zeros.flush();continue;}
   convolution(k,6);if(!A.empty())throw Stop("reverse_first_layer_NONZERO");convolution(k,7);if(A.empty())throw Stop("reproduction_unexpected_ZERO");
   phase=8;ic=-1;base_ready=0;clockcheck();std::ofstream o(std::string(argv[3])+"/first-nonzero-B-Hasse.txt");
   o<<"Q1BHV001\n"<<k<<' '<<str(D)<<' '<<18<<'\n'<<A.size()<<'\n';for(auto&[key,v]:A)o<<key/EBASE<<' '<<key%EBASE<<' '<<str(v)<<'\n';
   if(!o)throw std::runtime_error("output_IO_failure_no_second_full_prefix");
   std::ofstream report(std::string(argv[3])+"/native-report.json");report<<"{\"all_passed\":true,\"B_valuation\":"<<k<<",\"first_nonzero_terms\":"<<A.size()<<",\"peak_owned\":"<<peak<<",\"root_and_inner_products\":"<<products<<",\"scalar_products\":"<<scalings<<",\"updates\":"<<updates<<",\"cache_terms\":"<<ns<<",\"reverse_first_layer_zero\":true,\"reproduced_nonzero\":true,\"symbolic_b1_power18\":true,\"b1_NOT_assumed_FULLU_unit\":true}\n";return 0;
  }throw Stop("all19_zero_contradicts_nonzero_input_no_retry");
 }catch(const Stop&e){
  prefix(std::string(argv[3])+"/Q1BHP001-prefix.txt",e.what(),D,L,activebase,in,R,P,N,S,A);
  std::ofstream report(std::string(argv[3])+"/native-report.json");report<<"{\"all_passed\":false,\"prefix_only\":true,\"reason\":\""<<e.what()<<"\",\"peak_owned\":"<<peak<<",\"root_and_inner_products\":"<<products<<",\"scalar_products\":"<<scalings<<"}\n";return 10;
 }
}
