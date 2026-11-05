// Fixed primitive B, FOUR lex50ad divisions. --help does no coefficient work.
#include <flint/fmpz.h>
#include <sys/resource.h>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
struct Z{fmpz_t v;Z(){fmpz_init(v);}Z(const Z&o){fmpz_init(v);fmpz_set(v,o.v);}Z(Z&&o)noexcept{fmpz_init(v);fmpz_swap(v,o.v);}~Z(){fmpz_clear(v);}};
using Poly=std::map<long long,Z>;struct Stop:std::runtime_error{using std::runtime_error::runtime_error;};
static constexpr size_t OWN=52826,META=16000;static constexpr long long MA=1048576,ME=1024;
static constexpr uint64_t SEED=14695981039346656037ULL;
static size_t nr=0,nq=0,nf=0,nn=0,peak=0,terms[4]{};
static unsigned long long products=0,scalars=0,updates=0;
static uint64_t input_digest=SEED,reread_digest=SEED,stage_digest[4]={SEED,SEED,SEED,SEED};
static int phase=0,pass=0,bc=-1,fh=0,source_i=-1,active_ready=0;
static long long aqkey=0,fkey=0,vqkey=0,vfkey=0;
static std::chrono::steady_clock::time_point started;
static std::string zs(const Z&z){char*p=fmpz_get_str(nullptr,10,z.v);std::string s(p);flint_free(p);return s;}
static uint64_t digest(uint64_t h,long long key,const Z&z){std::string s=std::to_string(key)+" "+zs(z)+"\n";for(unsigned char c:s){h^=c;h*=1099511628211ULL;}return h;}
static void bits(const Z&z){if(fmpz_bits(z.v)>518)throw Stop("518bit_before_commit");}
static void own(size_t extra=0){size_t v=2*(nr+nq+nf+nn+extra)+META;if(v>OWN)throw Stop("52826_owned_before_commit");peak=std::max(peak,v);}
static void clockcheck(){if(std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()>45)throw Stop("45s_soft_guard");}
static void decode(long long key,int&a,int&e,int&d){a=(int)(key/MA);long long x=key%MA;e=(int)(x/ME);d=(int)(x%ME);}
static void box(long long key,int j){int a,e,d;decode(key,a,e,d);if(a<0||a>16-j||e<0||e>33-2*j||d<0||d>18-j||a+e+d>50-3*j)throw Stop("stage_exact_support_guard");}
static void update(Poly&p,long long key,const Z&value,size_t&n,int j,bool factor=false){
 if(factor){int a,e,d;decode(key,a,e,d);if(a<0||a>fh||e<0||e>2*fh||d<0||d>fh||a+e+d>3*fh)throw Stop("fixed_Bpower_support");}else box(key,j);
 auto it=p.find(key);Z z;if(it==p.end())fmpz_set(z.v,value.v);else fmpz_add(z.v,it->second.v,value.v);bits(z);
 if(fmpz_is_zero(z.v)){if(it!=p.end()){p.erase(it);--n;}return;}
 if(it==p.end()){own(1);p.emplace(key,std::move(z));++n;}else fmpz_swap(it->second.v,z.v);++updates;
}
static void polyout(std::ostream&o,const Poly&p){o<<p.size()<<'\n';for(auto&[key,c]:p)o<<key<<' '<<zs(c)<<'\n';}
static void prefix(const std::string&path,const std::string&why,const Z&D,const Z&L,const Z&Lr,const Z&A,const Poly&R,const Poly&Q,const Poly&F,const Poly&N){
 std::ofstream o(path);o<<"Q1B4QP001\n"<<phase<<' '<<pass<<' '<<bc<<' '<<fh<<' '<<source_i<<' '<<active_ready<<' '<<aqkey<<' '<<fkey<<' '<<vqkey<<' '<<vfkey<<' '<<products<<' '<<scalars<<' '<<updates<<' '<<peak<<' '<<input_digest<<' '<<reread_digest<<'\n'<<why<<'\n';
 for(auto t:terms)o<<t<<' ';o<<'\n';for(auto h:stage_digest)o<<h<<' ';o<<'\n';o<<zs(D)<<' '<<zs(L)<<' '<<zs(Lr)<<' '<<zs(A)<<'\n';polyout(o,R);polyout(o,Q);polyout(o,F);polyout(o,N);
 if(!o)throw std::runtime_error("prefix_IO_failure_no_double_state");
}
int main(int argc,char**argv){
 if(argc==2&&std::string(argv[1])=="--help"){
  std::cout<<"--run INPUT OUTDIR BASE_RAW. Q1/B4 fixed6-term lex50ad integerdivision FOUR passes; originalD197/no b1inverse; R/Q only, original reread for ALL B4 multiplyback. owned52826/hard518/45soft60hard/2GB/RAW8M.\n";return 0;
 }
 if(argc!=5||std::string(argv[1])!="--run")return 2;
 if(std::stoull(argv[4])+6000000+200000>8000000){std::cerr<<"RAW closure fails BEFORE arithmetic\n";return 3;}
 rlimit ram{2000000000,2000000000},cpu{60,60};setrlimit(RLIMIT_AS,&ram);setrlimit(RLIMIT_CPU,&cpu);started=std::chrono::steady_clock::now();
 std::ifstream in(argv[2]);std::string magic,s;size_t n;int npass,shift;in>>magic>>n>>npass>>shift;if(magic!="Q1B4Q001"||n!=4811||npass!=4||shift!=0)return 4;
 Z D,L,Lr,active;in>>s;fmpz_set_str(D.v,s.c_str(),10);Poly R,Q,F,N;long long bkey[6]={MA+1,2*ME+1,ME+1,1,ME,0};int bcoef[6]={50,18,-9,-2,60,-40};
 try{
  long long prev=-1;
  for(size_t i=0;i<n;++i){source_i=(int)i;long long key;if(!(in>>key>>s))throw std::runtime_error("literal_input_IO_failure");Z c,ac,next;fmpz_set_str(c.v,s.c_str(),10);box(key,0);if(!fmpz_sgn(c.v)||key<=prev)throw Stop("literal_unique_order");
   fmpz_abs(ac.v,c.v);fmpz_add(next.v,L.v,ac.v);if(fmpz_bits(next.v)>301)throw Stop("input_L1_301_before_commit");own(1);R.emplace(key,c);++nr;fmpz_swap(next.v,L.v);prev=key;input_digest=digest(input_digest,key,c);
  }
  for(pass=0;pass<4;++pass){
   while(!R.empty()){
    phase=1;active_ready=0;bc=-1;auto lead=std::prev(R.end());int a,e,d;decode(lead->first,a,e,d);if(a<1||d<1)throw Stop("leading_monomial_not_divisible_ad");
    aqkey=lead->first-MA-1;box(aqkey,pass+1);if(!fmpz_divisible_ui(lead->second.v,50))throw Stop("integer_LC50_not_divisible");
    Z qc,back;fmpz_divexact_ui(qc.v,lead->second.v,50);bits(qc);fmpz_mul_ui(back.v,qc.v,50);++scalars;bits(back);if(!fmpz_equal(back.v,lead->second.v))throw Stop("LC50_backcheck");
    fmpz_swap(active.v,qc.v);active_ready=1;if(Q.find(aqkey)!=Q.end())throw Stop("nondecreasing_quotient_monomial");own(1);Q.emplace(aqkey,active);++nq;active_ready=2;terms[pass]=nq;stage_digest[pass]=digest(stage_digest[pass],aqkey,active);
    for(int j=0;j<6;++j){bc=j;if(((products+scalars)&4095)==0)clockcheck();Z v;fmpz_mul_si(v.v,active.v,-bcoef[j]);bits(v);update(R,aqkey+bkey[j],v,nr,pass);++products;bc=j+1;}
    active_ready=0;
   }
   phase=2;if(Q.empty())throw Stop("unexpected_zero_stage_quotient");if(pass<3){R.swap(Q);std::swap(nr,nq);}
  }
  // Q is final; no original vector survives. Build only fixed B^4.
  phase=3;Z one;fmpz_one(one.v);own(1);F.emplace(0,std::move(one));nf=1;
  for(fh=1;fh<=4;++fh){
   for(auto&[key,c]:F)for(int j=0;j<6;++j){fkey=key;bc=j;if(((products+scalars)&4095)==0)clockcheck();Z v;fmpz_mul_si(v.v,c.v,bcoef[j]);bits(v);update(N,key+bkey[j],v,nn,0,true);++products;}
   F.clear();nf=0;F.swap(N);std::swap(nf,nn);
  }
  // Re-read ONLY the immutable original4811 body, not a source oracle.
  phase=4;std::ifstream again(argv[2]);std::string m2,d2;size_t n2;int p2,s2;again>>m2>>n2>>p2>>s2>>d2;Z DD;fmpz_set_str(DD.v,d2.c_str(),10);
  if(m2!=magic||n2!=n||p2!=4||s2!=0||!fmpz_equal(DD.v,D.v))throw Stop("reread_header_D_changed");prev=-1;
  for(size_t i=0;i<n;++i){source_i=(int)i;long long key;if(!(again>>key>>s))throw std::runtime_error("reread_input_IO_failure");Z c,ac,next;fmpz_set_str(c.v,s.c_str(),10);box(key,0);if(!fmpz_sgn(c.v)||key<=prev)throw Stop("reread_unique_order");
   fmpz_abs(ac.v,c.v);fmpz_add(next.v,Lr.v,ac.v);if(fmpz_bits(next.v)>301)throw Stop("reread_L1_301_before_commit");own(1);R.emplace(key,c);++nr;fmpz_swap(next.v,Lr.v);prev=key;reread_digest=digest(reread_digest,key,c);
  }
  if(!fmpz_equal(Lr.v,L.v)||reread_digest!=input_digest)throw Stop("reread_original_mismatch");
  phase=5;
  for(auto&[qkey,qc]:Q)for(auto&[key,c]:F){vqkey=qkey;vfkey=key;if(((products+scalars)&4095)==0)clockcheck();Z v;fmpz_mul(v.v,qc.v,c.v);fmpz_neg(v.v,v.v);bits(v);update(R,qkey+key,v,nr,0);++products;}
  if(!R.empty())throw Stop("FULL_B4_multiplyback_NONZERO");phase=6;clockcheck();
  std::ofstream o(std::string(argv[3])+"/quotient-integer-common-D.txt");o<<"Q1B4QV001\n"<<zs(D)<<'\n'<<Q.size()<<'\n';for(auto&[key,c]:Q)o<<key<<' '<<zs(c)<<'\n';
  if(!o)throw std::runtime_error("output_IO_failure_no_second_fullprefix");
  std::ofstream r(std::string(argv[3])+"/native-report.json");r<<"{\"all_passed\":true,\"quotient_terms\":"<<Q.size()<<",\"all_four_stage_remainders_zero\":true,\"full_B4_multiplyback_zero\":true,\"original_input_reread_matched\":true,\"peak_owned\":"<<peak<<",\"integer_products\":"<<products<<",\"LC50_scalar_reversechecks\":"<<scalars<<",\"fixed_B4_terms\":"<<F.size()<<",\"stage_terms\":[";
  for(int j=0;j<4;++j){if(j)r<<',';r<<terms[j];}r<<"],\"stage_FNV64_telemetry_only\":[";for(int j=0;j<4;++j){if(j)r<<',';r<<stage_digest[j];}r<<"]}\n";return 0;
 }catch(const Stop&e){prefix(std::string(argv[3])+"/Q1B4QP001-prefix.txt",e.what(),D,L,Lr,active,R,Q,F,N);
  std::ofstream r(std::string(argv[3])+"/native-report.json");r<<"{\"all_passed\":false,\"prefix_only\":true,\"reason\":\""<<e.what()<<"\",\"peak_owned\":"<<peak<<",\"integer_products\":"<<products<<",\"LC50_scalar_reversechecks\":"<<scalars<<"}\n";return 10;
 }
}
