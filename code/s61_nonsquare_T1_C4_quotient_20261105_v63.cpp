// ONE fixed known-unit quotient. --help does not build powers or coefficients.
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
struct Rec{long long key;Z c;};using Poly=std::map<long long,Z>;
struct Stop:std::runtime_error{using std::runtime_error::runtime_error;};
static constexpr size_t OWN=81776,META=16000;
static constexpr long long MA=1048576,ME=1024;
static size_t inputN=0,rootN=0,constantN=0,qN=0,slabN=0,remN=0,peak=0;
static unsigned long long products=0,scalings=0,updates=0;
static int phase=0,qi=-1,ij=-1,rr=-1,ch=-1,ct=-1;
static long long transferkey=0,verifyqkey=0,verifyckey=0;
static std::chrono::steady_clock::time_point started;
static std::string zs(const Z&z){char*s=fmpz_get_str(nullptr,10,z.v);std::string r(s);flint_free(s);return r;}
static void bits(const Z&z){if(fmpz_bits(z.v)>518)throw Stop("518bit_before_commit_guard");}
static void clockcheck(){if(std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()>45)throw Stop("45s_soft_guard");}
static void own(size_t extra=0){size_t v=2*(inputN+rootN+constantN+qN+slabN+remN+extra)+META;peak=std::max(peak,v);if(v>OWN)throw Stop("81776_owned_before_commit_guard");}
static void update(Poly&p,long long key,const Z&delta,size_t&count){
 auto it=p.find(key);Z n;if(it==p.end())fmpz_set(n.v,delta.v);else fmpz_add(n.v,it->second.v,delta.v);bits(n);
 if(fmpz_is_zero(n.v)){if(it!=p.end()){p.erase(it);--count;}return;}
 if(it==p.end()){own(1);p.emplace(key,std::move(n));++count;}else fmpz_swap(it->second.v,n.v);++updates;
}
static void savepoly(std::ostream&o,const Poly&p){o<<p.size()<<'\n';for(auto&[key,c]:p)o<<key<<' '<<zs(c)<<'\n';}
static void prefix(const std::string&path,const std::string&why,const Z&D,const Z&L,const std::vector<Rec>&in,const std::vector<std::vector<Z>>&R,const std::vector<Poly>&C,const Poly&Q,const Poly&S,const Poly&M){
 std::ofstream o(path);o<<"T1C4QP001\n"<<phase<<' '<<qi<<' '<<ij<<' '<<rr<<' '<<ch<<' '<<ct<<' '<<transferkey<<' '<<verifyqkey<<' '<<verifyckey<<' '<<products<<' '<<scalings<<' '<<updates<<' '<<peak<<'\n'<<why<<'\n'<<zs(D)<<' '<<zs(L)<<'\n';
 o<<in.size()<<'\n';for(auto&r:in)o<<r.key<<' '<<zs(r.c)<<'\n';o<<R.size()<<'\n';for(auto&r:R){o<<r.size()<<'\n';for(auto&c:r)o<<zs(c)<<'\n';}
 o<<C.size()<<'\n';for(auto&p:C)savepoly(o,p);savepoly(o,Q);savepoly(o,S);savepoly(o,M);
 if(!o)throw std::runtime_error("prefix_output_device_failure");
}
static void decode(long long key,int&a,int&e,int&d){a=(int)(key/MA);long long v=key%MA;e=(int)(v/ME);d=(int)(v%ME);}
int main(int argc,char**argv){
 if(argc==2&&std::string(argv[1])=="--help"){
  std::cout<<"--run INPUT OUTDIR BASE_RAW. ONLY Q=T1/(a^3*C^4); full C4 multiplyback; owned81776,bits518,45soft60hard,RAW10610960.\n";return 0;
 }
 if(argc!=5||std::string(argv[1])!="--run")return 2;
 if(std::stoull(argv[4])+8410960+200000>10610960){std::cerr<<"RAW closure fails BEFORE arithmetic\n";return 3;}
 rlimit ram{2000000000,2000000000},cpu{60,60};setrlimit(RLIMIT_AS,&ram);setrlimit(RLIMIT_CPU,&cpu);started=std::chrono::steady_clock::now();
 std::ifstream f(argv[2]);std::string magic,str;size_t n;int am,shift,qm;f>>magic>>n>>am>>shift>>qm;
 if(magic!="T1C4Q001"||n!=7595||am!=20||shift!=3||qm!=16)return 4;
 Z D,L;f>>str;fmpz_set_str(D.v,str.c_str(),10);std::vector<Rec>in;std::vector<std::vector<Z>>R;std::vector<Poly>C;Poly Q,S,M;
 try{
  for(size_t i=0;i<n;++i){ij=(int)i;Rec r;if(!(f>>r.key>>str))throw std::runtime_error("literal_input_IO_failure");fmpz_set_str(r.c.v,str.c_str(),10);bits(r.c);
   int a,e,d;decode(r.key,a,e,d);if(a<3||a>23||e<0||e>41||d<0||d>18||a+e+d>61||fmpz_is_zero(r.c.v))throw Stop("literal_input_source_bounds");
   Z abs,next;fmpz_abs(abs.v,r.c.v);fmpz_add(next.v,L.v,abs.v);if(fmpz_bits(next.v)>363)throw Stop("literal_L1_363bit_guard");
   own(1);in.push_back(std::move(r));++inputN;fmpz_swap(L.v,next.v);
  }
  std::vector<Z>hundred(21);fmpz_one(hundred[0].v);for(int i=1;i<=20;++i){fmpz_mul_ui(hundred[i].v,hundred[i-1].v,100);++scalings;}
  unsigned long long choose[17];for(int h=0;h<=16;++h)choose[h]=(unsigned long long)(h+3)*(h+2)*(h+1)/6;
  phase=1;R.resize(17);for(int h=0;h<=16;++h){own(2*h+1);R[h].resize(2*h+1);rootN+=2*h+1;}fmpz_one(R[0][0].v);int root[3]={2,9,-18};
  for(int h=1;h<=16;++h)for(int r=0;r<(int)R[h-1].size();++r)for(int t=0;t<3;++t){ch=h;rr=r;ct=t;Z z,next;fmpz_mul_si(z.v,R[h-1][r].v,root[t]);++scalings;fmpz_add(next.v,R[h][r+t].v,z.v);bits(next);fmpz_swap(next.v,R[h][r+t].v);}
  for(int i=0;i<=16;++i){
   qi=i;phase=2;ij=rr=0;
   for(size_t j=0;j<in.size();++j){ij=(int)j;int a,e,d;decode(in[j].key,a,e,d);a-=3;if(a<i+4)continue;int h=a-i-4;rr=0;
    Z scale,base,z;fmpz_mul_ui(scale.v,hundred[16-h].v,choose[h]);fmpz_mul(base.v,in[j].c.v,scale.v);scalings+=2;bits(base);
    for(int r=0;r<(int)R[h].size();++r){rr=r;if((products&4095)==0)clockcheck();fmpz_mul(z.v,base.v,R[h][r].v);bits(z);++products;
     int ee=e+r;if(ee>73||d>18)throw Stop("partial_slab_Newton_bound");update(S,ee*64+d,z,slabN);
    }
   }
   // High-e/total-degree terms may cancel only after this WHOLE slab.
   phase=3;ij=rr=-1;
   for(auto&[key,c]:S){int e=(int)(key/64),d=(int)(key%64);if(e>33||i+e+d>50)throw Stop("completed_slab_outside_proved_Q_bounds");}
   while(!S.empty()){
    auto it=S.begin();transferkey=it->first;clockcheck();int e=(int)(it->first/64),d=(int)(it->first%64);
    update(Q,(long long)i*MA+e*ME+d,it->second,qN);S.erase(it);--slabN;
   }
  }
  if(Q.empty())throw Stop("unexpected_zero_quotient");
  // Release ALL R-power slots before allocating C^h/remainder.
  phase=4;R.clear();R.shrink_to_fit();rootN=0;C.resize(5);
  Z one;fmpz_one(one.v);update(C[0],0,one,constantN);
  int ckey[4]={64,2,1,0},ccoef[4]={100,18,-9,-2};
  for(int h=1;h<=4;++h)for(auto&[key,c]:C[h-1])for(int t=0;t<4;++t){ch=h;verifyckey=key;ct=t;Z z;fmpz_mul_si(z.v,c.v,ccoef[t]);++scalings;bits(z);update(C[h],key+ckey[t],z,constantN);}
  phase=5;
  for(size_t j=0;j<in.size();++j){ij=(int)j;clockcheck();Z z;fmpz_mul(z.v,in[j].c.v,hundred[20].v);++scalings;bits(z);update(M,in[j].key-3*MA,z,remN);}
  phase=6;ij=rr=-1;
  for(auto&[qkey,qc]:Q){verifyqkey=qkey;for(auto&[ck,cc]:C[4]){verifyckey=ck;if((products&4095)==0)clockcheck();
    int ca=(int)(ck/64),ce=(int)(ck%64);long long key=qkey+ca*MA+ce*ME;int a,e,d;decode(key,a,e,d);
    if(a>20||e>41||d>18||a+e+d>58)throw Stop("multiplyback_Newton_bounds");
    Z z;fmpz_mul(z.v,qc.v,cc.v);++products;fmpz_neg(z.v,z.v);bits(z);update(M,key,z,remN);
   }}
  if(!M.empty())throw Stop("FULL_C4_multiplyback_NONZERO");
  phase=7;clockcheck();Z Dq;fmpz_mul(Dq.v,D.v,hundred[20].v);bits(Dq);
  std::ofstream o(std::string(argv[3])+"/quotient-integer-common-D.txt");o<<"T1C4QV001\n"<<zs(Dq)<<'\n'<<Q.size()<<'\n';for(auto&[key,c]:Q)o<<key<<' '<<zs(c)<<'\n';
  if(!o)throw std::runtime_error("output_IO_failure_no_second_full_prefix");
  std::ofstream r(std::string(argv[3])+"/native-report.json");r<<"{\"all_passed\":true,\"quotient_terms\":"<<Q.size()<<",\"full_C4_multiplyback_zero\":true,\"all_completed_slabs_checked\":true,\"peak_owned\":"<<peak<<",\"polynomial_products\":"<<products<<",\"scaled_root_products\":"<<scalings<<",\"updates\":"<<updates<<"}\n";return 0;
 }catch(const Stop&e){
  prefix(std::string(argv[3])+"/T1C4QP001-prefix.txt",e.what(),D,L,in,R,C,Q,S,M);
  std::ofstream r(std::string(argv[3])+"/native-report.json");r<<"{\"all_passed\":false,\"prefix_only\":true,\"reason\":\""<<e.what()<<"\",\"peak_owned\":"<<peak<<",\"polynomial_products\":"<<products<<",\"scaled_root_products\":"<<scalings<<"}\n";return 10;
 }
}
