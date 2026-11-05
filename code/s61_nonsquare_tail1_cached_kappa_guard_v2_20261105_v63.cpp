// ONE fixed tail1 expression.  No arithmetic unless --run is explicitly used.
#include <flint/fmpz.h>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/resource.h>
struct Z {
 fmpz_t v;
 Z(){fmpz_init(v);} Z(const Z&o){fmpz_init(v);fmpz_set(v,o.v);}
 Z(Z&&o)noexcept{fmpz_init(v);fmpz_swap(v,o.v);}
 ~Z(){fmpz_clear(v);} Z&operator=(Z&&o)noexcept{fmpz_swap(v,o.v);return *this;}
};
using P=std::map<long long,Z>;
struct Node{int op=0,x=-1,y=-1;size_t sup=0;int bits=0;Z den;};
struct Stop:std::runtime_error{using std::runtime_error::runtime_error;};
static constexpr size_t OWN=91676, META=15100;
static constexpr long long MA=1048576,ME=1024;
static size_t terms=0,peak=0;static unsigned long long pairs=0,adds=0;
static int phase=0,atnode=-1,atpair=-1;static size_t ai=0,bi=0;
static std::chrono::steady_clock::time_point began;
static std::string zs(const Z&z){char*s=fmpz_get_str(nullptr,10,z.v);std::string r(s);flint_free(s);return r;}
static void clockcheck(){if(std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()>45)throw Stop("45s_soft_guard");}
static void capacity(size_t extra){size_t n=2*(terms+extra)+META;peak=std::max(peak,n);if(n>OWN)throw Stop("91676_owned_guard");}
static void keycheck(long long k){long long a=k/MA,r=k%MA;if(r<0){--a;r+=MA;}long long e=r/ME,d=r%ME;if(a<-128||a>128||e>128||d>128)throw Stop("packed_exponent_guard");}
static void commit(P&p,long long k,const Z&delta,int bitbound){
 keycheck(k);auto it=p.find(k);Z next;
 if(it==p.end())fmpz_set(next.v,delta.v);else fmpz_add(next.v,it->second.v,delta.v);
 if(fmpz_bits(next.v)>bitbound)throw Stop("coefficient_source_bit_guard");
 if(fmpz_is_zero(next.v)){if(it!=p.end()){p.erase(it);--terms;}return;}
 if(it==p.end()){capacity(1);p.emplace(k,std::move(next));++terms;}else fmpz_swap(it->second.v,next.v);
 ++adds;
}
static void release(P&p){terms-=p.size();p.clear();}
static void savepoly(std::ostream&o,const P&p){o<<p.size()<<'\n';for(auto&[k,c]:p)o<<k<<' '<<zs(c)<<'\n';}
static void prefix(const std::string&path,const std::vector<P>&v,const std::vector<int>&uses,const std::string&why){
 std::ofstream o(path);o<<"J1KP001\n"<<phase<<' '<<atnode<<' '<<atpair<<' '<<ai<<' '<<bi<<' '<<pairs<<' '<<adds<<' '<<peak<<'\n'<<why<<'\n';
 o<<v.size()<<'\n';for(size_t i=0;i<v.size();++i){o<<i<<' '<<(i<uses.size()?uses[i]:0)<<'\n';savepoly(o,v[i]);}
 if(!o)throw std::runtime_error("prefix_write_failed");
}
int main(int argc,char**argv){
 if(argc==2&&std::string(argv[1])=="--help"){
  std::cout<<"--run INPUT OUTPUT_DIRECTORY BASE_RAW_BYTES. Fixed tail1; J1K001; own<=91676, coeff<=359bits,45ssoft/60shard,9517864 totalRAW. No auto-resume.\n";return 0;
 }
 if(argc!=5||std::string(argv[1])!="--run"){std::cerr<<"Use --help. No arithmetic started.\n";return 2;}
 // 38288 live coefficient/key pairs, each at most109+10+3 bytes,
 // plus all155 headers/cursors: prefix<=4.70MB. No double buffer.
 if(std::stoull(argv[4])+4700000+100000>9517864){std::cerr<<"RAW closure fails before arithmetic\n";return 11;}
 rlimit ram{2000000000,2000000000};setrlimit(RLIMIT_AS,&ram);
 rlimit cpu{60,60};setrlimit(RLIMIT_CPU,&cpu);
 began=std::chrono::steady_clock::now();std::string out=argv[3];std::ifstream in(argv[2]);
 std::string magic;int nn,ne;in>>magic>>nn>>ne;if(magic!="J1K001"||nn!=155||ne!=149)return 3;
 std::vector<Node>N(nn);std::vector<P>V(nn+1);std::vector<int>order(ne),uses(nn+1);
 int col[4],kap[4];for(int&i:col)in>>i;for(int&i:kap)in>>i;
 Z D;std::string s;in>>s;fmpz_set_str(D.v,s.c_str(),10);
 for(int&i:order)in>>i;
 for(int i=0;i<nn;++i){int id;in>>id>>N[i].op>>N[i].x>>N[i].y>>N[i].sup>>N[i].bits>>s;
  if(id!=i)return 4;fmpz_set_str(N[i].den.v,s.c_str(),10);
 }
 int na;in>>na;
 try{
  for(int t=0;t<na;++t){int id;size_t m;in>>id>>m;phase=0;atnode=id;atpair=t;ai=bi=0;if(id<0||id>=nn||N[id].op!=0)throw Stop("literal_atom_id");
   for(size_t j=0;j<m;++j){long long k;ai=j;in>>k>>s;Z c;fmpz_set_str(c.v,s.c_str(),10);if(fmpz_bits(c.v)>N[id].bits)throw Stop("literal_atom_bits");commit(V[id],k,c,N[id].bits);}
   if(V[id].size()>N[id].sup)throw Stop("atom_support");
  }
  if(!in)throw Stop("input_truncated");
  for(int i:order)if(N[i].op){++uses[N[i].x];++uses[N[i].y];}
  for(int i:col)++uses[i];for(int i:kap)++uses[i];
  phase=1;
  for(int i:order){atnode=i;atpair=-1;ai=bi=0;if(!N[i].op)continue;
   int x=N[i].x,y=N[i].y;Z a,b,t,dcheck;
   if(N[i].op==3){
    fmpz_mul(dcheck.v,N[x].den.v,N[y].den.v);if(!fmpz_equal(dcheck.v,N[i].den.v))throw Stop("mul_D_identity");
    for(auto&[kx,cx]:V[x]){bi=0;for(auto&[ky,cy]:V[y]){
      if((pairs&4095)==0)clockcheck();fmpz_mul(t.v,cx.v,cy.v);commit(V[i],kx+ky,t,N[i].bits);++pairs;++bi;
     }++ai;
    }
   }else{
    if(!fmpz_divisible(N[i].den.v,N[x].den.v)||!fmpz_divisible(N[i].den.v,N[y].den.v))throw Stop("add_D_divexact");
    fmpz_divexact(a.v,N[i].den.v,N[x].den.v);fmpz_divexact(b.v,N[i].den.v,N[y].den.v);
    if(N[i].op==2)fmpz_neg(b.v,b.v);
    for(auto&[k,c]:V[x]){fmpz_mul(t.v,c.v,a.v);commit(V[i],k,t,N[i].bits);++ai;}
    bi=0;for(auto&[k,c]:V[y]){fmpz_mul(t.v,c.v,b.v);commit(V[i],k,t,N[i].bits);++bi;}
   }
   if(V[i].size()>N[i].sup)throw Stop("node_support");
   if(!--uses[x])release(V[x]);if(!--uses[y])release(V[y]);
   // -2 explicitly marks this node COMPLETE and its uses already released.
   atpair=-2;ai=bi=0;clockcheck();
  }
  atnode=nn;P&acc=V[nn];
  // Three deterministic passes: forward; reverse subtraction to literal zero;
  // forward reproduction for final publication. No second accumulator map.
  for(int pass=0;pass<3;++pass){phase=2+pass;
   for(int ii=0;ii<4;++ii){int j=pass==1?3-ii:ii;atpair=j;ai=bi=0;Z dd,scale,t,w;
    fmpz_mul(dd.v,N[kap[j]].den.v,N[col[j]].den.v);
    if(!fmpz_divisible(D.v,dd.v))throw Stop("dot_D_divexact");fmpz_divexact(scale.v,D.v,dd.v);if(pass==1)fmpz_neg(scale.v,scale.v);
    if(pass==1){
     for(auto ix=V[kap[j]].rbegin();ix!=V[kap[j]].rend();++ix){bi=0;for(auto iy=V[col[j]].rbegin();iy!=V[col[j]].rend();++iy){
       if((pairs&4095)==0)clockcheck();fmpz_mul(t.v,ix->second.v,iy->second.v);fmpz_mul(w.v,t.v,scale.v);commit(acc,ix->first+iy->first,w,359);++pairs;++bi;
      }++ai;}
    }else{
     for(auto&[kx,cx]:V[kap[j]]){bi=0;for(auto&[ky,cy]:V[col[j]]){
       if((pairs&4095)==0)clockcheck();fmpz_mul(t.v,cx.v,cy.v);fmpz_mul(w.v,t.v,scale.v);commit(acc,kx+ky,w,359);++pairs;++bi;
      }++ai;}
    }
   }
   if(pass==1&&!acc.empty())throw Stop("reverse_all_four_products_nonzero");
  }
  phase=5;clockcheck();
  std::ofstream o(out+"/integer-column4-and-result.txt");o<<"J1KV001\n";
  for(int j=0;j<4;++j){o<<col[j]<<' '<<zs(N[col[j]].den)<<'\n';savepoly(o,V[col[j]]);}
  o<<nn<<' '<<zs(D)<<'\n';savepoly(o,acc);
  // An output-device failure is NOT a mathematical guard and must not
  // create a second full prefix beside a partially written success file.
  if(!o)throw std::runtime_error("output_device_failure_no_second_full_prefix");
  std::ofstream r(out+"/native-report.json");r<<"{\"all_passed\":true,\"zero\":"<<(acc.empty()?"true":"false")<<",\"result_terms\":"<<acc.size()<<",\"peak_owned\":"<<peak<<",\"pairs\":"<<pairs<<",\"adds\":"<<adds<<",\"reverse_all_four_products_zero\":true,\"three_passes_complete\":true}\n";
  return 0;
 }catch(const Stop&e){
  prefix(out+"/J1KP001-prefix.txt",V,uses,e.what());
  std::ofstream r(out+"/native-report.json");r<<"{\"all_passed\":false,\"prefix_only\":true,\"reason\":\""<<e.what()<<"\",\"peak_owned\":"<<peak<<",\"pairs\":"<<pairs<<"}\n";return 10;
 }
}
