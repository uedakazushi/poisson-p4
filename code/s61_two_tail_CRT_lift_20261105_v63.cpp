// ONE fixed40-prime candidate. No values before the separate preparation audit/GO.
#include "s61_two_tail_CRT_state_20261105_v63.hpp"
#include <sys/resource.h>
#include <algorithm>
#include <iostream>
#include <sstream>
using namespace tt;
static void guardInput(State&s,unsigned k,uint32_t key){unsigned a,e,d;xyz(key,a,e,d);unsigned A=k?12:16,E=k?25:32,Dd=k?14:18,T=k?38:49,W=k?35:42;if(a>A||e>E||d>Dd||a+e+d>T||2*a+e>W)throw Stop("original_input_box");}
static void readInput(State&s,unsigned k,const std::string&path){
 s.u[MODE]=LOAD;s.u[WHICH]=k;s.u[IDX]=0;std::ifstream f(path);std::string magic,Dstr,cstr;size_t n;f>>magic>>Dstr>>n;
 if(magic!=(k?"Q1B4QV001":"T0C4QV001")||n!=(k?2317:4500))throw Stop("input_header");
 fmpz_set_str(s.z[TMP].v,Dstr.c_str(),10);if(fmpz_sgn(s.z[TMP].v)<=0||fmpz_bits(s.z[TMP].v)!=197)throw Stop("input_D197");
 if(k){if(!fmpz_equal(s.z[D].v,s.z[TMP].v))throw Stop("two_Ds_differ");}else fmpz_set(s.z[D].v,s.z[TMP].v);
 uint32_t prev=0;bool first=true;unsigned dl=k?INDEG2:INDEG0,lc=k?LC2:LC0,ls=k?L12:L10;s.u[dl]=0;
 for(size_t i=0;i<n;++i){s.u[IDX]=i;s.u[INT_STAGE]=0;uint32_t key;if(!(f>>key>>cstr))throw Stop("input_truncated");s.u[ACTIVE_KEY]=key;guardInput(s,k,key);if(!first&&key<=prev)throw Stop("input_order");first=false;prev=key;
  fmpz_set_str(s.z[ACTIVE].v,cstr.c_str(),10);s.u[INT_STAGE]=1;if(!fmpz_sgn(s.z[ACTIVE].v))throw Stop("input_zero_record");s.bit(s.z[ACTIVE],301);
  fmpz_abs(s.z[TMP].v,s.z[ACTIVE].v);fmpz_add(s.z[PENDING].v,s.z[ls].v,s.z[TMP].v);s.u[INT_STAGE]=2;s.bit(s.z[PENDING],k?273:298);
  Rec r;r.key=key;fmpz_set(r.c.v,s.z[ACTIVE].v);s.in[k].push_back(std::move(r));fmpz_swap(s.z[ls].v,s.z[PENDING].v);s.u[IDX]=i+1;
  unsigned E=kron(key);if(E>=s.u[dl]){s.u[dl]=E;fmpz_set(s.z[lc].v,s.z[ACTIVE].v);}if((i&1023)==0)s.clock();
 }
 if(s.u[dl]!=(k?7975:10544))throw Stop("input_degree");if(f>>cstr)throw Stop("input_extra");
}
static void seedCRT(State&s,const std::string&path){
 s.u[MODE]=SEED;fmpz_gcd(s.z[GAMMA].v,s.z[LC0].v,s.z[LC2].v);fmpz_abs(s.z[GAMMA].v,s.z[GAMMA].v);if(fmpz_bits(s.z[GAMMA].v)!=226)throw Stop("gamma226");
 s.own(3900);s.crt.resize(3900);fmpz_set_ui(s.z[MODULUS].v,P[0]);uint32_t gm=fmpz_fdiv_ui(s.z[GAMMA].v,P[0]);if(!gm)throw Stop("seed_gamma_notunit");
 std::ifstream f(path);unsigned degree;f>>degree;if(degree!=3899)throw Stop("seed_degree");unsigned E,c,prev=0,count=0;bool first=true;
 while(f>>E>>c){s.u[IDX]=count;if(E>3899||!c||c>=P[0]||(!first&&E<=prev))throw Stop("seed_body");first=false;prev=E;unsigned a=E/646,r=E%646,e=r/19,d=r%19;
  if(a>6||e>9||d>4||a+e+d>11||2*a+e>13||(2*a+e==13&&(a!=6||e!=1||d!=4)))throw Stop("seed_support");
  fmpz_set_ui(s.crt[E].v,s.mul(gm,c,P[0]));++count;
 }if(count!=94||prev!=3899||fmpz_fdiv_ui(s.crt.back().v,P[0])!=gm)throw Stop("seed_94_or_lead");
 s.u[ACCEPTED]=1;s.ledger[0]={P[0],1,3899,s.u[PRIMEOPS],1,15,fmpz_fdiv_ui(s.z[LC0].v,P[0]),fmpz_fdiv_ui(s.z[LC2].v,P[0])};
}
static bool onePrime(State&s,unsigned pi){
 s.u[PI]=pi;s.u[MODE]=FIELD;s.u[PRIMEOPS]=0;s.u[EUPASS]=0;s.u[EUST]=100;unsigned p=P[pi];
 release(s.f);release(s.g);s.ledger[pi][0]=p;s.ledger[pi][6]=fmpz_fdiv_ui(s.z[LC0].v,p);s.ledger[pi][7]=fmpz_fdiv_ui(s.z[LC2].v,p);
 if(!s.ledger[pi][6]||!s.ledger[pi][7]||!fmpz_fdiv_ui(s.z[GAMMA].v,p)||!fmpz_fdiv_ui(s.z[D].v,p)){s.ledger[pi][1]=2;return false;}
 s.own(18521);std::vector<uint32_t>(10545).swap(s.f);std::vector<uint32_t>(7976).swap(s.g);s.own();
 for(unsigned k=0;k<2;++k){s.u[WHICH]=k;auto&v=k?s.g:s.f;for(size_t i=0;i<s.in[k].size();++i){s.u[IDX]=i;auto&r=s.in[k][i];v[kron(r.key)]=fmpz_fdiv_ui(r.c.v,p);s.u[IDX]=i+1;if((i&1023)==0)s.clock();}}
 if(!s.f.back()||!s.g.back())throw Stop("prime_input_degree_drop_after_LCcheck");s.u[EUST]=0;
 while(!s.g.empty()){
  s.u[EUST]=0;s.clock();s.u[IL]=s.g.back();s.u[IL]=s.inverse(s.u[IL],p);s.u[EUST]=1;
  while(!s.f.empty()&&s.f.size()>=s.g.size()){
   s.u[SHIFT]=s.f.size()-s.g.size();s.u[FACTOR]=s.mul(s.f.back(),s.u[IL],p);s.u[JNEXT]=0;s.u[EUST]=2;
   for(size_t j=0;j<s.g.size();++j){s.u[JNEXT]=j;if((s.u[MODOPS]&4095)==0)s.clock();unsigned v=s.mul(s.u[FACTOR],s.g[j],p);s.f[s.u[SHIFT]+j]=(s.f[s.u[SHIFT]+j]+p-v)%p;s.u[JNEXT]=j+1;}
   s.u[EUST]=3;trim(s.f);s.u[EUST]=1;
  }s.f.swap(s.g);++s.u[EUPASS];s.u[EUST]=0;s.clock();
 }
 if(s.f.empty())throw Stop("unexpected_zero_gcd");unsigned degree=s.f.size()-1;s.ledger[pi][2]=degree;s.ledger[pi][3]=s.u[PRIMEOPS];
 if(degree<3899){s.ledger[pi][1]=4;throw Stop("gcd_degree_below_fixed3899_no_reset");}
 if(degree>3899){s.ledger[pi][1]=3;return false;}
 s.u[MODE]=MONIC;s.u[NORM_STAGE]=0;s.u[NORMI]=0;s.u[IL]=s.f.back();s.u[IL]=s.inverse(s.u[IL],p);s.u[NORM_STAGE]=1;
 for(size_t i=0;i<s.f.size();++i){s.u[NORMI]=i;if((i&1023)==0)s.clock();s.f[i]=s.mul(s.f[i],s.u[IL],p);s.u[NORMI]=i+1;}
 s.u[MODE]=CRT;s.u[NORM_STAGE]=0;s.u[CRSUB]=0;s.u[CRINV]=fmpz_fdiv_ui(s.z[MODULUS].v,p);s.u[CRINV]=s.inverse(s.u[CRINV],p);s.u[NORM_STAGE]=1;
 s.bigreserve();fmpz_mul_ui(s.z[NEWMOD].v,s.z[MODULUS].v,p);++s.u[BIGOPS];s.bit(s.z[NEWMOD],600);
 s.u[NORM_STAGE]=2;unsigned gm=fmpz_fdiv_ui(s.z[GAMMA].v,p);
 for(size_t i=0;i<3900;++i){s.u[CRTI]=i;s.u[CRSUB]=0;if((i&255)==0)s.clock();s.u[CRRES]=s.mul(gm,s.f[i],p);s.u[CRSUB]=1;
  unsigned cm=fmpz_fdiv_ui(s.crt[i].v,p);s.u[CRU]=s.mul((s.u[CRRES]+p-cm)%p,s.u[CRINV],p);s.u[CRSUB]=2;
  s.bigreserve();fmpz_mul_ui(s.z[PENDING].v,s.z[MODULUS].v,s.u[CRU]);fmpz_add(s.z[PENDING].v,s.z[PENDING].v,s.crt[i].v);++s.u[BIGOPS];s.u[CRSUB]=3;s.bit(s.z[PENDING],600);
  fmpz_swap(s.crt[i].v,s.z[PENDING].v);s.u[CRSUB]=0;s.u[CRTI]=i+1;
 }fmpz_swap(s.z[MODULUS].v,s.z[NEWMOD].v);s.u[CRSUB]=4;++s.u[ACCEPTED];s.ledger[pi][1]=1;s.ledger[pi][3]=s.u[PRIMEOPS];s.ledger[pi][4]=s.u[ACCEPTED];s.ledger[pi][5]=fmpz_bits(s.z[MODULUS].v);return true;
}
static void normalizeCandidate(State&s,const std::string&out){
 release(s.f);release(s.g);s.u[MODE]=BALANCE;fmpz_mul_ui(s.z[TMP].v,s.z[GAMMA].v,2);if(fmpz_cmp(s.z[MODULUS].v,s.z[TMP].v)<=0)throw Stop("CRT_modulus_not_above_2gamma");
 fmpz_fdiv_q_2exp(s.z[TMP].v,s.z[MODULUS].v,1);
 for(size_t i=0;i<3900;++i){s.u[IDX]=i;s.u[NORM_STAGE]=0;if((i&255)==0)s.clock();if(fmpz_cmp(s.crt[i].v,s.z[TMP].v)>0){fmpz_sub(s.z[PENDING].v,s.crt[i].v,s.z[MODULUS].v);s.u[NORM_STAGE]=1;s.bit(s.z[PENDING],600);fmpz_swap(s.crt[i].v,s.z[PENDING].v);}s.u[IDX]=i+1;}
 if(!fmpz_equal(s.crt.back().v,s.z[GAMMA].v))throw Stop("CRT_leading_not_exact_gamma");
 s.u[MODE]=CONTENT_GCD;fmpz_zero(s.z[CONTENT].v);
 for(size_t i=0;i<3900;++i){s.u[IDX]=i;s.u[NORM_STAGE]=0;if((i&255)==0)s.clock();fmpz_gcd(s.z[PENDING].v,s.z[CONTENT].v,s.crt[i].v);s.u[NORM_STAGE]=1;s.bit(s.z[PENDING],600);fmpz_swap(s.z[CONTENT].v,s.z[PENDING].v);s.u[IDX]=i+1;}
 if(fmpz_sgn(s.z[CONTENT].v)<=0)throw Stop("zero_integer_content");saveCRT(s,out+"/CRT-balanced.bin");s.u[HAS_CRT]=1;
 s.u[MODE]=PRIMITIVE;
 for(size_t i=0;i<3900;++i){s.u[IDX]=i;s.u[NORM_STAGE]=0;if((i&255)==0)s.clock();if(!fmpz_divisible(s.crt[i].v,s.z[CONTENT].v))throw Stop("content_not_exact_integer_divisor");fmpz_divexact(s.z[PENDING].v,s.crt[i].v,s.z[CONTENT].v);s.u[NORM_STAGE]=1;s.bigreserve();fmpz_mul(s.z[TMP].v,s.z[PENDING].v,s.z[CONTENT].v);++s.u[BIGOPS];s.u[NORM_STAGE]=2;s.bit(s.z[TMP],600);if(!fmpz_equal(s.z[TMP].v,s.crt[i].v))throw Stop("content_exact_reversecheck");fmpz_swap(s.crt[i].v,s.z[PENDING].v);s.u[IDX]=i+1;}
 if(fmpz_sgn(s.crt.back().v)<=0)throw Stop("primitive_leading_not_positive");
 s.u[MODE]=CANDIDATE;s.u[NORM_STAGE]=0;unsigned count=0,A=0,E=0,Dd=0,T=0;
 for(size_t i=0;i<3900;++i){s.u[IDX]=i;if(!fmpz_sgn(s.crt[i].v))continue;unsigned a=i/646,r=i%646,e=r/19,d=r%19;
  if(a>6||e>9||d>4||a+e+d>11||2*a+e>13||(2*a+e==13&&(a!=6||e!=1||d!=4)))throw Stop("selected_candidate_outside_fixed_shape_no_growth");
  ++count;A=std::max(A,a);E=std::max(E,e);Dd=std::max(Dd,d);T=std::max(T,a+e+d);
 }if(count>207||A!=6||E!=9||Dd!=4||T!=11)throw Stop("candidate_degree_profile_or207");
 s.u[NORM_STAGE]=1;for(size_t i=0;i<3900;++i){s.u[IDX]=i;if(!fmpz_sgn(s.crt[i].v))continue;unsigned a=i/646,r=i%646,e=r/19,d=r%19;uint64_t key=weighted(a*MA+e*ME+d);s.own(2);s.H.emplace(key,s.crt[i]);s.u[IDX]=i+1;}
 if(s.H.empty()||s.H.rbegin()->first!=weighted(6*MA+ME+4))throw Stop("selected_factor_LM6e1d4");
 fmpz_set(s.z[LC].v,s.H.rbegin()->second.v);fmpz_zero(s.z[NORM].v);
 s.u[NORM_STAGE]=2;for(auto&[key,c]:s.H){s.u[HKEY]=key;fmpz_abs(s.z[TMP].v,c.v);fmpz_add(s.z[PENDING].v,s.z[NORM].v,s.z[TMP].v);s.bit(s.z[PENDING],608);fmpz_swap(s.z[NORM].v,s.z[PENDING].v);}
 saveMap(s.H,out+"/factor.bin","TTFAC001",2,1024*1024,600);s.u[HAS_FACTOR]=1;s.u[NORM_STAGE]=3;
 std::ofstream txt(out+"/factor-sparse.txt");txt<<s.H.size()<<'\n';for(auto&[key,c]:s.H){unsigned a,e,d;xyz(key%WRAD,a,e,d);txt<<a<<' '<<e<<' '<<d<<' '<<zs(c)<<'\n';}if(!txt)throw Stop("factor_text_IO_failure");
 s.u[MODE]=HEIGHT;fmpz_sub(s.z[PENDING].v,s.z[NORM].v,s.z[LC].v);unsigned r=0;
 if(fmpz_sgn(s.z[PENDING].v)>0){int d=int(fmpz_bits(s.z[PENDING].v))-int(fmpz_bits(s.z[LC].v));r=std::max(0,d);fmpz_mul_2exp(s.z[TMP].v,s.z[LC].v,r);if(fmpz_cmp(s.z[PENDING].v,s.z[TMP].v)>0)++r;}
 s.u[REXP]=r;unsigned bh=fmpz_bits(s.z[NORM].v),fl=fmpz_bits(s.z[LC].v)-1;
 for(unsigned k=0;k<2;++k){int b=int(fmpz_bits(s.z[k?L12:L10].v))-int(fl)+5+int((k?22:29)*r);unsigned qb=std::max(0,b),rb=std::max<unsigned>(fmpz_bits(s.z[k?L12:L10].v),qb+bh)+1;s.u[k?QB2:QB0]=qb;s.u[k?RB2:RB0]=rb;}
 if(s.u[RB0]>8192||s.u[RB2]>8192)throw Stop("conditional8192_height_gate_failed_no_division");
 std::vector<Z>().swap(s.crt);if(s.f.capacity()||s.g.capacity()||s.crt.capacity())throw Stop("field_CRT_not_released_before_integer_maps");s.own();
}
static void workBox(uint64_t key,unsigned which,bool quotient){
 unsigned a,e,d;xyz(key%WRAD,a,e,d);uint64_t w=key/WRAD;if(w!=2*a+e||e>=1024||d>=1024)throw Stop("weighted_key_alias");
 unsigned A=which?12:16,Dd=which?14:18,W=which?35:42;if(quotient){A-=6;Dd-=4;W-=13;}if(a>A||d>Dd||w>W)throw Stop("false_candidate_safe_partial_support");
}
static void addR(State&s,uint64_t key){
 workBox(key,s.u[WHICH],false);auto it=s.R.find(key);if(it==s.R.end())fmpz_set(s.z[PENDING].v,s.z[TMP].v);else fmpz_add(s.z[PENDING].v,it->second.v,s.z[TMP].v);s.u[INT_STAGE]=3;s.bit(s.z[PENDING],8192);
 if(!fmpz_sgn(s.z[PENDING].v)){if(it!=s.R.end())s.R.erase(it);}else if(it==s.R.end()){s.own(2);Z c;fmpz_set(c.v,s.z[PENDING].v);s.R.emplace(key,std::move(c));}else fmpz_swap(it->second.v,s.z[PENDING].v);s.u[INT_STAGE]=1;
}
static void loadR(State&s,unsigned k,bool verify){
 s.u[MODE]=verify?VERIFY_LOAD:INT_LOAD;s.u[WHICH]=k;if(!s.R.empty())throw Stop("R_not_empty_before_load");
 for(size_t i=0;i<s.in[k].size();++i){s.u[IDX]=i;auto&r=s.in[k][i];uint64_t key=weighted(r.key);s.u[ACTIVE_KEY]=key;workBox(key,k,false);s.own(2);s.R.emplace(key,r.c);s.u[IDX]=i+1;if((i&1023)==0)s.clock();}
}
static void integerIdentity(State&s,unsigned k,const std::string&out){
 loadR(s,k,false);s.u[MODE]=INT_DIV;s.u[ACTIVE_READY]=0;uint64_t lead=s.H.rbegin()->first;
 while(!s.R.empty()){
  s.u[INT_STAGE]=0;s.u[ACTIVE_READY]=0;s.clock();auto it=std::prev(s.R.end());unsigned a,e,d;xyz(it->first%WRAD,a,e,d);if(a<6||e<1||d<4)throw Stop("candidate_not_dividing_leading_monomial");
  uint64_t qkey=it->first-lead;s.u[ACTIVE_KEY]=qkey;workBox(qkey,k,true);if(!fmpz_divisible(it->second.v,s.z[LC].v))throw Stop("primitive_candidate_not_integer_LCdivisor");
  fmpz_divexact(s.z[ACTIVE].v,it->second.v,s.z[LC].v);s.u[ACTIVE_READY]=1;s.bit(s.z[ACTIVE],s.u[k?QB2:QB0]);
  s.bigreserve();fmpz_mul(s.z[TMP].v,s.z[ACTIVE].v,s.z[LC].v);++s.u[BIGOPS];s.u[INT_STAGE]=2;s.bit(s.z[TMP],8192);if(!fmpz_equal(s.z[TMP].v,it->second.v))throw Stop("integer_LC_reversecheck");
  if(s.Q.count(qkey))throw Stop("quotient_order_duplicate");s.own(2);s.Q.emplace(qkey,s.z[ACTIVE]);s.u[ACTIVE_READY]=2;s.u[HINDEX]=0;
  for(auto&[hkey,hc]:s.H){s.u[HKEY]=hkey;s.u[INT_STAGE]=0;if((s.u[BIGOPS]&255)==0)s.clock();s.bigreserve();fmpz_mul(s.z[TMP].v,s.z[ACTIVE].v,hc.v);fmpz_neg(s.z[TMP].v,s.z[TMP].v);++s.u[BIGOPS];s.u[INT_STAGE]=1;s.bit(s.z[TMP],8192);addR(s,qkey+hkey);++s.u[HINDEX];}
  s.u[ACTIVE_READY]=0;
 }
 loadR(s,k,true);s.u[MODE]=VERIFY_MUL;
 for(auto&[qkey,qc]:s.Q){s.u[QKEY]=qkey;s.u[HINDEX]=0;for(auto&[hkey,hc]:s.H){s.u[HKEY]=hkey;s.u[INT_STAGE]=0;if((s.u[BIGOPS]&255)==0)s.clock();s.bigreserve();fmpz_mul(s.z[TMP].v,qc.v,hc.v);fmpz_neg(s.z[TMP].v,s.z[TMP].v);++s.u[BIGOPS];s.u[INT_STAGE]=1;s.bit(s.z[TMP],8192);addR(s,qkey+hkey);++s.u[HINDEX];}}
 if(!s.R.empty())throw Stop("ALL_multivariate_multiplyback_NONZERO");s.u[k?VERIFIED2:VERIFIED0]=1;
 s.u[INT_STAGE]=99;saveMap(s.Q,out+(k?"/cofactor-N2.bin":"/cofactor-N0.bin"),"TTCOF001",k,(uint64_t(k?1309:3300)*1040)+4096,8192);s.Q.clear();s.own();
}
static void report(State&s,const std::string&out,bool ok,const std::string&why,bool pref){
 std::ofstream r(out+"/native-report.json");r<<"{\"all_passed\":"<<(ok?"true":"false")<<",\"common_char0_factor_proved\":"<<(ok?"true":"false")<<",\"prefix_saved\":"<<(pref?"true":"false")<<",\"reason\":\""<<why<<"\",\"mode\":"<<s.u[MODE]<<",\"completed_primes\":"<<s.u[ACCEPTED]<<",\"factor_terms\":"<<s.H.size()<<",\"expected_Kronecker_degree\":3899,\"expected_LM\":[6,1,4],\"candidate_LM_checked\":"<<(s.u[HAS_FACTOR]?"true":"false")<<",\"verified_originals\":["<<s.u[VERIFIED0]<<','<<s.u[VERIFIED2]<<"],\"height_bits\":["<<s.u[RB0]<<','<<s.u[RB2]<<"],\"modular_products\":"<<s.u[MODOPS]<<",\"bigint_products\":"<<s.u[BIGOPS]<<",\"peak_owned\":"<<s.u[PEAK]<<",\"no_irreducibility_maximality_Fitting_or_fullrank_claim\":true}\n";if(!r)throw Stop("native_report_IO_failure");
}
int main(int argc,char**argv){
 if(argc==2&&std::string(argv[1])=="--help"){std::cout<<"--run N0 N2 SEED OUTDIR BASE_RAW; fixed40/reuse32003,degree3899/noreset,primitivecandidate/conditional8192 weighteddivision/ALL2identities;60000owned/6Bmod/2.1Mbig/45soft60hard/2GB/RAW32MiB. --prefix-fixture PATH uses synthetic IO data ONLY.\n";return 0;}
 if(argc==3&&std::string(argv[1])=="--prefix-fixture"){
  State s;s.started=std::chrono::steady_clock::now();s.u[MODE]=FIXTURE;s.in[0].reserve(1);Rec r;r.key=0;fmpz_set_si(r.c.v,-3);s.in[0].push_back(std::move(r));s.crt.resize(3);fmpz_one(s.crt[1].v);fmpz_mul_2exp(s.crt[1].v,s.crt[1].v,599);fmpz_add_ui(s.crt[1].v,s.crt[1].v,7);Z x;fmpz_one(x.v);fmpz_mul_2exp(x.v,x.v,8191);fmpz_add_ui(x.v,x.v,7);s.R.emplace(weighted(16*MA+10*ME+18),std::move(x));fmpz_one(s.z[PENDING].v);fmpz_mul_2exp(s.z[PENDING].v,s.z[PENDING].v,8192);s.own();savePrefix(s,argv[2],"SYNTHETIC_IO_ONLY_NO_TARGET_ARITHMETIC");return 0;
 }
 if(argc!=7||std::string(argv[1])!="--run")return 2;uint64_t base=std::stoull(argv[6]);if(base>3*1024*1024||base+23ULL*1024*1024>32ULL*1024*1024)return 3;
 rlimit ram{2000000000,2000000000},cpu{60,60};setrlimit(RLIMIT_AS,&ram);setrlimit(RLIMIT_CPU,&cpu);State s;s.started=std::chrono::steady_clock::now();s.paths={argv[2],argv[3],argv[4]};std::string out=argv[5];bool attemptedPrefix=false;
 try{
  s.own(2*6817);s.in[0].reserve(4500);s.in[1].reserve(2317);readInput(s,0,argv[2]);readInput(s,1,argv[3]);seedCRT(s,argv[4]);std::ofstream log(out+"/prime-ledger.jsonl");
  for(unsigned pi=1;pi<40;++pi){onePrime(s,pi);log<<"{\"p\":"<<P[pi]<<",\"status\":"<<s.ledger[pi][1]<<",\"degree\":"<<s.ledger[pi][2]<<",\"prime_products\":"<<s.u[PRIMEOPS]<<",\"Mbits\":"<<fmpz_bits(s.z[MODULUS].v)<<"}\n";log.flush();if(!log)throw Stop("prime_ledger_IO_failure");release(s.f);release(s.g);s.u[PI]=pi+1;s.clock();}
  normalizeCandidate(s,out);integerIdentity(s,0,out);integerIdentity(s,1,out);s.u[MODE]=DONE;report(s,out,true,"complete_exact_TWO_multivariate_identities",false);return 0;
 }catch(const std::exception&e){
  std::string why=e.what();for(char&c:why)if(c=='"'||c=='\\'||c=='\n')c='_';bool saved=false;
  if(!attemptedPrefix){attemptedPrefix=true;try{savePrefix(s,out+"/TTCRT001-prefix.bin",why);saved=true;}catch(const std::exception&io){std::cerr<<"PREFIX_IO_FAILURE_NO_SECOND_STATE: "<<io.what()<<'\n';}}
  try{report(s,out,false,why,saved);}catch(const std::exception&io){std::cerr<<"REPORT_IO_FAILURE_NO_RERUN: "<<io.what()<<'\n';}return 10;
 }
}
