// Verify given sparse monic modular gcd by exact divisions and quotient coprimality.
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <stdexcept>
#include <chrono>
#include <sys/resource.h>
using V=std::vector<unsigned>;const unsigned P=32003;V A,B,R,T;std::vector<std::pair<int,unsigned>> H;unsigned long long ops=0;int phase=0,cur=-1,pass=0;auto begun=std::chrono::steady_clock::now();
void own(size_t extra=0){if(A.capacity()+B.capacity()+R.capacity()+T.capacity()+2*H.capacity()+extra+4096>30000)throw std::runtime_error("30000_owned_guard");}
unsigned mul(unsigned a,unsigned b){if(++ops>150000000)throw std::runtime_error("150M_product_guard");return(a*b)%P;}
void tick(){if(std::chrono::duration<double>(std::chrono::steady_clock::now()-begun).count()>45)throw std::runtime_error("45s_soft_guard");}
unsigned inv(unsigned a){unsigned v=1;for(unsigned n=P-2;n;n>>=1,a=mul(a,a))if(n&1)v=mul(v,a);return v;}
void trim(V&v){while(!v.empty()&&!v.back())v.pop_back();}
void load(const char*path,const std::string&magic,int count,int degree){own(degree+1);R.resize(degree+1);std::ifstream f(path);std::string m,D,c;int n;f>>m>>D>>n;if(m!=magic||n!=count||D!="164144011096127878172550148707469312481250000000000000000000")throw std::runtime_error("header");
 for(int i=0;i<count;++i){long long k;f>>k>>c;if(!f)throw std::runtime_error("body");unsigned v=0;bool neg=c[0]=='-';for(size_t j=neg;j<c.size();++j)v=(10*v+c[j]-'0')%P;int a=k>>20,e=(k&1048575)>>10,d=k&1023;int w=646*a+19*e+d;if(w>degree)throw std::runtime_error("degree");R[w]=neg&&v?P-v:v;}
 if(!R.back())throw std::runtime_error("degree_dropped");
}
void divide(V&dest){int hd=H.back().first;if(H.back().second!=1)throw std::runtime_error("H_not_monic");own(R.size()-hd);T.assign(R.size()-hd,0);
 for(cur=(int)R.size()-1;cur>=hd;--cur){unsigned c=R[cur];T[cur-hd]=c;if(c)for(auto[h,v]:H)R[cur-hd+h]=(R[cur-hd+h]+P-mul(c,v))%P;if((cur&1023)==0)tick();}
 for(auto x:R)if(x)throw std::runtime_error("nonzero_division_remainder");dest.swap(T);V{}.swap(R);V{}.swap(T);own();
}
int main(int ac,char**av){if(ac==2&&std::string(av[1])=="--help"){std::cout<<"--run N0 N2 H OUT; exact sparse monic H division twice then coprime quotient gcd; 30k owned/150M products/45soft60hard/2GB.\n";return 0;}if(ac!=6)return 2;rlimit ram{2000000000,2000000000},cpu{60,60};setrlimit(RLIMIT_AS,&ram);setrlimit(RLIMIT_CPU,&cpu);begun=std::chrono::steady_clock::now();
 try{
  {std::ifstream f(av[3]);int degree,k;unsigned v;f>>degree;if(degree!=3899)throw std::runtime_error("H_degree");while(f>>k>>v){if(!v||v>=P||(!H.empty()&&k<=H.back().first))throw std::runtime_error("H_order");H.emplace_back(k,v);}if(H.size()!=94||H.back().first!=3899||H.back().second!=1)throw std::runtime_error("H_terms");own();}
  phase=1;load(av[1],"T0C4QV001",4500,10544);divide(A);if(A.size()!=6646)throw std::runtime_error("Q0_degree");phase=2;load(av[2],"Q1B4QV001",2317,7975);divide(B);if(B.size()!=4077)throw std::runtime_error("Q2_degree");phase=3;
  while(!B.empty()){unsigned il=inv(B.back());while(A.size()>=B.size()&&!A.empty()){cur=A.size()-B.size();unsigned c=mul(A.back(),il);for(size_t j=0;j<B.size();++j)A[cur+j]=(A[cur+j]+P-mul(c,B[j]))%P;trim(A);if((ops&8191)<B.size())tick();}A.swap(B);++pass;tick();}
  if(A.size()!=1||!A[0])throw std::runtime_error("quotients_not_coprime");std::ofstream r(std::string(av[4])+"/native-report.json");r<<"{\"all_passed\":true,\"both_original_remainders_zero\":true,\"quotient_degrees\":[6645,4076],\"quotients_coprime\":true,\"given_monic_H_is_exact_modular_gcd\":true,\"H_degree\":3899,\"H_terms\":94,\"products\":"<<ops<<",\"quotient_Euclid_passes\":"<<pass<<"}\n";return 0;
 }catch(const std::exception&e){std::ofstream f(std::string(av[4])+"/prefix.txt");f<<"GCDQUOTIENTPREFIX1\n"<<phase<<' '<<cur<<' '<<ops<<' '<<pass<<'\n';for(V*v:{&A,&B,&R,&T}){f<<v->size()<<'\n';for(auto x:*v)f<<x<<' ';f<<'\n';}std::ofstream r(std::string(av[4])+"/native-report.json");r<<"{\"all_passed\":false,\"prefix_saved\":true,\"reason\":\""<<e.what()<<"\"}\n";return 10;}
}
