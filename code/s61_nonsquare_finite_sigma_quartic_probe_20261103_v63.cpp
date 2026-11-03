// One bounded finite-Sigma quartic rank probe over the existing GF(32003).
#include <flint/nmod_mat.h>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
using namespace std;
constexpr uint32_t p=32003;constexpr int d=28;
using Small=array<uint32_t,d*d>;using Param=array<int,3>;using Mon=array<int,6>;
Small identity(){Small a{};for(int i=0;i<d;i++)a[i*d+i]=1;return a;}
Small product(const Small&a,const Small&b){Small c{};for(int i=0;i<d;i++)for(int j=0;j<d;j++){uint64_t z=0;for(int k=0;k<d;k++)z+=(uint64_t)a[i*d+k]*b[k*d+j];c[i*d+j]=z%p;}return c;}
void addscaled(Small&out,const Small&in,uint32_t c){for(int k=0;k<d*d;k++)out[k]=(out[k]+(uint64_t)c*in[k])%p;}
uint32_t power(uint32_t a,uint32_t n){uint32_t b=1;while(n){if(n&1)b=(uint64_t)b*a%p;a=(uint64_t)a*a%p;n>>=1;}return b;}
void monomials(vector<Mon>&out,Mon m,int i,int degree){if(i==5){m[i]=degree;out.push_back(m);return;}for(int k=0;k<=degree;k++){m[i]=k;monomials(out,m,i+1,degree-k);}}
void u32(ostream&f,uint32_t x){char b[4];for(int k=0;k<4;k++)b[k]=(char)((x>>(8*k))&255);f.write(b,4);}
int main(int argc,char**argv){
 assert(argc==5);auto start=chrono::steady_clock::now();ifstream ops(argv[1]),raw(argv[2]);int D,prime;ops>>D>>prime;assert(D==d&&prime==(int)p);
 Small MP{},ME{};for(auto*A:{&MP,&ME})for(auto&v:*A){ops>>v;assert(v<p);}assert(product(MP,ME)==product(ME,MP));
 int nr,nc,nv;raw>>nr>>nc>>nv;assert(nr==24&&nc==56&&nv==3);
 using Poly=map<Param,uint32_t>;vector<Poly>literal(nr*nc);int ma=0,me=0,mp=0;long terms=0;
 for(auto&q:literal){int n;raw>>n;assert(n>=0);terms+=n;for(int k=0;k<n;k++){uint32_t c;Param ex;raw>>c>>ex[0]>>ex[1]>>ex[2];assert(c<p&&ex[0]>=0&&ex[1]>=0&&ex[2]>=0);q[ex]=(q[ex]+c)%p;ma=max(ma,ex[0]);me=max(me,ex[1]);mp=max(mp,ex[2]);}for(auto it=q.begin();it!=q.end();)if(!it->second)it=q.erase(it);else++it;}
 assert(terms==4824);vector<Mon>m2,m3,m4;monomials(m2,Mon{},0,2);monomials(m3,Mon{},0,3);monomials(m4,Mon{},0,4);assert(m2.size()==21&&m3.size()==56&&m4.size()==126);
 map<Mon,int>idx2,idx3,idx4;for(int i=0;i<21;i++)idx2[m2[i]]=i;for(int i=0;i<56;i++)idx3[m3[i]]=i;for(int i=0;i<126;i++)idx4[m4[i]]=i;
 vector<Poly>quads(4*21);for(int i=0;i<4;i++)for(int j=0;j<56;j++){Mon m=m3[j];const auto&q=literal[(6*i)*56+j];if(m[0]){m[0]--;quads[i*21+idx2.at(m)]=q;}else assert(q.empty());}
 // Independently check all24 old rows before using the extracted four quadrics.
 for(int i=0;i<4;i++)for(int k=0;k<6;k++){vector<Poly>rebuilt(56);for(int j=0;j<21;j++){Mon m=m2[j];m[k]++;rebuilt[idx3.at(m)]=quads[i*21+j];}for(int j=0;j<56;j++)assert(rebuilt[j]==literal[(6*i+k)*56+j]);}
 Small I=identity(),temp{},MA{};addscaled(temp,MP,20);addscaled(temp,ME,3);addscaled(temp,I,20);MA=product(MP,temp);for(auto&x:MA)x=(uint64_t)x*(p-power(180,p-2))%p;
 vector<Small>AP(ma+1),EP(me+1),PP(mp+1);AP[0]=EP[0]=PP[0]=I;for(int i=1;i<=ma;i++)AP[i]=product(AP[i-1],MA);for(int i=1;i<=me;i++)EP[i]=product(EP[i-1],ME);for(int i=1;i<=mp;i++)PP[i]=product(PP[i-1],MP);
 map<Param,Small>cache;vector<Small>Qops(4*21);for(int i=0;i<84;i++)for(auto&term:quads[i]){auto ex=term.first;if(!cache.count(ex))cache[ex]=product(product(AP[ex[0]],EP[ex[1]]),PP[ex[2]]);addscaled(Qops[i],cache.at(ex),term.second);}
 constexpr slong rows=84*d,cols=126*d;nmod_mat_t A;nmod_mat_init(A,rows,cols,p);
 for(int i=0;i<4;i++)for(int shift=0;shift<21;shift++)for(int qmon=0;qmon<21;qmon++){
  Mon m{};for(int k=0;k<6;k++)m[k]=m2[shift][k]+m2[qmon][k];int col=idx4.at(m);const auto&V=Qops[i*21+qmon];int row=i*21+shift;
  for(int r=0;r<d;r++)for(int c=0;c<d;c++)nmod_mat_entry(A,row*d+r,col*d+c)=V[r*d+c];
 }
 uint64_t fnv=1469598103934665603ULL;ofstream block(argv[3],ios::binary);u32(block,rows);u32(block,cols);u32(block,p);
 for(slong i=0;i<rows;i++)for(slong j=0;j<cols;j++){uint32_t v=nmod_mat_entry(A,i,j);u32(block,v);fnv^=v;fnv*=1099511628211ULL;}block.close();assert(block.good());
 double build=chrono::duration<double>(chrono::steady_clock::now()-start).count();cout<<"BUILT_ACTUAL_MODULE 2352 3528 32003 RAW24_RECONSTRUCTION_PASS SEC "<<build<<endl;
 auto rankstart=chrono::steady_clock::now();slong rank=nmod_mat_rank(A);double ranksec=chrono::duration<double>(chrono::steady_clock::now()-rankstart).count();
 ofstream result(argv[4]);result<<"{\"prime\":32003,\"finite_algebra_dimension\":28,\"raw_source_reconstruction_all24_rows\":true,\"original_quartic_shape\":[84,126],\"regular_representation_shape\":[2352,3528],\"rank\":"<<rank<<",\"expected_rank\":2184,\"full_matrix_entries\":"<<rows*cols<<",\"full_matrix_FNV64\":"<<fnv<<",\"build_seconds\":"<<build<<",\"rank_seconds\":"<<ranksec<<",\"single_nmod_mat_rank_call\":true}\n";
 cout<<"EXACT_GF_RANK "<<rank<<" EXPECTED 2184 SEC "<<ranksec<<endl;nmod_mat_clear(A);
}
