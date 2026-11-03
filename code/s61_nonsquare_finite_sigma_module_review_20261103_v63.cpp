// Independent standard-library modular arithmetic: actual finite-Sigma module and both inverse products.
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>
using namespace std;
constexpr uint32_t p=32003;constexpr int d=28;
using Small=array<uint32_t,d*d>;using Exp=array<int,3>;
Small identity(){Small a{};for(int i=0;i<d;i++)a[i*d+i]=1;return a;}
Small product(const Small&a,const Small&b){Small c{};for(int i=0;i<d;i++)for(int j=0;j<d;j++){uint64_t z=0;for(int k=0;k<d;k++)z+=(uint64_t)a[i*d+k]*b[k*d+j];c[i*d+j]=z%p;}return c;}
void addscaled(Small&out,const Small&in,uint32_t c){for(int k=0;k<d*d;k++)out[k]=(out[k]+(uint64_t)c*in[k])%p;}
uint32_t power(uint32_t a,uint32_t n){uint32_t out=1;while(n){if(n&1)out=(uint64_t)out*a%p;a=(uint64_t)a*a%p;n>>=1;}return out;}
struct Term{uint32_t c;Exp ex;};
int main(int argc,char**argv){
 assert(argc==7);auto started=chrono::steady_clock::now();ifstream ops(argv[1]),raw(argv[2]),polys(argv[3]),choose(argv[4]),inverse(argv[5]);
 int D,prime;ops>>D>>prime;assert(D==d&&prime==(int)p);Small MP{},ME{};for(auto*A:{&MP,&ME})for(auto&v:*A){ops>>v;assert(v<p);}assert(product(MP,ME)==product(ME,MP));
 int nr,nc,nvars;raw>>nr>>nc>>nvars;assert(nr==24&&nc==56&&nvars==3);vector<vector<Term>>literal(nr*nc);set<Exp>keys;int ma=0,me=0,mp=0;long nterms=0;
 for(auto&poly:literal){int n;raw>>n;for(int j=0;j<n;j++){Term t;raw>>t.c>>t.ex[0]>>t.ex[1]>>t.ex[2];assert(t.c<p);assert(t.ex[0]>=0&&t.ex[1]>=0&&t.ex[2]>=0);ma=max(ma,t.ex[0]);me=max(me,t.ex[1]);mp=max(mp,t.ex[2]);keys.insert(t.ex);poly.push_back(t);nterms++;}}
 int np,nv;polys>>np>>nv;assert(np==4&&nv==2);vector<vector<Term>>models(4);
 for(auto&poly:models){int n;polys>>n;for(int j=0;j<n;j++){Term t;t.ex[0]=0;polys>>t.c>>t.ex[2]>>t.ex[1];assert(t.c<p);mp=max(mp,t.ex[2]);me=max(me,t.ex[1]);poly.push_back(t);}}
 Small I=identity(),aMat{},tmp{};addscaled(tmp,MP,20);addscaled(tmp,ME,3);addscaled(tmp,I,20);aMat=product(MP,tmp);for(auto&c:aMat)c=(uint64_t)c*(p-power(180,p-2))%p;
 vector<Small>AP(ma+1),EP(me+1),PP(mp+1);AP[0]=EP[0]=PP[0]=I;for(int i=1;i<=ma;i++)AP[i]=product(AP[i-1],aMat);for(int i=1;i<=me;i++)EP[i]=product(EP[i-1],ME);for(int i=1;i<=mp;i++)PP[i]=product(PP[i-1],MP);
 Small zero{};array<Small,4>modelops{};for(int i=0;i<4;i++)for(auto&t:models[i])addscaled(modelops[i],product(PP[t.ex[2]],EP[t.ex[1]]),t.c);
 assert(modelops[0]==zero&&modelops[1]==zero);assert(product(modelops[2],modelops[3])==I&&product(modelops[3],modelops[2])==I);
 map<Exp,Small>monomial;for(auto&ex:keys)monomial[ex]=product(product(PP[ex[2]],EP[ex[1]]),AP[ex[0]]);
 int rows=24*d,cols=56*d;vector<uint32_t>C(rows*cols);uint64_t checksum=1469598103934665603ULL;
 for(int blockcol=0;blockcol<56;blockcol++)for(int blockrow=0;blockrow<24;blockrow++){
  Small val{};for(auto&t:literal[blockrow*56+blockcol])addscaled(val,monomial.at(t.ex),t.c);
  for(int i=0;i<d;i++)for(int j=0;j<d;j++)C[(blockrow*d+i)*cols+blockcol*d+j]=val[i*d+j];
 }
 for(uint32_t c:C){checksum^=c;checksum*=1099511628211ULL;}
 int nsel;choose>>nsel;assert(nsel==rows);vector<int>selected(rows);set<int>distinct;for(auto&c:selected){choose>>c;assert(c>=0&&c<cols);distinct.insert(c);}assert((int)distinct.size()==rows);
 int nInv,invPrime;inverse>>nInv>>invPrime;assert(nInv==rows&&invPrime==(int)p);vector<uint32_t>S(rows*rows),Inv(rows*rows),St(rows*rows),It(rows*rows);
 for(int i=0;i<rows;i++)for(int j=0;j<rows;j++){S[i*rows+j]=C[i*cols+selected[j]];inverse>>Inv[i*rows+j];assert(Inv[i*rows+j]<p);}
 for(int i=0;i<rows;i++)for(int j=0;j<rows;j++){St[j*rows+i]=S[i*rows+j];It[j*rows+i]=Inv[i*rows+j];}
 // Independent integer-dot-product checks in both orders; no rank routine or field-matrix library.
 for(int i=0;i<rows;i++)for(int j=0;j<rows;j++){
  uint64_t v=0,u=0;for(int k=0;k<rows;k++){v+=(uint64_t)S[i*rows+k]*It[j*rows+k];u+=(uint64_t)Inv[i*rows+k]*St[j*rows+k];}
  assert(v%p==(uint32_t)(i==j));assert(u%p==(uint32_t)(i==j));
 }
 ofstream out(argv[6]);out<<"{\"all_passed\":true,\"prime\":32003,\"algebra_dimension\":28,\"literal_raw_rows\":24,\"literal_raw_columns\":56,\"literal_terms\":"<<nterms<<",\"full_regular_representation_entries\":"<<C.size()<<",\"full_module_shape\":[672,1568],\"regular_operators_commute\":true,\"N_M_operator_relations_zero\":true,\"U_Uinverse_both_products_identity\":true,\"selected_square_both_inverse_products_identity\":true,\"inverse_identity_entries_checked\":"<<2LL*rows*rows<<",\"full_matrix_FNV64\":"<<checksum<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"}\n";
}
