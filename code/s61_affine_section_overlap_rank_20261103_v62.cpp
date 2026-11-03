#include <flint/nmod_mat.h>
#include <flint/ulong_extras.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;
constexpr int P=32003,EQS=134;
int N;
struct Mat {nmod_mat_t m;Mat(int r,int c){nmod_mat_init(m,r,c,P);}~Mat(){nmod_mat_clear(m);}ulong& at(int i,int j){return nmod_mat_entry(m,i,j);}};
template<size_t D> using Mon=array<int,D>;
template<size_t D> int key(const Mon<D>&a){int k=0;for(int i=0;i<D;i++)k|=a[i]<<(5*i);return k;}
template<size_t D> bool grev(const Mon<D>&a,const Mon<D>&b){for(int v=N-1;v>=0;v--){int x=count(a.begin(),a.end(),v),y=count(b.begin(),b.end(),v);if(x!=y)return x<y;}return false;}

// Compute the complete quartic row rank through cubic reductions and their
// overlaps. All rows are multiples of the original 134 input quadrics.
// A caller-supplied quotient lower bound may terminate a full-rank search:
// attaining that bound proves equality, while a failure only gives a lower
// bound on row rank. No Groebner-basis or flatness assumption is used.
int main(int argc,char**argv){
 assert(argc==4);ifstream in(argv[1]);ofstream out(argv[2]);int lower=stoi(argv[3]);
 int eq,p;in>>N>>eq>>p;assert(2<=N&&N<=32&&eq==EQS&&p==P&&lower>=0);
 auto start=chrono::steady_clock::now();auto tick=[&](const char*s){cerr<<s<<" "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"s\n";};
 vector<Mon<2>> pairs;vector<Mon<3>> cubes;vector<Mon<4>> quads;
 for(int i=0;i<N;i++)for(int j=i;j<N;j++)pairs.push_back({i,j});
 for(int i=0;i<N;i++)for(int j=i;j<N;j++)for(int k=j;k<N;k++)cubes.push_back({i,j,k});
 for(int i=0;i<N;i++)for(int j=i;j<N;j++)for(int k=j;k<N;k++)for(int l=k;l<N;l++)quads.push_back({i,j,k,l});
 sort(cubes.begin(),cubes.end(),grev<3>);sort(quads.begin(),quads.end(),grev<4>);
 int n2=pairs.size(),n3=cubes.size(),n4=quads.size();vector<int> cm(1<<15,-1),qm(1<<20,-1);
 for(int i=0;i<n3;i++)cm[key(cubes[i])]=i;for(int i=0;i<n4;i++)qm[key(quads[i])]=i;
 vector<vector<ulong>> F(EQS,vector<ulong>(n2));for(auto&row:F)for(auto&v:row){in>>v;assert(in&&v<P);}string extra;assert(!(in>>extra));
 Mat red(EQS*N,n3);for(int r=0;r<EQS;r++)for(int j=0;j<N;j++)for(int k=0;k<n2;k++)if(F[r][k]){
  Mon<3>a{pairs[k][0],pairs[k][1],j};sort(a.begin(),a.end());red.at(r*N+j,cm[key(a)])=F[r][k];}
 int rank=nmod_mat_rref(red.m);vector<int> piv,free,pivotrow(n3,-1);
 for(int i=0;i<rank;i++){int j=0;while(!red.at(i,j))j++;piv.push_back(j);pivotrow[j]=i;}
 for(int j=0;j<n3;j++)if(pivotrow[j]<0)free.push_back(j);int nc=free.size();
 Mat R(rank,nc);for(int i=0;i<rank;i++)for(int j=0;j<nc;j++){R.at(i,j)=(P-red.at(i,free[j]))%P;if(free[j]<piv[i])assert(!R.at(i,j));}
 cerr<<"CUBIC rank="<<rank<<" free="<<nc<<"\n";tick("CUBIC_REDUCED");
 vector<int> canon(n4,-1),var(n4,-1),standard(n4,-1),stdmons;vector<vector<pair<int,int>>> alternatives(n4);
 for(int m=0;m<n4;m++){
  auto a=quads[m];vector<pair<int,int>> choices;
  for(int removed=0;removed<4;removed++){
   Mon<3>b;int h=0;for(int k=0;k<4;k++)if(k!=removed)b[h++]=a[k];int row=pivotrow[cm[key(b)]];
   if(row>=0&&find(choices.begin(),choices.end(),make_pair(row,a[removed]))==choices.end())choices.push_back({row,a[removed]});}
  if(choices.empty()){standard[m]=stdmons.size();stdmons.push_back(m);}else{canon[m]=choices[0].first;var[m]=choices[0].second;alternatives[m]=choices;}
 }
 int d=stdmons.size();assert(d>=lower&&d<12000);cerr<<"QUARTIC_RESIDUAL "<<d<<"\n";tick("QUARTIC_ENUMERATED");
 vector<vector<int>> product(N,vector<int>(nc));
 for(int v=0;v<N;v++)for(int j=0;j<nc;j++){auto c=cubes[free[j]];Mon<4>a{c[0],c[1],c[2],v};sort(a.begin(),a.end());product[v][j]=qm[key(a)];}
 vector<uint16_t> NF(size_t(n4)*d);vector<uint64_t> accum(d);
 for(int m=n4-1;m>=0;m--){auto dest=&NF[size_t(m)*d];if(standard[m]>=0){dest[standard[m]]=1;continue;}fill(accum.begin(),accum.end(),0);
  for(int j=0;j<nc;j++)if(ulong c=R.at(canon[m],j)){int k=product[var[m]][j];assert(k>m);auto src=&NF[size_t(k)*d];for(int h=0;h<d;h++)accum[h]+=c*src[h];}
  for(int h=0;h<d;h++)dest[h]=accum[h]%P;
 }tick("QUARTIC_NORMAL_FORMS");
 vector<vector<ulong>> echelon(d);vector<array<int,4>> selected;int brank=0;ulong pivotproduct=1;
 for(int m=0;m<n4&&brank<d-lower;m++)for(auto choice:alternatives[m]){
  int row=choice.first,v=choice.second;if(row==canon[m]&&v==var[m])continue;
  fill(accum.begin(),accum.end(),0);for(int j=0;j<nc;j++)if(ulong c=R.at(row,j)){auto src=&NF[size_t(product[v][j])*d];for(int h=0;h<d;h++)accum[h]+=c*src[h];}
  vector<ulong> value(d);for(int h=0;h<d;h++)value[h]=(NF[size_t(m)*d+h]+P-accum[h]%P)%P;
  for(int h=0;h<d;h++)if(ulong c=value[h]){if(echelon[h].empty()){
    ulong inverse=n_invmod(c,P);echelon[h].resize(d);for(int k=h;k<d;k++)echelon[h][k]=value[k]*inverse%P;
    selected.push_back({row,v,h,int(c)});pivotproduct=pivotproduct*c%P;brank++;break;
   }else for(int k=h;k<d;k++)value[k]=(value[k]+P-c*echelon[h][k]%P)%P;}
  if(brank==d-lower)break;
 }tick("OVERLAP_RANK");
 out<<"N "<<N<<" RANK3 "<<rank<<" FREE3 "<<nc<<" RESIDUAL4 "<<d<<" OVERLAP_RANK "<<brank<<" RANK4 "<<n4-d+brank<<" QUOTIENT4_BOUND "<<d-brank<<" PIVOT_PRODUCT "<<pivotproduct<<"\n";
 out<<"CUBIC_PIVOTS";for(int a:piv)out<<' '<<a;out<<"\n";
 for(auto a:selected)out<<a[0]<<' '<<a[1]<<' '<<a[2]<<' '<<a[3]<<'\n';out.flush();assert(out.good());
 flint_cleanup();
}
