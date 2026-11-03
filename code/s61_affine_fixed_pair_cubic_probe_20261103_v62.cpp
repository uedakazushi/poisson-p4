#include <flint/nmod_mat.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
constexpr int N=44,P=32003,EQS=134;
using Cube=array<int,3>;
int key(Cube a){return a[0]|(a[1]<<6)|(a[2]<<12);}
bool grev(Cube a,Cube b){
 int x=count(a.begin(),a.end(),0),y=count(b.begin(),b.end(),0);if(x!=y)return x<y;
 for(int v=N-1;v>=1;v--){x=count(a.begin(),a.end(),v);y=count(b.begin(),b.end(),v);if(x!=y)return x<y;}return false;
}
int main(int argc,char**argv){
 assert(argc==4);ifstream in(argv[1]);ofstream out(argv[2]);ofstream dualout(argv[3]);int n,eq,p;in>>n>>eq>>p;assert(n==N&&eq==EQS&&p==P);
 auto start=chrono::steady_clock::now();vector<array<int,2>>pairs;vector<Cube>cubes;
 for(int i=0;i<N;i++)for(int j=i;j<N;j++)pairs.push_back({i,j});
 for(int i=0;i<N;i++)for(int j=i;j<N;j++)for(int k=j;k<N;k++)cubes.push_back({i,j,k});
 sort(cubes.begin(),cubes.end(),grev);int n3=cubes.size();vector<int>cm(1<<18,-1);for(int i=0;i<n3;i++)cm[key(cubes[i])]=i;
 nmod_mat_t M;nmod_mat_init(M,EQS*N,n3,P);vector<vector<ulong>>F(EQS,vector<ulong>(pairs.size()));
 for(int r=0;r<EQS;r++)for(int j=0;j<(int)pairs.size();j++){auto pair=pairs[j];ulong c;in>>c;assert(in&&c<P);F[r][j]=c;if(c)for(int k=0;k<N;k++){
  Cube a{pair[0],pair[1],k};sort(a.begin(),a.end());nmod_mat_entry(M,r*N+k,cm[key(a)])=c;}}
 string extra;assert(!(in>>extra));cerr<<"FULL_CUBIC_MATRIX_READY "<<EQS*N<<' '<<n3<<'\n';
 int rank=nmod_mat_rref(M);vector<int>pr(n3,-1),free;for(int r=0;r<rank;r++){int j=0;while(!nmod_mat_entry(M,r,j))j++;pr[j]=r;}
 for(int j=0;j<n3;j++)if(pr[j]<0)free.push_back(j);nmod_mat_t W;nmod_mat_init(W,N,free.size(),P);
 for(int k=0;k<N;k++){
  Cube a{0,0,k};sort(a.begin(),a.end());int j=cm[key(a)];
  if(pr[j]>=0){for(int h=0;h<(int)free.size();h++)nmod_mat_entry(W,k,h)=(P-nmod_mat_entry(M,pr[j],free[h]))%P;}
  else{auto it=find(free.begin(),free.end(),j);assert(it!=free.end());nmod_mat_entry(W,k,it-free.begin())=1;}
 }
 bool t3zero=true;for(int h=0;h<(int)free.size();h++)if(nmod_mat_entry(W,0,h))t3zero=false;
 int image=nmod_mat_rank(W);assert(image==N);
 nmod_mat_t WR,B,BI,D,left,right;nmod_mat_init(WR,N,free.size(),P);nmod_mat_set(WR,W);assert(nmod_mat_rref(WR)==N);
 vector<int>chosen;for(int i=0;i<N;i++){int j=0;while(!nmod_mat_entry(WR,i,j))j++;chosen.push_back(j);}
 nmod_mat_init(B,N,N,P);nmod_mat_init(BI,N,N,P);for(int i=0;i<N;i++)for(int j=0;j<N;j++)nmod_mat_entry(B,i,j)=nmod_mat_entry(W,i,chosen[j]);assert(nmod_mat_inv(BI,B));
 nmod_mat_init(D,n3,N,P);for(int j=0;j<N;j++)for(int l=0;l<N;l++)nmod_mat_entry(D,free[chosen[j]],l)=nmod_mat_entry(BI,j,l);
 nmod_mat_init(left,rank,N,P);nmod_mat_init(right,rank,N,P);for(int j=0;j<n3;j++)if(pr[j]>=0)for(int k=0;k<N;k++)nmod_mat_entry(left,pr[j],k)=nmod_mat_entry(M,pr[j],free[chosen[k]]);
 nmod_mat_mul(right,left,BI);for(int j=0;j<n3;j++)if(pr[j]>=0)for(int l=0;l<N;l++)nmod_mat_entry(D,j,l)=(P-nmod_mat_entry(right,pr[j],l))%P;
 for(int k=0;k<N;k++){Cube a{0,0,k};sort(a.begin(),a.end());for(int l=0;l<N;l++)assert(nmod_mat_entry(D,cm[key(a)],l)==ulong(k==l));}
 for(int r=0;r<EQS;r++)for(int k=0;k<N;k++){
  vector<uint64_t>acc(N);for(int j=0;j<(int)pairs.size();j++)if(ulong c=F[r][j]){Cube a{pairs[j][0],pairs[j][1],k};sort(a.begin(),a.end());int m=cm[key(a)];for(int l=0;l<N;l++)acc[l]+=c*nmod_mat_entry(D,m,l);}
  for(auto c:acc)assert(c%P==0);
 }
 dualout<<"N "<<N<<" DUAL_COLUMNS "<<N<<" PRIME "<<P<<"\n";
 for(int j=0;j<n3;j++){bool nonzero=false;for(int l=0;l<N;l++)if(nmod_mat_entry(D,j,l))nonzero=true;if(nonzero){for(int a:cubes[j])dualout<<a<<' ';for(int l=0;l<N;l++)dualout<<nmod_mat_entry(D,j,l)<<(l==N-1?'\n':' ');}}
 dualout.flush();assert(dualout.good());
 out<<"N "<<N<<" ORIGINAL134_CUBIC_RANK "<<rank<<" CUBIC_FREE "<<free.size()<<" T2_LINEAR_IMAGE_RANK "<<image<<" T3_ZERO "<<t3zero<<" DUAL_IDENTITIES_VERIFIED 1\n";
 cerr<<"FIXED_PAIR_CUBIC_PROBE_COMPLETE "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"s\n";out.flush();assert(out.good());
}
