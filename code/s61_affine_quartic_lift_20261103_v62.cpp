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
#ifndef S61_LIFT_N
#define S61_LIFT_N 28
#endif
constexpr int P=32003,N=S61_LIFT_N,EQS=134;
struct Mat {nmod_mat_t m;Mat(int r,int c){nmod_mat_init(m,r,c,P);}~Mat(){nmod_mat_clear(m);}ulong& at(int i,int j){return nmod_mat_entry(m,i,j);}};
template<size_t D> using Mon=array<int,D>;
template<size_t D> int key(const Mon<D>& a){int k=0;for(int i=0;i<D;i++)k|=a[i]<<(5*i);return k;}
template<size_t D> bool grev(const Mon<D>&a,const Mon<D>&b){for(int v=N-1;v>=0;v--){int x=count(a.begin(),a.end(),v),y=count(b.begin(),b.end(),v);if(x!=y)return x<y;}return false;}
vector<int> selectrows(Mat& a,int rank){Mat tmp(a.m->r,a.m->c);nmod_mat_set(tmp.m,a.m);vector<slong> perm(a.m->r);iota(perm.begin(),perm.end(),0);assert(nmod_mat_lu(perm.data(),tmp.m,0)==rank);return {perm.begin(),perm.begin()+rank};}
int main(int argc,char**argv){
 assert(argc==4);ifstream in(argv[1]);ofstream out(argv[2]);int targetvar=stoi(argv[3]);assert(targetvar>=0&&targetvar<N);
 int nn,eq,p;in>>nn>>eq>>p;assert(nn==N&&eq==EQS&&p==P);auto start=chrono::steady_clock::now();
 auto tick=[&](const char* label){cerr<<label<<" "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"s\n";};
 vector<Mon<2>> pairs;vector<Mon<3>> cubes;vector<Mon<4>> quads;
 int pairmap[N][N];
 for(int i=0;i<N;i++)for(int j=i;j<N;j++){pairmap[i][j]=pairmap[j][i]=pairs.size();pairs.push_back({i,j});}
 for(int i=0;i<N;i++)for(int j=i;j<N;j++)for(int k=j;k<N;k++)cubes.push_back({i,j,k});
 for(int i=0;i<N;i++)for(int j=i;j<N;j++)for(int k=j;k<N;k++)for(int l=k;l<N;l++)quads.push_back({i,j,k,l});
 sort(cubes.begin(),cubes.end(),grev<3>);sort(quads.begin(),quads.end(),grev<4>);
 const int n2=pairs.size(),n3=cubes.size(),n4=quads.size();vector<int> cm(1<<15,-1),qm(1<<20,-1);
 for(int i=0;i<n3;i++)cm[key(cubes[i])]=i;for(int i=0;i<n4;i++)qm[key(quads[i])]=i;
 vector<vector<ulong>> F(EQS,vector<ulong>(n2));for(auto&row:F)for(auto&v:row){in>>v;assert(in&&v<P);}string extra;assert(!(in>>extra));
 Mat M3(EQS*N,n3);for(int r=0;r<EQS;r++)for(int j=0;j<N;j++)for(int k=0;k<n2;k++)if(F[r][k]){
  Mon<3> a{pairs[k][0],pairs[k][1],j};sort(a.begin(),a.end());M3.at(r*N+j,cm[key(a)])=F[r][k];}
 vector<int> piv,free;int rank;
 {Mat red(EQS*N,n3);nmod_mat_set(red.m,M3.m);rank=nmod_mat_rref(red.m);for(int i=0;i<rank;i++){int j=0;while(!red.at(i,j))j++;piv.push_back(j);}}
 assert(rank>0);vector<int> pivotrow(n3,-1);for(int i=0;i<rank;i++)pivotrow[piv[i]]=i;
 for(int j=0;j<n3;j++)if(pivotrow[j]<0)free.push_back(j);const int nc=free.size();assert(nc>0 && nc<1000);tick("CUBIC_RANK");
 auto rows=selectrows(M3,rank);Mat inv(rank,rank);ulong det3;
 {Mat square(rank,rank);for(int i=0;i<rank;i++)for(int j=0;j<rank;j++)square.at(i,j)=M3.at(rows[i],piv[j]);det3=nmod_mat_det(square.m);assert(det3);assert(nmod_mat_inv(inv.m,square.m));}
 Mat tail(rank,nc),R(rank,nc);for(int i=0;i<rank;i++)for(int j=0;j<nc;j++)tail.at(i,j)=M3.at(rows[i],free[j]);nmod_mat_mul(R.m,inv.m,tail.m);nmod_mat_neg(R.m,R.m);
 for(int i=0;i<rank;i++)for(int j=0;j<nc;j++)if(free[j]<piv[i])assert(!R.at(i,j));tick("CUBIC_REDUCTION");
 vector<int> canon(n4,-1),var(n4,-1),standard(n4,-1),stdmons;vector<vector<pair<int,int>>> alternatives(n4);
 for(int m=0;m<n4;m++){
  auto a=quads[m];vector<pair<int,int>> choices;
  for(int removed=0;removed<4;removed++){
   Mon<3>b;int h=0;for(int k=0;k<4;k++)if(k!=removed)b[h++]=a[k];int row=pivotrow[cm[key(b)]];
   if(row>=0&&find(choices.begin(),choices.end(),make_pair(row,a[removed]))==choices.end())choices.push_back({row,a[removed]});}
  if(choices.empty()){standard[m]=stdmons.size();stdmons.push_back(m);}else{canon[m]=choices[0].first;var[m]=choices[0].second;alternatives[m]=choices;}
 }
 int d=stdmons.size();assert(d>0 && d<2000);vector<vector<int>> product(N,vector<int>(nc));
 for(int v=0;v<N;v++)for(int j=0;j<nc;j++){auto c=cubes[free[j]];Mon<4>a{c[0],c[1],c[2],v};sort(a.begin(),a.end());product[v][j]=qm[key(a)];}
 vector<uint16_t> NF(size_t(n4)*d);vector<uint64_t> accum(d);
 for(int m=n4-1;m>=0;m--){auto dest=&NF[size_t(m)*d];if(standard[m]>=0){dest[standard[m]]=1;continue;}fill(accum.begin(),accum.end(),0);
  for(int j=0;j<nc;j++)if(ulong c=R.at(canon[m],j)){int k=product[var[m]][j];assert(k>m);auto src=&NF[size_t(k)*d];for(int h=0;h<d;h++)accum[h]+=c*src[h];}
  for(int h=0;h<d;h++)dest[h]=accum[h]%P;
 }tick("QUARTIC_NORMAL_FORMS");
 vector<vector<ulong>> echelon(d);vector<pair<int,int>> selected;Mat B(d,d);int brank=0;
 for(int m=0;m<n4&&brank<d;m++)for(auto choice:alternatives[m]){
  int row=choice.first,v=choice.second;if(row==canon[m]&&v==var[m])continue;
  fill(accum.begin(),accum.end(),0);for(int j=0;j<nc;j++)if(ulong c=R.at(row,j)){auto src=&NF[size_t(product[v][j])*d];for(int h=0;h<d;h++)accum[h]+=c*src[h];}
  vector<ulong> raw(d),value(d);for(int h=0;h<d;h++)raw[h]=value[h]=(NF[size_t(m)*d+h]+P-accum[h]%P)%P;
  for(int h=0;h<d;h++)if(ulong c=value[h]){if(echelon[h].empty()){
    ulong inverse=n_invmod(c,P);echelon[h].resize(d);for(int k=h;k<d;k++)echelon[h][k]=value[k]*inverse%P;
    selected.push_back(choice);for(int k=0;k<d;k++)B.at(brank,k)=raw[k];brank++;break;
   }else for(int k=h;k<d;k++)value[k]=(value[k]+P-c*echelon[h][k]%P)%P;}
  if(brank==d)break;
 }
 assert(brank==d);ulong det4=nmod_mat_det(B.m);assert(det4);Mat Binv(d,d);assert(nmod_mat_inv(Binv.m,B.m));tick("QUARTIC_FULL_RANK");
 Mon<4> tar{targetvar,targetvar,targetvar,targetvar};int target=qm[key(tar)];Mat tnf(1,d),lambda(1,d);for(int k=0;k<d;k++)tnf.at(0,k)=NF[size_t(target)*d+k];nmod_mat_mul(lambda.m,tnf.m,Binv.m);
 vector<ulong> poly(n4);poly[target]=1;Mat C(N,rank);
 for(int h=0;h<d;h++)if(ulong c=lambda.at(0,h)){
  int row=selected[h].first,v=selected[h].second;auto a=cubes[piv[row]];Mon<4>b{a[0],a[1],a[2],v};sort(b.begin(),b.end());int m=qm[key(b)];poly[m]=(poly[m]+P-c)%P;
  for(int j=0;j<nc;j++)poly[product[v][j]]=(poly[product[v][j]]+c*R.at(row,j))%P;
  C.at(v,row)=(C.at(v,row)+c)%P;
 }
 for(int m=0;m<n4;m++)if(ulong c=poly[m]){
  assert(canon[m]>=0);int row=canon[m],v=var[m];poly[m]=0;
  for(int j=0;j<nc;j++)poly[product[v][j]]=(poly[product[v][j]]+c*R.at(row,j))%P;
  C.at(v,row)=(C.at(v,row)+c)%P;
 }
 assert(all_of(poly.begin(),poly.end(),[](ulong c){return c==0;}));
 Mat original(N,rank);nmod_mat_mul(original.m,C.m,inv.m);vector<vector<ulong>> lift(EQS,vector<ulong>(n2));
 for(int v=0;v<N;v++)for(int i=0;i<rank;i++){int r=rows[i]/N,k=rows[i]%N;int j=pairmap[v][k];lift[r][j]=(lift[r][j]+original.at(v,i))%P;}
 vector<uint64_t> verification(n4);int terms=0;
 for(int r=0;r<EQS;r++)for(int j=0;j<n2;j++)if(ulong c=lift[r][j]){terms++;for(int k=0;k<n2;k++)if(F[r][k]){
  Mon<4>a{pairs[j][0],pairs[j][1],pairs[k][0],pairs[k][1]};sort(a.begin(),a.end());verification[qm[key(a)]]+=c*F[r][k];}}
 for(int m=0;m<n4;m++)assert(verification[m]%P==ulong(m==target));tick("ORIGINAL134_IDENTITY_VERIFIED");
 out<<"RANK3 "<<rank<<" FREE3 "<<nc<<" DET3 "<<det3<<" RANK4 "<<n4<<" DET4 "<<det4<<" TARGET "<<targetvar<<" TERMS "<<terms<<"\n";
 out<<"CUBIC_ROWS";for(int r:rows)out<<' '<<r;out<<"\nCUBIC_PIVOTS";for(int r:piv)out<<' '<<r;out<<"\nQUARTIC_ALTERNATIVES";for(auto q:selected)out<<' '<<q.first<<' '<<q.second;out<<'\n';
 for(int r=0;r<EQS;r++)for(int j=0;j<n2;j++)if(lift[r][j])out<<r<<' '<<pairs[j][0]<<' '<<pairs[j][1]<<' '<<lift[r][j]<<'\n';out.flush();assert(out.good());
 flint_cleanup();
}
