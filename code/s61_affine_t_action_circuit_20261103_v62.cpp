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
#include <random>
#include <vector>
using namespace std;
#ifndef S61_LIFT_N
#define S61_LIFT_N 29
#endif
constexpr int P=32003,N=S61_LIFT_N,EQS=134;
struct Mat {nmod_mat_t m;Mat(int r,int c){nmod_mat_init(m,r,c,P);}~Mat(){nmod_mat_clear(m);}ulong& at(int i,int j){return nmod_mat_entry(m,i,j);}};
template<size_t D> using Mon=array<int,D>;
template<size_t D> int key(const Mon<D>& a){int k=0;for(int i=0;i<D;i++)k|=a[i]<<(5*i);return k;}
template<size_t D> bool grev(const Mon<D>&a,const Mon<D>&b){for(int v=N-1;v>=0;v--){int x=count(a.begin(),a.end(),v),y=count(b.begin(),b.end(),v);if(x!=y)return x<y;}return false;}
vector<int> selectrows(Mat& a,int rank){Mat tmp(a.m->r,a.m->c);nmod_mat_set(tmp.m,a.m);vector<slong> perm(a.m->r);iota(perm.begin(),perm.end(),0);assert(nmod_mat_lu(perm.data(),tmp.m,0)==rank);return {perm.begin(),perm.begin()+rank};}
int main(int argc,char**argv){
 assert(argc==4||argc==5);ifstream in(argv[1]);ofstream out(argv[2]);int trials=stoi(argv[3]);assert(trials>=3&&trials<=30);
 int nn,eq,p,nb;in>>nn>>eq>>p>>nb;assert(nn==N&&N==29&&eq==EQS&&p==P&&nb==17);auto start=chrono::steady_clock::now();
 auto tick=[&](const char* label){cerr<<label<<" "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"s\n";};
 vector<Mon<2>> pairs;vector<Mon<3>> cubes;vector<Mon<4>> quads;
 int pairmap[N][N];
 for(int i=0;i<N;i++)for(int j=i;j<N;j++){pairmap[i][j]=pairmap[j][i]=pairs.size();pairs.push_back({i,j});}
 for(int i=0;i<N;i++)for(int j=i;j<N;j++)for(int k=j;k<N;k++)cubes.push_back({i,j,k});
 for(int i=0;i<N;i++)for(int j=i;j<N;j++)for(int k=j;k<N;k++)for(int l=k;l<N;l++)quads.push_back({i,j,k,l});
 sort(cubes.begin(),cubes.end(),grev<3>);sort(quads.begin(),quads.end(),grev<4>);
 const int n2=pairs.size(),n3=cubes.size(),n4=quads.size();vector<int> cm(1<<15,-1),qm(1<<20,-1);
 for(int i=0;i<n3;i++)cm[key(cubes[i])]=i;for(int i=0;i<n4;i++)qm[key(quads[i])]=i;
 vector<vector<ulong>> F(EQS,vector<ulong>(n2));for(auto&row:F)for(auto&v:row){in>>v;assert(in&&v<P);}vector<vector<ulong>> cross(EQS,vector<ulong>(N*nb)),base(EQS,vector<ulong>(nb*(nb+1)/2));
 for(auto&row:cross)for(auto&v:row){in>>v;assert(in&&v<P);}for(auto&row:base)for(auto&v:row){in>>v;assert(in&&v<P);}
 vector<ulong> knownbase(nb);for(auto&v:knownbase){in>>v;assert(in&&v<P);}string extra;assert(!(in>>extra));
 Mat quad(EQS,n2);for(int r=0;r<EQS;r++)for(int j=0;j<n2;j++)quad.at(r,j)=F[r][j];assert(nmod_mat_rref(quad.m)==EQS);
 vector<int> piv2,free2,pr2(n2,-1);for(int r=0;r<EQS;r++){int j=0;while(!quad.at(r,j))j++;piv2.push_back(j);pr2[j]=r;}for(int j=0;j<n2;j++)if(pr2[j]<0)free2.push_back(j);
 Mat square2(EQS,EQS),inv2(EQS,EQS),R2(EQS,free2.size());for(int r=0;r<EQS;r++)for(int j=0;j<EQS;j++)square2.at(r,j)=F[r][piv2[j]];assert(nmod_mat_inv(inv2.m,square2.m));
 for(int r=0;r<EQS;r++)for(int j=0;j<(int)free2.size();j++)R2.at(r,j)=(P-quad.at(r,free2[j]))%P;
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
 const int g=1+N+free2.size()+nc;
 assert(g==1031);
 // Optional, reproducible scratch cache for subsequent symbolic reductions.
 // The compact certificate does not require committing these dense matrices.
 if(argc==5){
  ofstream cache(argv[4],ios::binary);uint32_t header[]={0x53363154,N,nb,EQS,n2,n3,rank,nc};cache.write((char*)header,sizeof(header));
  auto save=[&](Mat&a){for(int i=0;i<a.m->r;i++)for(int j=0;j<a.m->c;j++){uint16_t c=a.at(i,j);cache.write((char*)&c,sizeof(c));}};
  save(inv2);save(R2);save(inv);save(R);cache.flush();assert(cache.good());
 }
 vector<int> free2index(n2,-1),free3index(n3,-1);for(int j=0;j<(int)free2.size();j++)free2index[free2[j]]=j;for(int j=0;j<nc;j++)free3index[free[j]]=j;
 auto add=[](ulong&a,ulong b){a=(a+b)%P;};
 auto sub=[](ulong&a,ulong b){a=(a+P-b%P)%P;};
 auto makeinput=[&](const vector<ulong>&u){
  vector<vector<ulong>> value(5);value[0].resize(1);value[1].resize(N);value[2].resize(n2);value[3].resize(n3);value[4].resize(n4);
  value[1][0]=u[0];
  for(int i=0;i<N;i++)value[2][pairmap[0][i]]=u[1+i];
  for(int j=0;j<(int)free2.size();j++){auto b=pairs[free2[j]];Mon<3>a{0,b[0],b[1]};sort(a.begin(),a.end());value[3][cm[key(a)]]=u[1+N+j];}
  for(int j=0;j<nc;j++){auto b=cubes[free[j]];Mon<4>a{0,b[0],b[1],b[2]};sort(a.begin(),a.end());value[4][qm[key(a)]]=u[1+N+free2.size()+j];}
  return value;
 };
 // This is a linear operator in u, with polynomial dependence on all17 Y.
 // Every reduction subtracts a recorded original-equation multiple, so its
 // correctness is a formal identity, not an inference from sampled Y.
 auto apply=[&](const vector<ulong>&u,const vector<ulong>&Y,int trial){
  auto input=makeinput(u);auto value=input;
  vector<vector<ulong>> lin(EQS,vector<ulong>(N));vector<ulong> constant(EQS),ypairs;
  for(int i=0;i<nb;i++)for(int j=i;j<nb;j++)ypairs.push_back(Y[i]*Y[j]%P);
  for(int r=0;r<EQS;r++){
   for(int i=0;i<N;i++)for(int b=0;b<nb;b++)add(lin[r][i],cross[r][i*nb+b]*Y[b]);
   for(int b=0;b<(int)ypairs.size();b++)add(constant[r],base[r][b]*ypairs[b]);
  }
  Mat tnf(1,d),lambda(1,d);
  for(int m=0;m<n4;m++)if(ulong c=value[4][m])for(int h=0;h<d;h++)tnf.at(0,h)=(tnf.at(0,h)+c*NF[size_t(m)*d+h])%P;
  nmod_mat_mul(lambda.m,tnf.m,Binv.m);Mat C(N,rank);
  for(int h=0;h<d;h++)if(ulong c=lambda.at(0,h)){
   int row=selected[h].first,v=selected[h].second;auto a=cubes[piv[row]];Mon<4>b{a[0],a[1],a[2],v};sort(b.begin(),b.end());int m=qm[key(b)];sub(value[4][m],c);
   for(int j=0;j<nc;j++)add(value[4][product[v][j]],c*R.at(row,j));add(C.at(v,row),c);
  }
  for(int m=0;m<n4;m++)if(ulong c=value[4][m]){
   assert(canon[m]>=0);int row=canon[m],v=var[m];value[4][m]=0;
   for(int j=0;j<nc;j++)add(value[4][product[v][j]],c*R.at(row,j));add(C.at(v,row),c);
  }
  assert(all_of(value[4].begin(),value[4].end(),[](ulong c){return c==0;}));
  Mat original(N,rank);nmod_mat_mul(original.m,C.m,inv.m);vector<vector<ulong>> H(EQS,vector<ulong>(n2));
  for(int v=0;v<N;v++)for(int i=0;i<rank;i++){int r=rows[i]/N,k=rows[i]%N;add(H[r][pairmap[v][k]],original.at(v,i));}
  for(int r=0;r<EQS;r++)for(int j=0;j<n2;j++)if(ulong c=H[r][j]){
   sub(value[2][j],c*constant[r]);
   for(int i=0;i<N;i++){Mon<3>a{pairs[j][0],pairs[j][1],i};sort(a.begin(),a.end());sub(value[3][cm[key(a)]],c*lin[r][i]);}
  }
  Mat C3(1,rank),original3(1,rank);for(int i=0;i<rank;i++){
   ulong c=value[3][piv[i]];C3.at(0,i)=c;value[3][piv[i]]=0;for(int j=0;j<nc;j++)add(value[3][free[j]],c*R.at(i,j));}
  nmod_mat_mul(original3.m,C3.m,inv.m);vector<vector<ulong>> L(EQS,vector<ulong>(N));
  for(int i=0;i<rank;i++)add(L[rows[i]/N][rows[i]%N],original3.at(0,i));
  for(int r=0;r<EQS;r++)for(int k=0;k<N;k++)if(ulong c=L[r][k]){
   sub(value[1][k],c*constant[r]);for(int i=0;i<N;i++)sub(value[2][pairmap[k][i]],c*lin[r][i]);}
  Mat C2(1,EQS),original2(1,EQS);for(int i=0;i<EQS;i++){
   ulong c=value[2][piv2[i]];C2.at(0,i)=c;value[2][piv2[i]]=0;for(int j=0;j<(int)free2.size();j++)add(value[2][free2[j]],c*R2.at(i,j));}
  nmod_mat_mul(original2.m,C2.m,inv2.m);vector<ulong>D(EQS);for(int r=0;r<EQS;r++){
   D[r]=original2.at(0,r);sub(value[0][0],D[r]*constant[r]);for(int i=0;i<N;i++)sub(value[1][i],D[r]*lin[r][i]);}
  vector<ulong> result(g);result[0]=value[0][0];for(int i=0;i<N;i++)result[1+i]=value[1][i];
  for(int j=0;j<(int)free2.size();j++)result[1+N+j]=value[2][free2[j]];
  for(int j=0;j<nc;j++)result[1+N+free2.size()+j]=value[3][free[j]];
  // Independently re-expand all original134 equations, checking every
  // unknown-monomial coefficient in degrees0,...,4 for this call.
  vector<vector<uint64_t>> check(5);check[0].resize(1);check[1].resize(N);check[2].resize(n2);check[3].resize(n3);check[4].resize(n4);
  for(int r=0;r<EQS;r++){
   check[0][0]+=D[r]*constant[r];for(int i=0;i<N;i++)check[1][i]+=D[r]*lin[r][i];for(int j=0;j<n2;j++)check[2][j]+=D[r]*F[r][j];
   for(int k=0;k<N;k++)if(ulong c=L[r][k]){
    check[1][k]+=c*constant[r];for(int i=0;i<N;i++)check[2][pairmap[k][i]]+=c*lin[r][i];
    for(int j=0;j<n2;j++){Mon<3>a{k,pairs[j][0],pairs[j][1]};sort(a.begin(),a.end());check[3][cm[key(a)]]+=c*F[r][j];}}
   for(int j=0;j<n2;j++)if(ulong c=H[r][j]){
    check[2][j]+=c*constant[r];for(int i=0;i<N;i++){Mon<3>a{pairs[j][0],pairs[j][1],i};sort(a.begin(),a.end());check[3][cm[key(a)]]+=c*lin[r][i];}
    for(int k=0;k<n2;k++)if(F[r][k]){Mon<4>a{pairs[j][0],pairs[j][1],pairs[k][0],pairs[k][1]};sort(a.begin(),a.end());check[4][qm[key(a)]]+=c*F[r][k];}}
  }
  for(int deg=0;deg<=4;deg++)for(int m=0;m<(int)check[deg].size();m++)assert(check[deg][m]%P==(input[deg][m]+P-value[deg][m])%P);
  out<<"TRIAL "<<trial<<" COEFFICIENT_IDENTITY_PASS 1\nY";for(auto c:Y)out<<' '<<c;out<<"\nINPUT";for(auto c:u)out<<' '<<c;out<<"\nOUTPUT";for(auto c:result)out<<' '<<c;out<<"\n";
  if(trial==trials-1){
   out<<"SAVED_COFACTOR_TRIAL "<<trial<<"\n";
   for(int r=0;r<EQS;r++){
    if(D[r])out<<"D "<<r<<' '<<D[r]<<'\n';
    for(int k=0;k<N;k++)if(L[r][k])out<<"L "<<r<<' '<<k<<' '<<L[r][k]<<'\n';
    for(int j=0;j<n2;j++)if(H[r][j])out<<"H "<<r<<' '<<pairs[j][0]<<' '<<pairs[j][1]<<' '<<H[r][j]<<'\n';
   }
  }
  return result;
 };
 out<<"GENERATORS "<<g<<" BASE 17 RANK3 "<<rank<<" FREE3 "<<nc<<" DET3 "<<det3<<" RESIDUAL4 "<<d<<" DET4 "<<det4<<" TRIALS "<<trials<<"\n";
 out<<"QUADRATIC_PIVOTS";for(int r:piv2)out<<' '<<r;out<<"\nCUBIC_ROWS";for(int r:rows)out<<' '<<r;out<<"\nCUBIC_PIVOTS";for(int r:piv)out<<' '<<r;out<<"\nQUARTIC_ALTERNATIVES";for(auto q:selected)out<<' '<<q.first<<' '<<q.second;out<<"\n";
 mt19937 rng(621103);vector<ulong>Y(nb),u(g);
 for(int trial=0;trial<trials;trial++){
  for(auto&c:u)c=rng()%P;for(auto&c:Y)c=rng()%P;
  if(trial==0){fill(Y.begin(),Y.end(),0);fill(u.begin(),u.end(),0);u[0]=1;}
  if(trial==1)Y=knownbase;
  if(trial==2)fill(Y.begin(),Y.end(),0);
  auto result=apply(u,Y,trial);if(trial==0){for(int i=0;i<g;i++)assert(result[i]==ulong(i==1));}
 }
 tick("GLOBAL_T_ACTION_CIRCUIT_TESTED");out.flush();assert(out.good());flint_cleanup();
}
