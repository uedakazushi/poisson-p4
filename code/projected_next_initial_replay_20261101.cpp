// Exact filtered replay with only the next homogeneous degree projected.
// The supplied projection must be the quotient by an already certified
// initial ideal in that degree. No lower filtered tail is changed.
#include <algorithm>
#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>
#include <gmpxx.h>
using namespace std; using Q=mpq_class; using Row=vector<pair<int,Q>>;
struct Trace{int pivot,kind,source,var;};
int main(int argc,char**argv){
 auto started=chrono::steady_clock::now();
 if(argc!=6)return 1;
 ifstream in(argv[1]);int p,n,N,F;in>>p>>n>>N>>F;
 vector<int>deg(N);vector<vector<int>>shift(n,vector<int>(N));
 for(int&i:deg)in>>i;for(auto&r:shift)for(int&i:r)in>>i;
 int last=*max_element(deg.begin(),deg.end());
 int offset=0;while(deg[offset]<2)offset++;
 ifstream pi(argv[4]);int h,raw;pi>>h>>raw;
 vector<Row>proj(raw);for(auto&r:proj){int c;pi>>c;for(int j=0;j<c;j++){int k;Q a;pi>>k>>a;r.emplace_back(k,a);}}
 vector<vector<int>>topshift(n,vector<int>(N,-1));
 for(int j=0;j<n;j++)for(int k=0;k<N;k++)if(deg[k]==last)pi>>topshift[j][k];
 if(!pi){cerr<<"bad projection\n";return 2;}
 int M=N+h;Q product;
 auto add=[&](vector<Q>&dst,int k,const Q&a,const Q&b){mpq_mul(product.get_mpq_t(),a.get_mpq_t(),b.get_mpq_t());mpq_add(dst[k].get_mpq_t(),dst[k].get_mpq_t(),product.get_mpq_t());};
 ifstream cf(argv[2]);vector<Row>init(F);
 for(int i=0;i<F;i++){
  vector<Q>r(M);
  for(int k=offset;k<N;k++){cf>>r[k];}
  for(int k=0;k<raw;k++){Q a;cf>>a;if(a!=0)for(auto&[j,b]:proj[k])add(r,N+j,a,b);}
  for(int k=offset;k<M;k++)if(r[k]!=0)init[i].emplace_back(k,std::move(r[k]));
 }
 if(!cf){cerr<<"bad coefficients\n";return 2;}
 cerr<<"input ready "<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"\n";
 ifstream tr(argv[3]);vector<Trace>trace;vector<int>bid(N,-1);int id;
 while(tr>>id){Trace t;int c;tr>>t.pivot>>t.kind>>t.source>>t.var>>c;for(int j=0;j<c;j++){int dummy;tr>>dummy;}bid[t.pivot]=id;trace.push_back(t);}
 vector<Row>B(trace.size());vector<int>state(trace.size());int done=0;long long reductions=0;
 auto multiply=[&](const Row&src,int j,vector<Q>&r){
  for(auto&[k,c]:src){if(k>=N)continue;if(shift[j][k]>=0)r[shift[j][k]]+=c;else if(topshift[j][k]>=0)for(auto&[l,a]:proj[topshift[j][k]])add(r,N+l,c,a);}
 };
 function<void(int)> replay=[&](int i){
  if(state[i]==2)return;if(state[i]){cerr<<"cycle\n";exit(3);}state[i]=1;
  auto&t=trace[i];vector<Q>r(M);
  if(t.kind==0)for(auto&[k,c]:init[t.source])r[k]=c;
  else{replay(t.source);multiply(B[t.source],t.var,r);}
  for(int k=0;k<N;k++)if(r[k]!=0&&bid[k]>=0&&bid[k]<i){int j=bid[k];replay(j);Q a=r[k];for(auto&[l,c]:B[j]){mpq_mul(product.get_mpq_t(),a.get_mpq_t(),c.get_mpq_t());mpq_sub(r[l].get_mpq_t(),r[l].get_mpq_t(),product.get_mpq_t());}reductions++;}
  int first=0;while(first<M&&r[first]==0)first++;
  if(first!=t.pivot){cerr<<"pivot mismatch "<<i<<" "<<first<<" vs "<<t.pivot<<"\n";exit(4);}
  Q a=r[first];for(int k=first;k<M;k++)if(r[k]!=0)B[i].emplace_back(k,r[k]/a);
  state[i]=2;done++;if(done%100==0)cerr<<"replayed "<<done<<" terms "<<B[i].size()<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"\n";
 };
 for(int i=0;i<(int)trace.size();i++)replay(i);
 long long checks=0,nonzero=0;ofstream out(argv[5]);if(!out){cerr<<"bad output\n";return 6;}
 auto check=[&](vector<Q>&r){
  for(int k=0;k<N;k++)if(r[k]!=0){if(bid[k]<0){cerr<<"new lower residue "<<k<<"\n";exit(5);}Q a=r[k];for(auto&[l,c]:B[bid[k]]){mpq_mul(product.get_mpq_t(),a.get_mpq_t(),c.get_mpq_t());mpq_sub(r[l].get_mpq_t(),r[l].get_mpq_t(),product.get_mpq_t());}}
  Row rem;for(int k=N;k<M;k++)if(r[k]!=0)rem.emplace_back(k-N,r[k]);
  if(!rem.empty()){out<<checks<<" "<<rem.size();for(auto&[k,c]:rem)out<<" "<<k<<" "<<c;out<<"\n";nonzero++;}
  checks++;
  if(checks%500==0)cerr<<"closure "<<checks<<" nonzero "<<nonzero<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"\n";
 };
 for(auto&src:init){vector<Q>r(M);for(auto&[k,c]:src)r[k]=c;check(r);}
 for(int i=0;i<(int)B.size();i++)for(int j=0;j<n;j++){vector<Q>r(M);multiply(B[i],j,r);check(r);}
 cerr<<"PROJECTED SEVENTH CHECKS "<<checks<<" NONZERO "<<nonzero<<" reductions "<<reductions<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"\n";
 return nonzero?7:0;
}
