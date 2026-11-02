#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <deque>
// Sparse multiplication maps retain scalar normal-form tails.
using namespace std;using Row=vector<pair<int,int>>;
struct Node{Row row;int kind,source,var;};
Row readrow(istream&in){int l;in>>l;Row r;while(l--){int i,c;in>>i>>c;r.emplace_back(i,c);}return r;}
int main(int argc,char**argv){
 ifstream in(argv[1]);ofstream trace(argv[2]);int p,n,N,F;in>>p>>n>>N>>F;
 vector<int>deg(N),bid(N,-1);for(int&i:deg)in>>i;
 vector<vector<Row>>shift(n,vector<Row>(N));for(auto&v:shift)for(auto&r:v)r=readrow(in);
 vector<Row>base(N);const int cutoff=*max_element(deg.begin(),deg.end());
 vector<deque<Node>>queue(cutoff+1);for(int j=0;j<F;j++){Row r=readrow(in);if(!r.empty())queue[deg[r[0].first]].push_back({move(r),0,j,-1});}
 auto power=[&](int a,int b){long long out=1;while(b){if(b&1)out=out*a%p;a=(long long)a*a%p;b>>=1;}return(int)out;};
 int rank=0;for(int stage=0;stage<=cutoff;stage++)while(!queue[stage].empty()){
  Node src=move(queue[stage].front());queue[stage].pop_front();vector<int>r(N),red;
  for(auto[i,c]:src.row)r[i]=c;
  for(int i=0;i<N;i++)if(r[i]&&!base[i].empty()){int a=r[i];red.push_back(bid[i]);for(auto[k,c]:base[i]){r[k]=(r[k]-(long long)a*c)%p;if(r[k]<0)r[k]+=p;}}
  int pi=0;while(pi<N&&!r[pi])pi++;if(pi==N)continue;
  if(deg[pi]>stage){queue[deg[pi]].push_back(move(src));continue;}
  int iv=power(r[pi],p-2);Row b;for(int i=pi;i<N;i++)if(r[i])b.emplace_back(i,(long long)r[i]*iv%p);
  base[pi]=b;bid[pi]=rank;trace<<rank<<" "<<pi<<" "<<src.kind<<" "<<src.source<<" "<<src.var<<" "<<red.size();for(int id:red)trace<<" "<<id;trace<<"\n";
  for(int j=0;j<n;j++){vector<int>w(N);for(auto[k,c]:b)for(auto[l,d]:shift[j][k])w[l]=(w[l]+(long long)c*d)%p;Row v;for(int k=0;k<N;k++)if(w[k])v.emplace_back(k,w[k]);if(!v.empty()){if(stage==cutoff)return 2;queue[stage+1].push_front({move(v),1,rank,j});}}
  rank++;
 }
 ofstream order(argv[3]);for(int i=0;i<N;i++)if(bid[i]>=0)order<<i<<" "<<bid[i]<<"\n";
 if(!trace||!order)return 3;
 cout<<rank<<"\n";
}
