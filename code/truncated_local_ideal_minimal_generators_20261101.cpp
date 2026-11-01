#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
using Row=vector<pair<int,int>>;
int main(int argc,char**argv){
 ifstream in(argv[1]);int p,n,N,F;in>>p>>n>>N>>F;
 vector<int>degree(N);for(int&i:degree)in>>i;
 vector<vector<int>>shift(n,vector<int>(N));for(auto&r:shift)for(int&i:r)in>>i;
 ifstream initial(argv[2]);int count;initial>>count;vector<Row>rows(count),basis(N);
 for(auto&r:rows){int len;initial>>len;r.resize(len);for(auto&[i,c]:r)initial>>i>>c;}
 auto powmod=[&](int a,int b){long long v=1;while(b){if(b&1)v=v*a%p;a=(long long)a*a%p;b>>=1;}return int(v);};
 auto add=[&](const Row&src){
  vector<int>r(N);for(auto[i,c]:src)r[i]=c;
  for(int i=0;i<N;i++)if(r[i]&&!basis[i].empty()){
   int a=r[i];for(auto[k,c]:basis[i]){r[k]=(r[k]-(long long)a*c)%p;if(r[k]<0)r[k]+=p;}}
  int pivot=0;while(pivot<N&&!r[pivot])pivot++;if(pivot==N)return false;
  int inv=powmod(r[pivot],p-2);for(int i=pivot;i<N;i++)if(r[i])basis[pivot].emplace_back(i,(long long)r[i]*inv%p);
  return true;
 };
 vector<int>selected;
 for(int d=2;d<=*max_element(degree.begin(),degree.end());d++){
  for(int i=0;i<N;i++)if(degree[i]==d-1&&!basis[i].empty())for(int j=0;j<n;j++){
   Row mult;for(auto[k,c]:basis[i])if(shift[j][k]>=0)mult.emplace_back(shift[j][k],c);add(mult);}
  int start=selected.size();for(int i=0;i<count;i++)if(degree[rows[i][0].first]==d&&add(rows[i]))selected.push_back(i);
  cerr<<"degree "<<d<<" new generators "<<selected.size()-start<<"\n";
 }
 ofstream out(argv[3]);for(int i:selected)out<<i<<"\n";
}
