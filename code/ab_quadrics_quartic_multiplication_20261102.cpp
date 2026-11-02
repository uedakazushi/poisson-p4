// Rank lower bounds for all degree-four multiples of the 274 rational AB quadrics.
// Integer products are bounded by 32003^2 below INT_MAX; every stored row is normalized.
// The early exit only stops after an independently known upper bound is reached.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include <unordered_map>
using namespace std;using Row=vector<pair<int,int>>;
constexpr int P=32003;int inv(int x){int a=x,b=P-2,y=1;while(b){if(b&1)y=y*a%P;a=a*a%P;b>>=1;}return y;}
int key4(int a,int b,int c,int d){array<int,4>x={a,b,c,d};sort(x.begin(),x.end());return x[0]+64*(x[1]+64*(x[2]+64*x[3]));}
struct Group{int want=-1,cols=-1,rank=0;vector<int> mons;vector<Row> rows;};
int main(){int n;cin>>n;assert(n==45);vector<int>x(n),y(n);for(int i=0;i<n;i++)cin>>x[i]>>y[i];int z;cin>>z;vector<vector<array<int,3>>>quad(z);for(auto &r:quad){int k;cin>>k;while(k--){array<int,3>e;cin>>e[0]>>e[1]>>e[2];r.push_back(e);}}
 map<pair<int,int>,Group>g;int t;cin>>t;while(t--){int a,b,c,r;cin>>a>>b>>c>>r;auto &h=g[{a,b}];h.cols=c;h.want=c-r;}
 for(int i=0;i<n;i++)for(int j=i;j<n;j++)for(int k=j;k<n;k++)for(int l=k;l<n;l++)g.at({x[i]+x[j]+x[k]+x[l],y[i]+y[j]+y[k]+y[l]}).mons.push_back(key4(i,j,k,l));
 unordered_map<int,int>index;for(auto &[w,h]:g){assert((int)h.mons.size()==h.cols);sort(h.mons.begin(),h.mons.end());for(int i=0;i<h.cols;i++)index[h.mons[i]]=i;}
 for(auto &r:quad){int a=r[0][0],b=r[0][1];for(int k=0;k<n;k++)for(int l=k;l<n;l++){
  auto &h=g.at({x[a]+x[b]+x[k]+x[l],y[a]+y[b]+y[k]+y[l]});Row row;for(auto e:r)row.push_back({index.at(key4(e[0],e[1],k,l)),e[2]});sort(row.begin(),row.end());h.rows.push_back(std::move(row));}}
 int total=0,count=0;for(auto &item:g){auto &w=item.first;auto &h=item.second;vector<Row>B(h.cols);vector<int>pending;auto store=[&](Row r){int t=inv(r[0].second);for(auto &v:r)v.second=v.second*t%P;B[r[0].first]=std::move(r);h.rank++;};
 for(int i=0;i<(int)h.rows.size();i++){auto&r=h.rows[i];if(B[r[0].first].empty())store(r);else pending.push_back(i);}
 vector<int>d(h.cols);for(int id:pending){if(h.rank==h.want)break;fill(d.begin(),d.end(),0);auto&r=h.rows[id];for(auto [i,c]:r)d[i]=c;
 for(int j=r[0].first;j<h.cols;j++){int v=d[j];if(!v)continue;if(B[j].empty()){Row nr;for(int k=j;k<h.cols;k++)if(d[k])nr.push_back({k,d[k]});store(std::move(nr));break;}
 for(auto [k,c]:B[j]){int u=(d[k]-v*c)%P;d[k]=u<0?u+P:u;}}}
 assert(h.rank<=h.want);cout<<w.first<<' '<<w.second<<' '<<h.cols<<' '<<h.rows.size()<<' '<<h.rank<<' '<<h.want<<'\n';total+=h.rank;count++;if(count%25==0)cerr<<count<<" groups; rank "<<total<<'\n';}
 cerr<<"TOTAL "<<total<<" GROUPS "<<count<<endl;}
