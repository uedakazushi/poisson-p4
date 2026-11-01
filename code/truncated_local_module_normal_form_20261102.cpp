// Reduce sparse query rows by a previously verified full finite local basis.
// Args: variable-count max-degree full-basis query-input answer-output metadata-output.
// Monomial order is increasing degree, then combinations-with-replacement order.
// This is a normal-form evaluator, not a new ideal-membership or basis certificate.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
using namespace std;using Row=vector<pair<int,int>>;
int main(int argc,char**argv){if(argc!=7)return 1;int n=stoi(argv[1]),D=stoi(argv[2]);vector<uint64_t>keys;vector<int>deg;auto rec=[&](auto&&self,int left,int first,uint64_t key,int d)->void{if(!left){keys.push_back(key);deg.push_back(d);return;}for(int j=first;j<n;j++)self(self,left-1,j,key+(1ULL<<(4*j)),d);};for(int d=0;d<=D;d++)rec(rec,d,0,0,d);int N=keys.size();unordered_map<uint64_t,int>ix;ix.reserve(2*N);for(int i=0;i<N;i++)ix[keys[i]]=i;ifstream in(argv[3]);int count,NN;in>>NN>>count;if(NN!=N)return 2;vector<Row>B(count);vector<int>lead(N,-1),red(N,-1);for(int j=0;j<count;j++){int len;in>>len;for(int h=0;h<len;h++){int i,c;in>>i>>c;B[j].emplace_back(i,c);}if(B[j].empty()||B[j][0].second!=1)return 3;int l=B[j][0].first;if(lead[l]<0)lead[l]=j;}if(!in)return 4;for(int i=0;i<N;i++){if(lead[i]>=0){red[i]=lead[i];continue;}for(int j=0;j<n;j++)if((keys[i]>>(4*j))&15){int q=ix.at(keys[i]-(1ULL<<(4*j)));if(red[q]>=0){red[i]=red[q];break;}}}
ifstream query(argv[4]);ofstream out(argv[5]);if(!query||!out)return 5;int Q;query>>Q;out<<Q<<"\n";unordered_map<uint64_t,Row>cache;
for(int q=0;q<Q;q++){int nn;query>>nn;vector<int>v(N);for(int j=0;j<nn;j++){int i,c;query>>i>>c;v[i]=c;}for(int i=0;i<N;i++)if(v[i]&&red[i]>=0){int b=red[i],h=B[b][0].first;uint64_t mul=keys[i]-keys[h];auto ck=(uint64_t(b)<<32)|uint32_t(ix.at(mul));auto it=cache.find(ck);if(it==cache.end()){Row r;for(auto[j,c]:B[b])if(deg[j]+deg[i]-deg[h]<=D)r.emplace_back(ix.at(keys[j]+mul),c);it=cache.emplace(ck,std::move(r)).first;}int c=32003-v[i];for(auto[j,a]:it->second)v[j]=(uint32_t(v[j])+uint32_t(c)*uint32_t(a))%32003U;}int len=0;for(int c:v)if(c)len++;out<<len;for(int i=0;i<N;i++)if(v[i])out<<" "<<i<<" "<<v[i];out<<"\n";if(q%500==0)cerr<<"query "<<q<<" cache "<<cache.size()<<"\n";}out.flush();if(!query||!out)return 6;ofstream meta(argv[6]);meta<<"columns "<<N<<" basis "<<count<<" standard "<<count_if(red.begin(),red.end(),[](int x){return x<0;})<<"\n";meta.flush();if(!meta)return 7;return 0;}
