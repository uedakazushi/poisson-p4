#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
using namespace std;using Row=vector<pair<int,int>>;
int main(int argc,char**argv){if(argc!=5)return 1;ifstream in(argv[1]),full(argv[2]),trace(argv[3]);ofstream out(argv[4]);if(!in||!full||!trace||!out)return 2;int p,n,N,F;in>>p>>n>>N>>F;vector<int>d(N);vector<vector<int>>shift(n,vector<int>(N));for(int&v:d)in>>v;for(auto&r:shift)for(int&v:r)in>>v;vector<Row>inputs(F);for(auto&r:inputs){int k;in>>k;for(int j=0;j<k;j++){int i,c;in>>i>>c;r.emplace_back(i,c);}}vector<int>unused(N);for(int&v:unused)in>>v;vector<uint64_t>key(N,UINT64_MAX);key[0]=0;for(int i=0;i<N;i++)for(int j=0;j<n;j++)if(shift[j][i]>=0)key[shift[j][i]]=key[i]+(1ULL<<(4*j));unordered_map<uint64_t,int>idx;for(int i=0;i<N;i++)idx[key[i]]=i;
int NN,B;full>>NN>>B;if(NN!=N)return 3;vector<Row>basis(B);for(auto&r:basis){int k;full>>k;for(int j=0;j<k;j++){int i,c;full>>i>>c;r.emplace_back(i,c);}}int BB;trace>>BB;if(BB!=B)return 4;unordered_map<uint64_t,Row>cache;vector<vector<int>>constant(B,vector<int>(F));long steps=0;auto add=[&](vector<int>&v,int id,int mul,int c){uint64_t ck=(uint64_t(id)<<32)|uint32_t(mul);auto it=cache.find(ck);if(it==cache.end()){Row row;for(auto[j,a]:basis[id])if(d[j]+d[mul]<=d.back())row.emplace_back(idx.at(key[j]+key[mul]),a);it=cache.emplace(ck,std::move(row)).first;}for(auto[j,a]:it->second)v[j]=(uint32_t(v[j])+uint32_t(c)*uint32_t(a))%32003U;};
for(int bi=0;bi<B;bi++){int kind,source,a,b,ma,mb,norm,count;trace>>kind>>source>>a>>b>>ma>>mb>>norm>>count;vector<int>v(N);if(kind==0){for(auto[j,c]:inputs[source])v[j]=c;constant[bi][source]=1;}else{if(a>=bi||b>=bi)return 5;add(v,a,ma,1);add(v,b,mb,32002);if(ma==0)constant[bi]=constant[a];if(mb==0)for(int j=0;j<F;j++)constant[bi][j]=(constant[bi][j]+32003-constant[b][j])%32003;}
for(int k=0;k<count;k++){int parent,mul,c;trace>>parent>>mul>>c;if(parent>=bi)return 6;add(v,parent,mul,c);if(mul==0)for(int j=0;j<F;j++)constant[bi][j]=(uint32_t(constant[bi][j])+uint32_t(c)*uint32_t(constant[parent][j]))%32003U;}steps+=count;vector<int>expected(N);for(auto[j,c]:basis[bi])expected[j]=c;for(int j=0;j<N;j++)if(uint64_t(v[j])*norm%32003!=expected[j]){cerr<<"bad "<<bi<<" "<<j<<"\n";return 7;}for(int&v:constant[bi])v=uint64_t(v)*norm%32003;}
out<<B<<" "<<F<<"\n";for(auto&r:constant){for(int v:r)out<<v<<" ";out<<"\n";}out.flush();if(!trace||!out)return 8;cout<<"TRACE PASS "<<B<<" steps "<<steps<<" cache "<<cache.size()<<"\n";
}
