// Finite local standard basis, optional quotient normal forms, and full-basis export.
// Same arithmetic/reduction as truncated_local_standard_basis_20261101.cpp.
// Args: input initial-output [query-input normal-form-output].
// The full filtered basis is additionally written to initial-output + ".full".
// Exact standard basis in F_32003[x]/(x)^(D+1), D<16, n<=15.
// Input/output agree with truncated_local_ideal_io_20261101.py.
// The input order is increasing total degree, then decreasing lexicographic
// exponent order. All reductions stay inside the finite truncated ring.
// Only homogeneous initial parts are exported; no assertion about untruncated
// higher-degree equations follows from this computation.
// Shifted polynomial rows are cached. Product and resolved-chain criteria
// discard S-pairs only after a standard representation is available.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace std;
using Row=vector<pair<int,int>>;
struct Pair { int a,b,d; uint64_t lcm; long serial; };
struct Compare { bool operator()(const Pair&a,const Pair&b)const{return a.d!=b.d?a.d>b.d:a.serial>b.serial;} };
int main(int argc,char**argv){
 if(argc<3)return 2;ifstream in(argv[1]);if(!in){cerr<<"input open failed\n";return 3;}
 int p,n,N,F;in>>p>>n>>N>>F;if(p!=32003||n>15)return 4;
 vector<int>deg(N);vector<vector<int>>shift(n,vector<int>(N));for(int&v:deg)in>>v;for(auto&r:shift)for(int&v:r)in>>v;
 vector<Row>input(F);for(auto&r:input){int l;in>>l;r.resize(l);for(auto&[i,c]:r)in>>i>>c;}vector<int>tv(N);for(int&v:tv)in>>v;if(!in)return 5;
 int D=*max_element(deg.begin(),deg.end());if(D>=16||deg[0]!=0)return 6;for(int i=1;i<N;i++)if(deg[i]<deg[i-1])return 6;
 vector<uint64_t>key(N,UINT64_MAX);key[0]=0;for(int i=0;i<N;i++)for(int j=0;j<n;j++)if(shift[j][i]>=0){uint64_t k=key[i]+(1ULL<<(4*j));int h=shift[j][i];if(h<=i||h>=N||deg[h]!=deg[i]+1)return 7;if(key[h]!=UINT64_MAX&&key[h]!=k)return 7;key[h]=k;}
 unordered_map<uint64_t,int>idx;idx.reserve(N*2);for(int i=0;i<N;i++){if(key[i]==UINT64_MAX)return 8;idx[key[i]]=i;}
 auto divides=[&](uint64_t a,uint64_t b){for(int j=0;j<n;j++)if(((a>>(4*j))&15)>((b>>(4*j))&15))return false;return true;};
 auto lcm=[&](uint64_t a,uint64_t b){uint64_t c=0;for(int j=0;j<n;j++)c|=max((a>>(4*j))&15,(b>>(4*j))&15)<<(4*j);return c;};
 auto degree=[&](uint64_t a){int d=0;for(int j=0;j<n;j++)d+=(a>>(4*j))&15;return d;};
 auto power=[](int a,int b){int o=1;while(b){if(b&1)o=(uint64_t)o*a%32003;a=(uint64_t)a*a%32003;b>>=1;}return o;};
 vector<Row>basis;vector<int>reducer(N,-1);priority_queue<Pair,vector<Pair>,Compare>pairs;long serial=0,processed=0;
 unordered_set<uint64_t>resolved;
 auto pairkey=[](int a,int b){return(uint64_t(min(a,b))<<32)|uint32_t(max(a,b));};
 unordered_map<uint64_t,Row>cache;
 auto axpy=[&](vector<int>&r,int bi,uint64_t mul,int c){
  const Row* b=&basis[bi];
  if(mul){uint64_t ck=(uint64_t(bi)<<32)|uint32_t(idx.at(mul));auto it=cache.find(ck);if(it==cache.end()){Row shifted;int md=degree(mul);for(auto[k,a]:*b)if(deg[k]+md<=D)shifted.emplace_back(idx.at(key[k]+mul),a);it=cache.emplace(ck,std::move(shifted)).first;}b=&it->second;}
  for(auto[k,a]:*b)r[k]=(uint32_t(r[k])+uint32_t(c)*uint32_t(a))%32003U;
 };
 auto reduce=[&](vector<int>&r){for(int i=0;i<N;i++)if(r[i]&&reducer[i]>=0){int bi=reducer[i];axpy(r,bi,key[i]-key[basis[bi][0].first],32003-r[i]);}};
 auto add=[&](vector<int>&r){reduce(r);int lead=0;while(lead<N&&!r[lead])lead++;if(lead==N)return false;int iv=power(r[lead],32001);Row b;for(int i=lead;i<N;i++)if(r[i])b.emplace_back(i,(uint64_t)r[i]*iv%32003);int bi=basis.size();
  for(int j=0;j<bi;j++){int h=basis[j][0].first;auto L=lcm(key[h],key[lead]);int d=degree(L);if(d<=D){if(d<deg[h]+deg[lead])pairs.push(Pair{j,bi,d,L,serial++});else resolved.insert(pairkey(j,bi));}}
  basis.push_back(std::move(b));for(int i=lead;i<N;i++)if(reducer[i]<0&&divides(key[lead],key[i]))reducer[i]=bi;
  if(basis.size()%25==0)cerr<<"basis "<<basis.size()<<" pairs "<<pairs.size()<<" processed "<<processed<<" degree "<<deg[lead]<<"\n";return true;
 };
 for(const auto&b:input){vector<int>r(N);for(auto[i,c]:b)r[i]=c;add(r);}
 long skipped=0;
 while(!pairs.empty()){auto a=pairs.top();pairs.pop();bool chain=false;for(int k=0;k<(int)basis.size();k++)if(k!=a.a&&k!=a.b&&divides(key[basis[k][0].first],a.lcm)&&resolved.count(pairkey(k,a.a))&&resolved.count(pairkey(k,a.b))){chain=true;break;}
  if(chain)skipped++;else{vector<int>r(N);axpy(r,a.a,a.lcm-key[basis[a.a][0].first],1);axpy(r,a.b,a.lcm-key[basis[a.b][0].first],32002);add(r);}
  resolved.insert(pairkey(a.a,a.b));processed++;if(processed%5000==0)cerr<<"processed "<<processed<<" skipped "<<skipped<<" pairs "<<pairs.size()<<" cache "<<cache.size()<<"\n";
 }
 ofstream out(argv[2]);if(!out){cerr<<"output open failed\n";return 9;}out<<basis.size()<<"\n";for(auto&b:basis){int d=deg[b[0].first],l=0;for(auto[i,c]:b)if(deg[i]==d)l++;out<<l;for(auto[i,c]:b)if(deg[i]==d)out<<" "<<i<<" "<<c;out<<"\n";}out.flush();if(!out)return 10;

 if(argc>=5){ifstream queries(argv[3]);ofstream answers(argv[4]);if(!queries||!answers)return 11;int count;queries>>count;answers<<count<<"\n";
 for(int qi=0;qi<count;qi++){int length;queries>>length;vector<int>v(N);for(int k=0;k<length;k++){int column,c;queries>>column>>c;v[column]=c;}reduce(v);int total=0;for(int x:v)if(x)total++;answers<<total;for(int k=0;k<N;k++)if(v[k])answers<<" "<<k<<" "<<v[k];answers<<"\n";}answers.flush();if(!queries||!answers)return 12;}
 {ofstream save(string(argv[2])+".full");if(!save)return 13;save<<N<<" "<<basis.size()<<"\n";for(auto&b:basis){save<<b.size();for(auto[k,c]:b)save<<" "<<k<<" "<<c;save<<"\n";}save.flush();if(!save)return 14;}
 int len=count(reducer.begin(),reducer.end(),-1);cout<<"basis "<<basis.size()<<" length "<<len<<" pairs "<<processed<<"\n";for(int d=1;d<=D+1;d++){int ct=0;for(int i=0;i<N;i++)if(deg[i]<d&&reducer[i]<0)ct++;cout<<"modm "<<d<<" length "<<ct<<"\n";}reduce(tv);cout<<"t";for(int i=0;i<N;i++)if(tv[i])cout<<" "<<i<<":"<<tv[i];cout<<"\n";
}
