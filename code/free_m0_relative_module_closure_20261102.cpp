#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
using namespace std;using Row=vector<pair<int,int>>;const int P=32003,N=2820,K=141;
int inv(int a){int b=P-2,r=1;for(;b;b>>=1,a=int(int64_t(a)*a%P))if(b&1)r=int(int64_t(r)*a%P);return r;}
int main(int argc,char**argv){if(argc!=5)return 1;auto st=chrono::steady_clock::now();ifstream af(argv[1]),rf(argv[2]),nf(argv[3]);if(!af||!rf||!nf)return 2;int p,nb,k,nv;af>>p>>nb>>k>>nv;if(p!=P||nb!=20||k!=K||nv!=8)return 3;vector<uint64_t>base,normal;vector<int>bd,nd;auto gen=[&](auto&&self,int left,int first,uint64_t key,int vars,vector<uint64_t>&out)->void{if(!left){out.push_back(key);return;}for(int j=first;j<vars;j++)self(self,left-1,j,key+(1ULL<<(4*j)),vars,out);};for(int d=0;d<=3;d++){gen(gen,d,0,0,3,base);bd.resize(base.size(),d);}for(int d=0;d<=8;d++){gen(gen,d,0,0,8,normal);nd.resize(normal.size(),d);}unordered_map<uint64_t,int>bi;for(int i=0;i<20;i++)bi[base[i]]=i;int prod[20][20];for(int a=0;a<20;a++)for(int b=0;b<20;b++)prod[a][b]=bd[a]+bd[b]<=3?bi.at(base[a]+base[b]):-1;
vector<vector<Row>>M(11,vector<Row>(N));for(int j=0;j<8;j++)for(int b=0;b<K;b++){int len;af>>len;Row r;for(int z=0;z<len;z++){int i,c;af>>i>>c;r.emplace_back(i,c);}for(int a=0;a<20;a++)for(auto[i,c]:r){int ab=prod[a][i/K];if(ab>=0)M[j][a*K+b].emplace_back(ab*K+i%K,c);}}
for(int j=0;j<3;j++)for(int a=0;a<20;a++){int ab=prod[a][j+1];if(ab>=0)for(int b=0;b<K;b++)M[8+j][a*K+b]={{ab*K+b,1}};}
rf>>p>>nb>>k>>nv;if(p!=P||nb!=20||k!=K||nv!=136)return 4;vector<Row>input(136);for(auto&r:input){int len;rf>>len;for(int z=0;z<len;z++){int i,c;rf>>i>>c;r.emplace_back(i,c);}}nf>>p>>nb>>k>>nv;vector<int>B(141);for(int&i:B)nf>>i;
auto dense=[](const Row&r){vector<int>v(N);for(auto[i,c]:r)v[i]=c;return v;};auto sparse=[](const vector<int>&v){Row r;for(int i=0;i<N;i++)if(v[i])r.emplace_back(i,v[i]);return r;};
auto apply=[&](int j,const Row&r){vector<uint64_t>v(N);for(auto[i,c]:r)for(auto[l,a]:M[j][i])v[l]+=uint64_t(c)*a;vector<int>w(N);for(int i=0;i<N;i++)w[i]=v[i]%P;return w;};
vector<Row>relations;vector<int>pivot(N,-1);auto reduce=[&](vector<int>&v){for(int i=0;i<N;i++)if(v[i]&&pivot[i]>=0){int c=P-v[i];for(auto[j,a]:relations[pivot[i]])v[j]=(uint64_t(v[j])+uint64_t(c)*a)%P;}};auto add=[&](vector<int>v){reduce(v);int lead=0;while(lead<N&&!v[lead])lead++;if(lead==N)return false;int a=inv(v[lead]);for(int&c:v)c=int(int64_t(c)*a%P);pivot[lead]=relations.size();relations.push_back(sparse(v));return true;};
for(int i=0;i<135;i++)add(dense(input[i]));cerr<<"RAW rank "<<relations.size()<<" sec "<<chrono::duration<double>(chrono::steady_clock::now()-st).count()<<"\n";
for(int i=0;i<8;i++)for(int j=i+1;j<8;j++){for(int b=0;b<K;b++){auto v=apply(i,M[j][b]);auto w=apply(j,M[i][b]);for(int z=0;z<N;z++)v[z]=(v[z]+P-w[z])%P;add(std::move(v));}cerr<<"COMM "<<i<<" "<<j<<" rank "<<relations.size()<<" sec "<<chrono::duration<double>(chrono::steady_clock::now()-st).count()<<"\n";}
for(size_t q=0;q<relations.size();q++){Row r=relations[q];for(int j=0;j<11;j++)add(apply(j,r));if(q%50==0)cerr<<"CLOSE "<<q<<" rank "<<relations.size()<<" sec "<<chrono::duration<double>(chrono::steady_clock::now()-st).count()<<"\n";}

int initialRank=relations.size();long checkedCutoffs=0;int addedCutoffs=0;map<uint64_t,Row>previous;previous[0]={{0,1}};
for(int d=1;d<=11;d++){vector<uint64_t>level;gen(gen,d,0,0,8,level);map<uint64_t,Row>next;long support=0;for(uint64_t key:level){int j=7;while(!((key>>(4*j))&15))j--;uint64_t parent=key-(1ULL<<(4*j));auto it=previous.find(parent);if(it==previous.end())return 10;auto v=apply(j,it->second);reduce(v);Row r=sparse(v);support+=r.size();if(d>=8){for(int a=0;a<20;a++)if(bd[a]+d==11){vector<int>w(N);for(auto[i,c]:r){int ab=prod[a][i/K];if(ab>=0)w[ab*K+i%K]=c;}checkedCutoffs++;if(add(std::move(w)))addedCutoffs++;}}next[key]=std::move(r);}previous=std::move(next);cerr<<"DEGREE "<<d<<" monomials "<<level.size()<<" total support "<<support<<" cutoff new "<<addedCutoffs<<" rank "<<relations.size()<<" sec "<<chrono::duration<double>(chrono::steady_clock::now()-st).count()<<"\n";}
for(size_t q=initialRank;q<relations.size();q++){Row r=relations[q];for(int j=0;j<11;j++)add(apply(j,r));}cerr<<"CUTOFF checked "<<checkedCutoffs<<" added "<<addedCutoffs<<" final rank "<<relations.size()<<"\n";
int rank=relations.size();cerr<<"CLOSED rank "<<rank<<" dim "<<N-rank<<"\n";auto t=dense(input[135]);reduce(t);Row T=sparse(t);cerr<<"T support "<<T.size()<<"\n";
map<uint64_t,Row>eval;eval[0]=T;unordered_map<uint64_t,int>bmap;for(int i=0;i<141;i++)bmap[normal[B[i]]]=i;for(int i=1;i<141;i++){uint64_t key=normal[B[i]];int j=7;while(!((key>>(4*j))&15))j--;uint64_t parent=key-(1ULL<<(4*j));if(!eval.count(parent))return 5;auto v=apply(j,eval.at(parent));reduce(v);eval[key]=sparse(v);}vector<uint64_t>sq(N);for(auto[index,c]:T){int a=index/K,b=index%K;for(auto[l,z]:eval.at(normal[B[b]])){int ab=prod[a][l/K];if(ab>=0)sq[ab*K+l%K]+=uint64_t(c)*z;}}vector<int>ts(N);for(int i=0;i<N;i++)ts[i]=sq[i]%P;reduce(ts);auto TS=sparse(ts);cerr<<"T2 support "<<TS.size()<<"\n";
// Independent low-base truncation of the complete presentation.
vector<Row>save=relations;vector<int>savep=pivot;for(int i=10*K;i<N;i++){vector<int>v(N);v[i]=1;add(std::move(v));}cerr<<"LOW A3 dim "<<N-relations.size()<<"\n";if(N-relations.size()!=1303)return 6;relations=save;pivot=savep;
ofstream out(argv[4]);if(!out)return 7;out<<P<<" "<<N<<" "<<rank<<"\n";for(auto&r:relations){out<<r.size();for(auto[i,c]:r)out<<" "<<i<<" "<<c;out<<"\n";}out<<T.size();for(auto[i,c]:T)out<<" "<<i<<" "<<c;out<<"\n";out<<TS.size();for(auto[i,c]:TS)out<<" "<<i<<" "<<c;out<<"\n";out<<141<<"\n";for(int b:B){auto &r=eval.at(normal[b]);out<<r.size();for(auto[i,c]:r)out<<" "<<i<<" "<<c;out<<"\n";}out.flush();if(!out)return 8;cerr<<"ALL MODULE CLOSURE PASS "<<chrono::duration<double>(chrono::steady_clock::now()-st).count()<<"\n";
}
