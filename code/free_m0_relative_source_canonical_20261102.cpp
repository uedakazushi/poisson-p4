#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
using namespace std; using Row=vector<pair<int,int>>;
int main(int argc,char**argv){if(argc!=4)return 1;auto st=chrono::steady_clock::now();ifstream f(argv[1],ios::binary),sel(argv[2]);if(!f||!sel)return 2;auto rd=[&](){int x;f.read((char*)&x,4);return x;};int p=rd(),D=rd(),nb=rd(),nrow=rd();if(p!=32003||D!=10||nb!=20||nrow!=89)return 3;vector<int>ids(nrow);for(int&i:ids)i=rd();vector<uint64_t>nk,bk;vector<int>nd,bd;auto gen=[&](auto&&self,int left,int first,uint64_t key,int vars,vector<uint64_t>&out)->void{if(!left){out.push_back(key);return;}for(int j=first;j<vars;j++)self(self,left-1,j,key+(1ULL<<(4*j)),vars,out);};for(int d=0;d<=10;d++){int old=nk.size();gen(gen,d,0,0,8,nk);nd.resize(nk.size(),d);}for(int d=0;d<=3;d++){gen(gen,d,0,0,3,bk);bd.resize(bk.size(),d);}unordered_map<uint64_t,int>ni,bi;for(int i=0;i<(int)nk.size();i++)ni[nk[i]]=i;for(int i=0;i<20;i++)bi[bk[i]]=i;vector<vector<Row>>H(20,vector<Row>(89));vector<int>size(20);for(int b=0;b<20;b++){if(rd()!=b)return 4;size[b]=rd();for(auto&r:H[b]){int len=rd();r.reserve(len);for(int j=0;j<len;j++){int i=rd(),c=rd();r.emplace_back(i,c);}}}if(!f)return 5;
vector<int>lead(89);for(int j=0;j<89;j++){lead[j]=H[0][j][0].first;if(H[0][j][0].second!=1)return 6;}vector<int>reducer(nk.size(),-1);auto divides=[](uint64_t a,uint64_t b){for(int j=0;j<8;j++)if(((a>>(4*j))&15)>((b>>(4*j))&15))return false;return true;};for(int k=0;k<(int)nk.size();k++)for(int j=0;j<89;j++)if(divides(nk[lead[j]],nk[k])){reducer[k]=j;break;}
int fcount;sel>>fcount;int discard;for(int j=0;j<fcount;j++)sel>>discard;sel>>discard;for(int j=0;j<discard;j++){int v;sel>>v;}

const int C=89;vector<vector<uint64_t>>v(20);for(int a=0;a<20;a++)v[a].resize(size[a]*C);for(int j=0;j<89;j++)v[0][lead[j]*C+j]=1;long long terms=0;long reductions=0;vector<vector<int>>constants(89,vector<int>(89));
for(int a=0;a<20;a++){for(int k=0;k<size[a];k++){uint64_t*src=v[a].data()+k*C;bool any=false;for(int c=0;c<C;c++){src[c]%=32003;any|=src[c]!=0;}if(!any||reducer[k]<0)continue;int j=reducer[k];uint64_t mul=nk[k]-nk[lead[j]];int md=nd[k]-nd[lead[j]];uint64_t coef[C];for(int c=0;c<C;c++){coef[c]=src[c];src[c]=0;}if(a==0&&mul==0)for(int c=0;c<C;c++)constants[c][j]=(constants[c][j]+coef[c])%32003;
for(int b=0;b<20;b++){if(bd[a]+bd[b]>3)continue;int ab=bi.at(bk[a]+bk[b]);for(auto[l,z]:H[b][j]){if(nd[l]+md+bd[ab]>10)break;if(b==0&&l==lead[j])continue;int target=ni.at(nk[l]+mul);uint64_t*dest=v[ab].data()+target*C;uint64_t scalar=32003-z;for(int c=0;c<C;c++)dest[c]+=scalar*coef[c];terms++;}}reductions++;}
int non=0;for(auto z:v[a])if(z%32003)non++;cerr<<"BLOCK "<<a<<" reduced "<<reductions<<" terms "<<terms<<" remaining "<<non<<" sec "<<chrono::duration<double>(chrono::steady_clock::now()-st).count()<<"\n";}
vector<int>B;for(int k=0;k<size[0];k++)if(reducer[k]<0)B.push_back(k);if(B.size()!=141)return 7;ofstream out(argv[3]);if(!out)return 8;out<<32003<<" "<<20<<" "<<89<<" "<<141<<"\n";for(int k:B)out<<k<<" ";out<<"\n";for(int j:lead)out<<j<<" ";out<<"\n";for(int a=0;a<20;a++)for(int j=0;j<89;j++){for(int k:B)out<<v[a][k*C+j]%32003<<" ";out<<"\n";}for(auto&r:constants){for(int x:r)out<<x<<" ";out<<"\n";}out.flush();if(!out)return 9;
cout<<"PASS reductions "<<reductions<<" terms "<<terms<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-st).count()<<"\n";
}
