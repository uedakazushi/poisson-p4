// Independent reversed-monomial regular-representation replay of Sigma quad78.
#include <flint/nmod_mat.h>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <vector>
using namespace std;const int D=28;const ulong prime=32003;using Ex=array<int,3>;using Mon=array<int,6>;
struct Small{nmod_mat_t x;Small(){nmod_mat_init(x,D,D,prime);}Small(const Small&o){nmod_mat_init(x,D,D,prime);nmod_mat_set(x,o.x);}Small&operator=(const Small&o){nmod_mat_set(x,o.x);return*this;}~Small(){nmod_mat_clear(x);}};
Small times(const Small&a,const Small&b){Small c;nmod_mat_mul(c.x,a.x,b.x);return c;}
Small one(){Small c;for(int i=0;i<D;i++)nmod_mat_entry(c.x,i,i)=1;return c;}
uint32_t powmod(uint32_t a,uint32_t n){uint64_t x=a,r=1;while(n){if(n&1)r=r*x%prime;x=x*x%prime;n>>=1;}return r;}
void monomials(vector<Mon>&out,Mon m,int start,int left){if(!left){out.push_back(m);return;}for(int i=start;i<6;i++){m[i]++;monomials(out,m,i,left-1);m[i]--;}}
uint32_t get32(istream&f){unsigned char b[4];f.read((char*)b,4);assert(f.good());return(uint32_t)b[0]|((uint32_t)b[1]<<8)|((uint32_t)b[2]<<16)|((uint32_t)b[3]<<24);}
int main(int argc,char**argv){assert(argc==5);auto start=chrono::steady_clock::now();ifstream ops(argv[1]),qfile(argv[2]);int dd,pp;ops>>dd>>pp;assert(dd==D&&pp==(int)prime);Small MP,ME;
for(auto*M:{&MP,&ME})for(int i=0;i<D;i++)for(int j=0;j<D;j++){ulong v;ops>>v;assert(v<prime);nmod_mat_entry(M->x,i,j)=v;}
Small eP=times(MP,ME),Pe=times(ME,MP);assert(nmod_mat_equal(eP.x,Pe.x));Small tmp;
for(int i=0;i<D;i++)for(int j=0;j<D;j++)nmod_mat_entry(tmp.x,i,j)=(20*nmod_mat_entry(MP.x,i,j)+3*nmod_mat_entry(ME.x,i,j)+(i==j?20:0))%prime;
Small MA=times(MP,tmp);ulong factor=prime-powmod(180,prime-2);for(int i=0;i<D;i++)for(int j=0;j<D;j++)nmod_mat_entry(MA.x,i,j)=nmod_mat_entry(MA.x,i,j)*factor%prime;
int nq,nmon,nvars;qfile>>nq>>nmon>>nvars;assert(nq==4&&nmon==21&&nvars==3);using Poly=map<Ex,ulong>;vector<Poly>Q(84);array<int,3>mx{};
for(auto&f:Q){int n;qfile>>n;for(int k=0;k<n;k++){Ex ex;ulong c;qfile>>c>>ex[0]>>ex[1]>>ex[2];assert(c<prime&&ex[0]>=0&&ex[1]>=0&&ex[2]>=0);f[ex]=c;for(int v=0;v<3;v++)mx[v]=max(mx[v],ex[v]);}}
vector<Small>A(mx[0]+1),E(mx[1]+1),P(mx[2]+1);A[0]=E[0]=P[0]=one();for(int i=1;i<=mx[0];i++)A[i]=times(A[i-1],MA);for(int i=1;i<=mx[1];i++)E[i]=times(E[i-1],ME);for(int i=1;i<=mx[2];i++)P[i]=times(P[i-1],MP);
map<Ex,Small>cache;vector<Small>V(84);for(int n=0;n<84;n++)for(auto term:Q[n]){const auto&ex=term.first;if(!cache.count(ex))cache.emplace(ex,times(times(A[ex[0]],E[ex[1]]),P[ex[2]]));auto&op=cache.at(ex);for(int i=0;i<D;i++)for(int j=0;j<D;j++)nmod_mat_entry(V[n].x,i,j)=(nmod_mat_entry(V[n].x,i,j)+term.second*nmod_mat_entry(op.x,i,j))%prime;}
vector<Mon>m2,m4;monomials(m2,Mon{},0,2);monomials(m4,Mon{},0,4);assert(m2.size()==21&&m4.size()==126);map<Mon,int>index;for(int i=0;i<21;i++)index[m2[i]]=i;
nmod_mat_t T;nmod_mat_init(T,2352,3528,prime);
for(int i=0;i<4;i++)for(int shift=0;shift<21;shift++)for(int target=0;target<126;target++){Mon beta{};bool valid=true;for(int v=0;v<6;v++){beta[v]=m4[target][v]-m2[shift][v];if(beta[v]<0)valid=false;}if(!valid)continue;const auto&op=V[i*21+index.at(beta)];for(int r=0;r<D;r++)for(int c=0;c<D;c++)nmod_mat_entry(T,(i*21+shift)*D+r,target*D+c)=nmod_mat_entry(op.x,r,c);}
// Compare every entry after our opposite monomial ordering. This checks both
// coordinate blocks and operator orientation, not merely the observed rank.
ifstream author(argv[3],ios::binary);assert(get32(author)==2352&&get32(author)==3528&&get32(author)==prime);uint64_t checked=0;
for(int authorrow=0;authorrow<2352;authorrow++){int qi=(authorrow/D)/21,shift=(authorrow/D)%21,r=authorrow%D;int ownrow=(qi*21+20-shift)*D+r;for(int authorcol=0;authorcol<3528;authorcol++){int target=authorcol/D,c=authorcol%D;ulong v=get32(author);assert(nmod_mat_entry(T,ownrow,(125-target)*D+c)==v);checked++;}}
char trailing;assert(!author.get(trailing));double buildsec=chrono::duration<double>(chrono::steady_clock::now()-start).count();auto rs=chrono::steady_clock::now();slong rank=nmod_mat_rank(T);double ranksec=chrono::duration<double>(chrono::steady_clock::now()-rs).count();assert(rank==2184);
ofstream out(argv[4]);out<<"{\"all_passed\":true,\"independent_FLINT_small_operator_arithmetic\":true,\"opposite_degree2_and4_monomial_order\":true,\"all_author_block_entries_permuted_exact\":"<<checked<<",\"own_shape\":[2352,3528],\"own_single_rank\":"<<rank<<",\"build_seconds\":"<<buildsec<<",\"rank_seconds\":"<<ranksec<<"}\n";cout<<"OWN_REVERSED_SOURCE_ALL_"<<checked<<"_ENTRIES_RANK_"<<rank<<"_PASS\n";nmod_mat_clear(T);}
