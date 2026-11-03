#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpq_mat.h>
#include <flint/nmod_mat.h>
#include <flint/nmod.h>
#include <array>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include <cassert>
#include <chrono>
using namespace std;using Ex=array<int,4>;struct T{int r,c;Ex e;string q;};
int main(int argc,char**argv){assert(argc==4);ifstream src(argv[1]),seed(argv[2]);int n,m,k;src>>n>>m>>k;assert(n==24&&m==56&&k==0);vector<Ex>mon;for(int a=0;a<=8;a++)for(int b=0;b<=8-a;b++)for(int c=0;c<=8-a-b;c++)mon.push_back({a,b,c,8-a-b-c});assert(mon.size()==165);map<Ex,int>ix;for(int i=0;i<165;i++)ix[mon[i]]=i;vector<Ex>cofs={Ex{0,0,0,1},Ex{0,0,1,0},Ex{0,1,0,0},Ex{1,0,0,0}};fmpq_mat_t C;fmpq_mat_init(C,3960,224);fmpq_t q;fmpq_init(q);
for(int i=0;i<n;i++)for(int j=0;j<m;j++){int z;src>>z;for(int t=0;t<z;t++){string qs;int a,b,c;src>>qs>>a>>b>>c;assert(a+b+c<=7&&!fmpq_set_str(q,qs.c_str(),10));fmpq_canonicalise(q);Ex base={a,b,c,7-a-b-c};for(int k=0;k<4;k++){Ex ee;for(int v=0;v<4;v++)ee[v]=base[v]+cofs[k][v];fmpq_set(fmpq_mat_entry(C,i*165+ix.at(ee),j*4+k),q);}}}
seed>>n>>m>>k;assert(n==24&&m==112&&k==0);vector<T>ts;vector<int>deg(m,-1);for(int i=0;i<n;i++)for(int j=0;j<m;j++){int z;seed>>z;for(int t=0;t<z;t++){T x;x.r=i;x.c=j;seed>>x.q>>x.e[0]>>x.e[1]>>x.e[2]>>x.e[3];int d=0;for(int v:x.e)d+=v;if(deg[j]<0)deg[j]=d;assert(deg[j]==d);ts.push_back(x);}}
vector<int>selected;for(int j=0;j<m;j++)if(deg[j]==8)selected.push_back(j);assert(selected.size()==3);map<int,int>dest;for(int j=0;j<3;j++)dest[selected[j]]=j;fmpq_mat_t B;fmpq_mat_init(B,3960,3);for(auto&t:ts)if(deg[t.c]==8){assert(!fmpq_set_str(q,t.q.c_str(),10));fmpq_canonicalise(q);fmpq_set(fmpq_mat_entry(B,t.r*165+ix.at(t.e),dest[t.c]),q);}
ulong p=32003;nmod_t mod;nmod_init(&mod,p);nmod_mat_t CT;nmod_mat_init(CT,224,3960,p);for(int i=0;i<3960;i++)for(int j=0;j<224;j++){auto z=fmpq_mat_entry(C,i,j);ulong den=fmpz_fdiv_ui(fmpq_denref(z),p);assert(den);nmod_mat_entry(CT,j,i)=nmod_mul(fmpz_fdiv_ui(fmpq_numref(z),p),nmod_inv(den,mod),mod);}assert(nmod_mat_rref(CT)==224);vector<int>piv;for(int j=0;j<224;j++){int i=0;while(i<3960&&!nmod_mat_entry(CT,j,i))i++;assert(i<3960);piv.push_back(i);}
fmpq_mat_t A,BS,X,Y;fmpq_mat_init(A,224,224);fmpq_mat_init(BS,224,3);fmpq_mat_init(X,224,3);fmpq_mat_init(Y,3960,3);for(int i=0;i<224;i++){for(int j=0;j<224;j++)fmpq_set(fmpq_mat_entry(A,i,j),fmpq_mat_entry(C,piv[i],j));for(int j=0;j<3;j++)fmpq_set(fmpq_mat_entry(BS,i,j),fmpq_mat_entry(B,piv[i],j));}auto t0=chrono::steady_clock::now();assert(fmpq_mat_solve_dixon(X,A,BS));fmpq_mat_mul(Y,C,X);int passed=0;vector<int>ok;for(int j=0;j<3;j++){bool eq=true;for(int i=0;i<3960;i++)if(!fmpq_equal(fmpq_mat_entry(Y,i,j),fmpq_mat_entry(B,i,j)))eq=false;if(eq){passed++;ok.push_back(selected[j]);}}
ofstream out(argv[3]);out<<"224 3\n";for(int i=0;i<224;i++)for(int j=0;j<3;j++){char* qs=fmpq_get_str(nullptr,10,fmpq_mat_entry(X,i,j));out<<qs<<'\n';flint_free(qs);}cout<<"{\"candidate_degree8_columns\":3,\"full_QQ_membership_verified_columns\":"<<passed<<",\"all_3960_by3_coefficient_identities\":"<<(passed==3?"true":"false")<<",\"elapsed_solve_verify_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-t0).count()<<",\"selected_candidate_columns\":[";for(int j=0;j<3;j++){if(j)cout<<',';cout<<selected[j];}cout<<"],\"selected_original_coefficient_rows\":[";for(int i=0;i<224;i++){if(i)cout<<',';cout<<piv[i];}cout<<"]}"<<endl;return passed==3?0:3;}
