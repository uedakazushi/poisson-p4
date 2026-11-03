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
int main(int argc,char**argv){assert(argc==4);ifstream src(argv[1]),seed(argv[2]);int n,m,k;src>>n>>m>>k;assert(n==24&&m==56&&k==0);vector<Ex>mon;for(int a=0;a<=7;a++)for(int b=0;b<=7-a;b++)for(int c=0;c<=7-a-b;c++)mon.push_back({a,b,c,7-a-b-c});assert(mon.size()==120);map<Ex,int>ix;for(int i=0;i<120;i++)ix[mon[i]]=i;fmpq_mat_t C;fmpq_mat_init(C,2880,56);fmpq_t q;fmpq_init(q);
for(int i=0;i<n;i++)for(int j=0;j<m;j++){int z;src>>z;for(int t=0;t<z;t++){string qs;int a,b,c;src>>qs>>a>>b>>c;assert(a+b+c<=7&&!fmpq_set_str(q,qs.c_str(),10));fmpq_canonicalise(q);fmpq_set(fmpq_mat_entry(C,i*120+ix.at({a,b,c,7-a-b-c}),j),q);}}
seed>>n>>m>>k;assert(n==24&&m==112&&k==0);vector<T>ts;vector<int>deg(m,-1);for(int i=0;i<n;i++)for(int j=0;j<m;j++){int z;seed>>z;for(int t=0;t<z;t++){T x;x.r=i;x.c=j;seed>>x.q>>x.e[0]>>x.e[1]>>x.e[2]>>x.e[3];int d=0;for(int v:x.e)d+=v;if(deg[j]<0)deg[j]=d;assert(deg[j]==d);ts.push_back(x);}}
vector<int>selected;for(int j=0;j<m;j++)if(deg[j]==7)selected.push_back(j);assert(selected.size()==44);map<int,int>dest;for(int j=0;j<44;j++)dest[selected[j]]=j;fmpq_mat_t B;fmpq_mat_init(B,2880,44);for(auto&t:ts)if(deg[t.c]==7){assert(!fmpq_set_str(q,t.q.c_str(),10));fmpq_canonicalise(q);fmpq_set(fmpq_mat_entry(B,t.r*120+ix.at(t.e),dest[t.c]),q);}
ulong p=32003;nmod_t mod;nmod_init(&mod,p);nmod_mat_t CT;nmod_mat_init(CT,56,2880,p);for(int i=0;i<2880;i++)for(int j=0;j<56;j++){auto z=fmpq_mat_entry(C,i,j);ulong den=fmpz_fdiv_ui(fmpq_denref(z),p);assert(den);nmod_mat_entry(CT,j,i)=nmod_mul(fmpz_fdiv_ui(fmpq_numref(z),p),nmod_inv(den,mod),mod);}assert(nmod_mat_rref(CT)==56);vector<int>piv;for(int j=0;j<56;j++){int i=0;while(i<2880&&!nmod_mat_entry(CT,j,i))i++;assert(i<2880);piv.push_back(i);}
fmpq_mat_t A,BS,X,Y;fmpq_mat_init(A,56,56);fmpq_mat_init(BS,56,44);fmpq_mat_init(X,56,44);fmpq_mat_init(Y,2880,44);for(int i=0;i<56;i++){for(int j=0;j<56;j++)fmpq_set(fmpq_mat_entry(A,i,j),fmpq_mat_entry(C,piv[i],j));for(int j=0;j<44;j++)fmpq_set(fmpq_mat_entry(BS,i,j),fmpq_mat_entry(B,piv[i],j));}auto t0=chrono::steady_clock::now();assert(fmpq_mat_solve_dixon(X,A,BS));fmpq_mat_mul(Y,C,X);int passed=0;vector<int>ok;for(int j=0;j<44;j++){bool eq=true;for(int i=0;i<2880;i++)if(!fmpq_equal(fmpq_mat_entry(Y,i,j),fmpq_mat_entry(B,i,j)))eq=false;if(eq){passed++;ok.push_back(selected[j]);}}
ofstream out(argv[3]);out<<"56 44\n";for(int i=0;i<56;i++)for(int j=0;j<44;j++){char* qs=fmpq_get_str(nullptr,10,fmpq_mat_entry(X,i,j));out<<qs<<'\n';flint_free(qs);}cout<<"{\"candidate_degree7_columns\":44,\"full_QQ_membership_verified_columns\":"<<passed<<",\"all_2880_by44_coefficient_identities\":"<<(passed==44?"true":"false")<<",\"elapsed_solve_verify_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-t0).count()<<",\"selected_candidate_columns\":[";for(int j=0;j<44;j++){if(j)cout<<',';cout<<selected[j];}cout<<"],\"selected_original_coefficient_rows\":[";for(int i=0;i<56;i++){if(i)cout<<',';cout<<piv[i];}cout<<"]}"<<endl;return passed==44?0:3;}
