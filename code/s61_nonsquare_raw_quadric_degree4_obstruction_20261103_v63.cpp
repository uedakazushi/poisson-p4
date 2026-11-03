#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#include <flint/fmpq.h>
#include <array>
#include <map>
#include <fstream>
#include <iostream>
#include <vector>
#include <cassert>
#include <chrono>
using namespace std;using Ex=array<int,3>;struct Term{Ex e;ulong c;};
vector<Ex>mons(int D){vector<Ex>v;for(int s=0;s<=D;s++)for(int a=0;a<=s;a++)for(int b=0;b<=s-a;b++)v.push_back({a,b,s-a-b});return v;}
ulong pw(ulong x,int n,nmod_t mod){ulong r=1;while(n--){r=nmod_mul(r,x,mod);}return r;}
int main(int argc,char**argv){assert(argc==3);int D=stoi(argv[2]);assert(D>=1&&D<=4);ulong prime=32003;nmod_t mod;nmod_init(&mod,prime);ifstream f(argv[1]);int n,m,k;f>>n>>m>>k;assert(n==24&&m==56&&k==0);vector<vector<Term>>A(n*m);fmpq_t c;fmpq_init(c);long terms=0;
for(auto&q:A){int z;f>>z;for(int i=0;i<z;i++){string cs;Term t;f>>cs>>t.e[0]>>t.e[1]>>t.e[2];assert(f&&!fmpq_set_str(c,cs.c_str(),10));fmpq_canonicalise(c);ulong den=fmpz_fdiv_ui(fmpq_denref(c),prime);assert(den);t.c=nmod_mul(fmpz_fdiv_ui(fmpq_numref(c),prime),nmod_inv(den,mod),mod);if(t.c){q.push_back(t);terms++;}}}
auto small=mons(D),big=mons(D+7);map<Ex,long>ids;for(long i=0;i<(long)big.size();i++)ids[big[i]]=i;long rows=n*big.size(),cols=m*small.size();nmod_mat_t E;nmod_mat_init(E,rows,cols+n,prime);auto started=chrono::steady_clock::now();
for(int i=0;i<n;i++)for(int j=0;j<m;j++)for(long z=0;z<(long)small.size();z++)for(auto&t:A[i*m+j]){Ex ee={t.e[0]+small[z][0],t.e[1]+small[z][1],t.e[2]+small[z][2]};long r=i*big.size()+ids.at(ee),cc=j*small.size()+z;nmod_mat_entry(E,r,cc)=nmod_add(nmod_mat_entry(E,r,cc),t.c,mod);}
long fac[]={1,1,2,6};for(int a=0;a<=3;a++)for(int b=0;b<=3-a;b++){int z=3-a-b;long integer=fac[3]/fac[a]/fac[b]/fac[z];ulong cf=integer;cf=nmod_mul(cf,pw(5,a,mod),mod);cf=nmod_mul(cf,pw(prime-6,b,mod),mod);cf=nmod_mul(cf,pw(prime-4,z,mod),mod);Ex ee={5,b,a};for(int i=0;i<n;i++)nmod_mat_entry(E,i*big.size()+ids.at(ee),cols+i)=cf;}
cout<<"COEFFICIENT_MATRIX_READY "<<rows<<" "<<cols<<"\n"<<flush;long rank=nmod_mat_rref(E);long rankC=0;vector<bool>sol(n,true);
for(long i=0;i<rank;i++){long pivot=0;while(pivot<cols+n&&!nmod_mat_entry(E,i,pivot))pivot++;assert(pivot<cols+n);if(pivot<cols)rankC++;else for(int t=0;t<n;t++)if(nmod_mat_entry(E,i,cols+t))sol[t]=false;}
int solved=0;for(bool b:sol)solved+=b;double secs=chrono::duration<double>(chrono::steady_clock::now()-started).count();cout<<"{\"prime\":"<<prime<<",\"cofactor_degree_bound\":"<<D<<",\"coefficient_matrix_rows\":"<<rows<<",\"unknown_coefficients_per_target\":"<<cols<<",\"rank_coefficient_matrix\":"<<rankC<<",\"rank_augmented_all24\":"<<rank<<",\"individually_solvable_targets\":"<<solved<<",\"all24_solvable\":"<<(solved==24?"true":"false")<<",\"elapsed_seconds\":"<<secs<<",\"GF_only\":true,\"QQ_claim\":false}"<<endl;
nmod_mat_clear(E);fmpq_clear(c);flint_cleanup();}
