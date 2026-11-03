#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#include <flint/fmpq.h>
#include <array>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <iostream>
#include <string>
#include <cassert>
#include <chrono>
using namespace std;const int d=28;const ulong p=32003;nmod_t mod;
struct Mat{nmod_mat_t m;Mat(){nmod_mat_init(m,d,d,p);}Mat(const Mat&a):Mat(){nmod_mat_set(m,a.m);}Mat(Mat&&a):Mat(){nmod_mat_swap(m,a.m);}Mat&operator=(Mat a){nmod_mat_swap(m,a.m);return *this;}~Mat(){nmod_mat_clear(m);}};
void addmul(Mat&dest,const Mat&src,ulong c){for(int i=0;i<d;i++)for(int j=0;j<d;j++)nmod_mat_entry(dest.m,i,j)=nmod_add(nmod_mat_entry(dest.m,i,j),nmod_mul(c,nmod_mat_entry(src.m,i,j),mod),mod);}
using Ex=array<int,3>;struct Term{Ex e;ulong c;};ulong rat(const string&s){fmpq_t q;fmpq_init(q);assert(!fmpq_set_str(q,s.c_str(),10));fmpq_canonicalise(q);ulong den=fmpz_fdiv_ui(fmpq_denref(q),p);assert(den);ulong z=nmod_mul(fmpz_fdiv_ui(fmpq_numref(q),p),nmod_inv(den,mod),mod);fmpq_clear(q);return z;}
int main(int argc,char**argv){assert(argc==5);auto start=chrono::steady_clock::now();nmod_init(&mod,p);ifstream matrices(argv[1]),raw(argv[2]),polys(argv[3]);int n,prime;matrices>>n>>prime;assert(n==d&&prime==p);Mat MP,ME,I,tmp; nmod_mat_one(I.m);for(Mat*A:{&MP,&ME})for(int i=0;i<d;i++)for(int j=0;j<d;j++){ulong v;matrices>>v;nmod_mat_entry(A->m,i,j)=v;}Mat pe,ep;nmod_mat_mul(pe.m,MP.m,ME.m);nmod_mat_mul(ep.m,ME.m,MP.m);assert(nmod_mat_equal(pe.m,ep.m));int nr,nc,flag;raw>>nr>>nc>>flag;assert(nr==24&&nc==56&&flag==0);vector<vector<Term>>A(nr*nc);set<Ex>keys;int maxA=0,maxE=0,maxP=0;long rawterms=0;
for(auto&q:A){int nt;raw>>nt;for(int k=0;k<nt;k++){string cs;Term t;raw>>cs>>t.e[0]>>t.e[1]>>t.e[2];t.c=rat(cs);q.push_back(t);keys.insert(t.e);maxA=max(maxA,t.e[0]);maxE=max(maxE,t.e[1]);maxP=max(maxP,t.e[2]);rawterms++;}}
int pn,pm,pk;polys>>pn>>pm>>pk;assert(pn==1&&pm==3&&pk==0);vector<vector<pair<array<int,2>,ulong>>>model(3);
for(auto&q:model){int nt;polys>>nt;for(int k=0;k<nt;k++){string cs;array<int,2>e;polys>>cs>>e[0]>>e[1];q.push_back({e,rat(cs)});maxP=max(maxP,e[0]);maxE=max(maxE,e[1]);}}
Mat aMat;addmul(tmp,MP,20);addmul(tmp,ME,3);addmul(tmp,I,20);nmod_mat_mul(aMat.m,MP.m,tmp.m);nmod_mat_scalar_mul(aMat.m,aMat.m,p-nmod_inv(180,mod));
vector<Mat>AP(maxA+1),EP(maxE+1),PP(maxP+1);nmod_mat_one(AP[0].m);nmod_mat_one(EP[0].m);nmod_mat_one(PP[0].m);for(int i=1;i<=maxA;i++)nmod_mat_mul(AP[i].m,AP[i-1].m,aMat.m);for(int i=1;i<=maxE;i++)nmod_mat_mul(EP[i].m,EP[i-1].m,ME.m);for(int i=1;i<=maxP;i++)nmod_mat_mul(PP[i].m,PP[i-1].m,MP.m);
vector<Mat>op(3);for(int j=0;j<3;j++)for(auto&t:model[j]){Mat z;nmod_mat_mul(z.m,PP[t.first[0]].m,EP[t.first[1]].m);addmul(op[j],z,t.second);}assert(nmod_mat_is_zero(op[0].m)&&nmod_mat_is_zero(op[1].m));Mat UI,Ucheck;assert(nmod_mat_inv(UI.m,op[2].m));nmod_mat_mul(Ucheck.m,UI.m,op[2].m);assert(nmod_mat_equal(Ucheck.m,I.m));
map<Ex,Mat>cache;for(auto ex:keys){Mat z,v;nmod_mat_mul(z.m,AP[ex[0]].m,EP[ex[1]].m);nmod_mat_mul(v.m,z.m,PP[ex[2]].m);cache.emplace(ex,move(v));}
const int rows=24*d,cols=56*d;nmod_mat_t C,Rref;nmod_mat_init(C,rows,cols,p);nmod_mat_init(Rref,rows,cols,p);
for(int i=0;i<24;i++)for(int j=0;j<56;j++){Mat block;for(auto&t:A[i*56+j])addmul(block,cache.at(t.e),t.c);for(int r=0;r<d;r++)for(int c=0;c<d;c++)nmod_mat_entry(C,i*d+r,j*d+c)=nmod_mat_entry(block.m,r,c);}
nmod_mat_set(Rref,C);long rank=nmod_mat_rref(Rref);assert(rank==rows);vector<int>selected;for(int i=0;i<rows;i++){int j=0;while(j<cols&&!nmod_mat_entry(Rref,i,j))j++;assert(j<cols);selected.push_back(j);}nmod_mat_t S,Inv,Check;nmod_mat_init(S,rows,rows,p);nmod_mat_init(Inv,rows,rows,p);nmod_mat_init(Check,rows,rows,p);for(int i=0;i<rows;i++)for(int j=0;j<rows;j++)nmod_mat_entry(S,i,j)=nmod_mat_entry(C,i,selected[j]);assert(nmod_mat_inv(Inv,S));nmod_mat_mul(Check,S,Inv);assert(nmod_mat_is_one(Check));ofstream out(argv[4]);out<<rows<<' '<<p<<'\n';for(int i=0;i<rows;i++){for(int j=0;j<rows;j++)out<<nmod_mat_entry(Inv,i,j)<<' ';out<<'\n';}
cout<<"{\"all_exact_native_checks_passed\":true,\"field_prime\":32003,\"algebra_dimension\":28,\"raw_source_rows\":24,\"raw_source_columns\":56,\"raw_source_terms\":"<<rawterms<<",\"all_raw_denominators_prime_units\":true,\"commuting_multiplication_matrices\":true,\"N_and_M_operator_relations_verified\":true,\"FULLU_operator_inverse_verified\":true,\"raw_block_matrix_shape\":[672,1568],\"raw_block_rank\":"<<rank<<",\"selected_square_inverse_identity_verified\":true,\"selected_columns\":[";for(int i=0;i<rows;i++){if(i)cout<<',';cout<<selected[i];}cout<<"],\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}"<<endl;
nmod_mat_clear(Check);nmod_mat_clear(Inv);nmod_mat_clear(S);nmod_mat_clear(Rref);nmod_mat_clear(C);}
