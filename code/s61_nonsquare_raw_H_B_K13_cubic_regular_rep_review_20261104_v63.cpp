// Fixed selected minor only: independent monomial reductions, PLU and FLINTdet.
#include <flint/nmod_mat.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;const uint64_t p=32003;const int n=221,D=13;
int main(){auto start=chrono::steady_clock::now();
 int pp,nr,nc,dd;cin>>pp>>nr>>nc>>dd;assert(pp==int(p)&&nr==17&&nc==24&&dd==13);
 vector<uint64_t>g(13),f(17*24*13);for(auto&x:g){cin>>x;assert(x<p);}for(auto&x:f){cin>>x;assert(x<p);}
 vector<int>cols(n),perm(n);vector<uint64_t>diag(n),lu(n*n);for(auto&x:cols)cin>>x;for(auto&x:perm)cin>>x;
 for(auto&x:diag){cin>>x;assert(x&&x<p);}for(auto&x:lu){cin>>x;assert(x<p);}uint64_t expected;cin>>expected;assert(expected&&expected<p);
 string extra;assert(!(cin>>extra));auto cp=cols,rp=perm;sort(cp.begin(),cp.end());sort(rp.begin(),rp.end());
 assert(cp.front()>=0&&cp.back()<312&&adjacent_find(cp.begin(),cp.end())==cp.end());for(int i=0;i<n;i++)assert(rp[i]==i);
 uint64_t red[25][13]={};for(int i=0;i<13;i++)red[i][i]=1;
 for(int k=13;k<25;k++)for(int i=0;i<13;i++)for(int j=0;j<13;j++)red[k][i]=(red[k][i]+p-(g[j]*red[k-13+j][i])%p)%p;
 nmod_mat_t M;nmod_mat_init(M,n,n,p);
 for(int i=0;i<17;i++)for(int j=0;j<n;j++){int block=cols[j]/13,power=cols[j]%13;
  for(int k=0;k<13;k++){uint64_t s=0;for(int m=0;m<13;m++)s+=f[(i*24+block)*13+m]*red[m+power][k];nmod_mat_entry(M,i*13+k,j)=s%p;}}
 int parity=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)parity^=perm[i]>perm[j];uint64_t pivotdet=1;
 for(int i=0;i<n;i++){assert(lu[i*n+i]==diag[i]);pivotdet=pivotdet*diag[i]%p;}
 if(parity)pivotdet=(p-pivotdet)%p;assert(pivotdet==expected);
 long checks=0;for(int i=0;i<n;i++)for(int j=0;j<n;j++){
  uint64_t sum=0;for(int k=0;k<=min(i,j);k++)sum+=(k==i?1:lu[i*n+k])*lu[k*n+j];
  assert(sum%p==nmod_mat_entry(M,perm[i],j));checks++;
 }
 // Free the supplied LU before the library's independent determinant copy.
 vector<uint64_t>().swap(lu);uint64_t actual=nmod_mat_det(M);assert(actual==expected&&actual==30542);nmod_mat_clear(M);
 cout<<"{\"all_passed\":true,\"prime\":32003,\"regular_shape\":[221,312],\"fixed_selected_minor_size\":221,\"all_PLU_entries\":"<<checks
 <<",\"independent_FLINT_fixed_determinant\":"<<actual<<",\"pivot_product_equals_det\":true,\"source_column_multiplication_not_point_evaluation\":true,\"explicit_scalar_peak_bound\":"<<(2*n*n+int(f.size())+25*13+3*n)
 <<",\"no_rank_search\":true,\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
}
