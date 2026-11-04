// Exact rational-function rank and a primitive polynomial one-dimensional kernel.
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <gmp.h>
#include <flint/fmpz_mpoly_q.h>
using namespace std;
fmpz_mpoly_ctx_t ctx;int variables;
struct R{
 fmpz_mpoly_q_t v;R(){fmpz_mpoly_q_init(v,ctx);}R(const R&p):R(){fmpz_mpoly_q_set(v,p.v,ctx);}
 R(R&&p):R(){fmpz_mpoly_q_swap(v,p.v,ctx);}R&operator=(R p){fmpz_mpoly_q_swap(v,p.v,ctx);return *this;}
 ~R(){fmpz_mpoly_q_clear(v,ctx);}bool zero()const{return fmpz_mpoly_q_is_zero(v,ctx);}
};
struct P{fmpz_mpoly_t v;P(){fmpz_mpoly_init(v,ctx);}P(const P&p):P(){fmpz_mpoly_set(v,p.v,ctx);}~P(){fmpz_mpoly_clear(v,ctx);}};
void readpoly(istream&in,fmpz_mpoly_t p){long n;in>>n;assert(in&&n>=0);fmpz_mpoly_zero(p,ctx);fmpz_t z;fmpz_init(z);vector<ulong>e(variables);for(long i=0;i<n;i++){string c;in>>c;for(auto&v:e)in>>v;assert(in&&fmpz_set_str(z,c.c_str(),10)==0);fmpz_mpoly_push_term_fmpz_ui(p,z,e.data(),ctx);}fmpz_mpoly_sort_terms(p,ctx);fmpz_mpoly_combine_like_terms(p,ctx);fmpz_clear(z);}
int main(int argc,char**argv){assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);int n,m;in>>variables>>n>>m;assert(variables==5&&n>0&&m>1);fmpz_mpoly_ctx_init(ctx,variables,ORD_DEGLEX);{
 vector<vector<R>> A(n,vector<R>(m));for(auto&row:A)for(auto&q:row){readpoly(in,fmpz_mpoly_q_numref(q.v));readpoly(in,fmpz_mpoly_q_denref(q.v));assert(!fmpz_mpoly_is_zero(fmpz_mpoly_q_denref(q.v),ctx));fmpz_mpoly_q_canonicalise(q.v,ctx);}string extra;assert(!(in>>extra));auto original=A;
 int rank=0;vector<int>pivots;R factor,product,sum;
 for(int col=0;col<m&&rank<n;col++){
  int p=rank;while(p<n&&A[p][col].zero())p++;if(p==n)continue;if(p!=rank)swap(A[p],A[rank]);
  for(int row=rank+1;row<n;row++)if(!A[row][col].zero()){
   fmpz_mpoly_q_div(factor.v,A[row][col].v,A[rank][col].v,ctx);
   for(int j=col+1;j<m;j++){fmpz_mpoly_q_mul(product.v,factor.v,A[rank][j].v,ctx);fmpz_mpoly_q_sub(A[row][j].v,A[row][j].v,product.v,ctx);}fmpz_mpoly_q_zero(A[row][col].v,ctx);
  }
  pivots.push_back(col);rank++;cerr<<"PIVOT "<<col<<" RANK "<<rank<<endl;
 }
 for(int row=rank;row<n;row++)for(int col=0;col<m;col++)assert(A[row][col].zero());assert(rank==m-1);
 int freecol=0;for(;freecol<m;freecol++){bool used=false;for(int p:pivots)used|=p==freecol;if(!used)break;}assert(freecol<m);
 vector<R> kernel(m);fmpz_mpoly_q_one(kernel[freecol].v,ctx);
 for(int row=rank-1;row>=0;row--){int p=pivots[row];fmpz_mpoly_q_zero(sum.v,ctx);for(int col=p+1;col<m;col++){fmpz_mpoly_q_mul(product.v,A[row][col].v,kernel[col].v,ctx);fmpz_mpoly_q_add(sum.v,sum.v,product.v,ctx);}fmpz_mpoly_q_div(kernel[p].v,sum.v,A[row][p].v,ctx);fmpz_mpoly_q_neg(kernel[p].v,kernel[p].v,ctx);}
 P common,gcd,quot;fmpz_mpoly_one(common.v,ctx);
 for(auto&q:kernel){assert(fmpz_mpoly_gcd(gcd.v,common.v,fmpz_mpoly_q_denref(q.v),ctx));assert(fmpz_mpoly_divides(quot.v,fmpz_mpoly_q_denref(q.v),gcd.v,ctx));fmpz_mpoly_mul(common.v,common.v,quot.v,ctx);}
 vector<P> polynomial(m);for(int col=0;col<m;col++){assert(fmpz_mpoly_divides(quot.v,common.v,fmpz_mpoly_q_denref(kernel[col].v),ctx));fmpz_mpoly_mul(polynomial[col].v,quot.v,fmpz_mpoly_q_numref(kernel[col].v),ctx);}
 fmpz_mpoly_zero(gcd.v,ctx);for(auto&p:polynomial)assert(fmpz_mpoly_gcd(gcd.v,gcd.v,p.v,ctx));assert(!fmpz_mpoly_is_zero(gcd.v,ctx));
 for(int col=0;col<m;col++){assert(fmpz_mpoly_divides(quot.v,polynomial[col].v,gcd.v,ctx));fmpz_mpoly_set(polynomial[col].v,quot.v,ctx);fmpz_mpoly_set(fmpz_mpoly_q_numref(kernel[col].v),quot.v,ctx);fmpz_mpoly_one(fmpz_mpoly_q_denref(kernel[col].v),ctx);}
 for(int row=0;row<n;row++){fmpz_mpoly_q_zero(sum.v,ctx);for(int col=0;col<m;col++){fmpz_mpoly_q_mul(product.v,original[row][col].v,kernel[col].v,ctx);fmpz_mpoly_q_add(sum.v,sum.v,product.v,ctx);}assert(sum.zero());}
 const char* names[]={"u","v","X","W","B"};out<<"{\"all_passed\":true,\"rank\":"<<rank<<",\"free_column\":"<<freecol<<",\"pivot_columns\":[";for(int i=0;i<rank;i++){if(i)out<<',';out<<pivots[i];}out<<"],\"primitive_polynomial_kernel\":[";
 for(int i=0;i<m;i++){if(i)out<<',';char*s=fmpz_mpoly_get_str_pretty(polynomial[i].v,names,ctx);out<<'"'<<s<<'"';flint_free(s);}out<<"]}\n";
 cerr<<"EXACT_ONE_KERNEL_AND_ALL_RESIDUAL_IDENTITIES_PASS\n";
 }fmpz_mpoly_ctx_clear(ctx);flint_cleanup();}
