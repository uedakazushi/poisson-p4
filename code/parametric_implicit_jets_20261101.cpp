// Exact parameter-polynomial implicit jets. Internal research engine.
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <flint/fmpq_mpoly.h>
using namespace std;
fmpq_mpoly_ctx_t ctx;
int nv;
struct Poly {
 fmpq_mpoly_struct p[1];
 Poly(){fmpq_mpoly_init(p,ctx);}
 Poly(const Poly&o):Poly(){fmpq_mpoly_set(p,o.p,ctx);}
 Poly(Poly&&o):Poly(){fmpq_mpoly_swap(p,o.p,ctx);}
 Poly&operator=(Poly o){fmpq_mpoly_swap(p,o.p,ctx);return *this;}
 ~Poly(){fmpq_mpoly_clear(p,ctx);}
 bool zero()const{return fmpq_mpoly_is_zero(p,ctx);}
};
void read(istream&in,Poly&a){
 long count;in>>count;assert(in&&count>=0);fmpq_t c;fmpq_init(c);vector<ulong>e(nv);
 for(long i=0;i<count;i++){string v;in>>v;assert(fmpq_set_str(c,v.c_str(),10)==0);fmpq_canonicalise(c);for(auto&x:e)in>>x;assert(in);fmpq_mpoly_push_term_fmpq_ui(a.p,c,e.data(),ctx);}
 fmpq_mpoly_sort_terms(a.p,ctx);fmpq_mpoly_combine_like_terms(a.p,ctx);fmpq_clear(c);
}
void write(ostream&out,const Poly&a){
 long count=fmpq_mpoly_length(a.p,ctx);out<<count<<'\n';fmpq_t c;fmpq_init(c);vector<ulong>e(nv);
 for(long i=0;i<count;i++){fmpq_mpoly_get_term_coeff_fmpq(c,a.p,i,ctx);char*v=fmpq_get_str(nullptr,10,c);out<<v;flint_free(v);fmpq_mpoly_get_term_exp_ui(e.data(),a.p,i,ctx);for(auto x:e)out<<' '<<x;out<<'\n';}
 fmpq_clear(c);assert(out);
}
int main(int argc,char**argv){
 assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);assert(in&&out);
 int F,np,n,D,r;in>>F>>np>>n>>D>>r;assert(in&&F>0&&n>0&&r>0);nv=np+n;
 fmpq_mpoly_ctx_init(ctx,nv,ORD_LEX);flint_set_num_threads(2);
 auto start=chrono::steady_clock::now();vector<int>rows(r);for(auto&i:rows)in>>i;
 Poly den,tmp,product,sum;read(in,den);
 vector<vector<Poly>>J(F,vector<Poly>(45)),T(45,vector<Poly>(n)),U(45,vector<Poly>(r));
 for(auto&M:{&J,&T,&U})for(auto&row:*M)for(auto&a:row)read(in,a);
 map<pair<int,int>,vector<pair<int,Poly>>>H;
 for(int row=0;row<F;row++){int count;in>>count;for(int k=0;k<count;k++){int i,j;in>>i>>j;Poly c;read(in,c);if(i>j)swap(i,j);H[{i,j}].emplace_back(row,std::move(c));}}
 assert(in);vector<vector<Poly>>W(D+1,vector<Poly>(45)),A(D+1,vector<Poly>(F));
 for(int i=0;i<45;i++)for(int j=0;j<n;j++)if(!T[i][j].zero()){
  fmpq_mpoly_gen(tmp.p,np+j,ctx);fmpq_mpoly_mul(product.p,tmp.p,T[i][j].p,ctx);fmpq_mpoly_add(W[1][i].p,W[1][i].p,product.p,ctx);
 }
 for(int d=2;d<=D;d++){
  vector<Poly>Q(F);
  for(auto&[ij,coefficients]:H){
   fmpq_mpoly_zero(sum.p,ctx);
   for(int a=1;a<d;a++)if(!W[a][ij.first].zero()&&!W[d-a][ij.second].zero()){
    fmpq_mpoly_mul(product.p,W[a][ij.first].p,W[d-a][ij.second].p,ctx);fmpq_mpoly_add(sum.p,sum.p,product.p,ctx);
   }
   if(sum.zero())continue;
   for(auto&[row,c]:coefficients){fmpq_mpoly_mul(product.p,c.p,sum.p,ctx);fmpq_mpoly_add(Q[row].p,Q[row].p,product.p,ctx);}
  }
  for(int i=0;i<45;i++)for(int j=0;j<r;j++)if(!U[i][j].zero()&&!Q[rows[j]].zero()){
   fmpq_mpoly_mul(product.p,U[i][j].p,Q[rows[j]].p,ctx);fmpq_mpoly_sub(W[d][i].p,W[d][i].p,product.p,ctx);
  }
  long terms=0;
  for(int i=0;i<F;i++){
   fmpq_mpoly_mul(A[d][i].p,den.p,Q[i].p,ctx);
   for(int j=0;j<45;j++)if(!J[i][j].zero()&&!W[d][j].zero()){
    fmpq_mpoly_mul(product.p,J[i][j].p,W[d][j].p,ctx);fmpq_mpoly_add(A[d][i].p,A[d][i].p,product.p,ctx);
   }
   terms+=fmpq_mpoly_length(A[d][i].p,ctx);
  }
  for(int i:rows)assert(A[d][i].zero());
  cerr<<"degree "<<d<<" residual terms "<<terms<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';
 }
 out<<F<<' '<<np<<' '<<n<<' '<<D<<'\n';
 for(int d=2;d<=D;d++){for(auto&a:W[d])write(out,a);for(auto&a:A[d])write(out,a);}
 cerr<<"PARAMETRIC EXACT JETS PASS "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';
}
