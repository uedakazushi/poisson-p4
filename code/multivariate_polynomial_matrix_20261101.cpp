// Exact scale*L-C*B over QQ in the supplied parameter variables.
#include <flint/fmpq_mpoly.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
using namespace std;
using Poly = fmpq_mpoly_struct;
struct Array {
 vector<Poly> a;
 const fmpq_mpoly_ctx_struct* ctx;
 Array(size_t n,const fmpq_mpoly_ctx_t c):a(n),ctx(c){for(auto&x:a)fmpq_mpoly_init(&x,ctx);}
 ~Array(){for(auto&x:a)fmpq_mpoly_clear(&x,ctx);}
 Poly* operator[](size_t i){return &a[i];}
};
void readpoly(istream&in,Poly*x,const fmpq_mpoly_ctx_t ctx,int n){
 long terms;in>>terms;fmpq_t c;fmpq_init(c);vector<ulong>ex(n);
 for(long k=0;k<terms;k++){string q;in>>q;if(fmpq_set_str(c,q.c_str(),10))throw runtime_error("bad rational");fmpq_canonicalise(c);for(auto&e:ex)in>>e;fmpq_mpoly_push_term_fmpq_ui(x,c,ex.data(),ctx);}
 fmpq_mpoly_sort_terms(x,ctx);fmpq_mpoly_combine_like_terms(x,ctx);fmpq_clear(c);if(!in)throw runtime_error("bad input");
}
void writepoly(ostream&out,Poly*x,const fmpq_mpoly_ctx_t ctx,int n){
 long len=fmpq_mpoly_length(x,ctx);out<<len<<'\n';fmpq_t c;fmpq_init(c);vector<ulong>ex(n);
 for(long k=0;k<len;k++){fmpq_mpoly_get_term_coeff_fmpq(c,x,k,ctx);char*q=fmpq_get_str(nullptr,10,c);out<<q;flint_free(q);fmpq_mpoly_get_term_exp_ui(ex.data(),x,k,ctx);for(auto e:ex)out<<' '<<e;out<<'\n';}fmpq_clear(c);
}
int main(int argc,char**argv){try{
 if(argc!=3)throw runtime_error("input output required");ifstream in(argv[1]);if(!in)throw runtime_error("input open");int n,nr,nk,nc;in>>n>>nr>>nk>>nc;fmpq_mpoly_ctx_t ctx;fmpq_mpoly_ctx_init(ctx,n,ORD_LEX);
 {Array L(nr*nc,ctx),C(nr*nk,ctx),B(nk*nc,ctx),temp(1,ctx),scale(1,ctx);vector<vector<int>> cn(nr),bn(nk);
 auto readmat=[&](Array&A,int rows,int cols,vector<vector<int>>*support){for(int i=0;i<rows;i++){int count;in>>count;for(int k=0;k<count;k++){int j;in>>j;readpoly(in,A[i*cols+j],ctx,n);if(support)(*support)[i].push_back(j);}}};
 readmat(L,nr,nc,nullptr);readmat(C,nr,nk,&cn);readmat(B,nk,nc,&bn);readpoly(in,scale[0],ctx,n);if(!in)throw runtime_error("truncated input");
 for(int i=0;i<nr;i++){for(int j=0;j<nc;j++)fmpq_mpoly_mul(L[i*nc+j],L[i*nc+j],scale[0],ctx);for(int k:cn[i])for(int j:bn[k]){fmpq_mpoly_mul(temp[0],C[i*nk+k],B[k*nc+j],ctx);fmpq_mpoly_sub(L[i*nc+j],L[i*nc+j],temp[0],ctx);}if(i%25==0)cerr<<"row "<<i<<"/"<<nr<<"\n";}
 ofstream out(argv[2]);if(!out)throw runtime_error("output open failed");out<<nr<<' '<<nc<<' '<<n<<'\n';long total=0;for(int i=0;i<nr;i++){vector<int>supp;for(int j=0;j<nc;j++)if(!fmpq_mpoly_is_zero(L[i*nc+j],ctx))supp.push_back(j);out<<supp.size()<<'\n';for(int j:supp){out<<j<<'\n';writepoly(out,L[i*nc+j],ctx,n);total+=fmpq_mpoly_length(L[i*nc+j],ctx);}}out.flush();if(!out)throw runtime_error("output write failed");cerr<<"terms "<<total<<"\n";
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();return 0;
}catch(const exception&e){cerr<<e.what()<<"\n";return 1;}}
