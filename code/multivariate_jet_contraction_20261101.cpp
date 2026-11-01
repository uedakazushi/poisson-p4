// Exact QQ multivariate coefficient contraction, with a supplied homogeneous top-degree projection.
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
 if(argc!=3)throw runtime_error("input output required");ifstream in(argv[1]);if(!in)throw runtime_error("input open");int n,coords,eq,fullcols,cols,degree;in>>n>>coords>>eq>>fullcols>>cols>>degree;int nw=degree-1,npairs=degree/2;if(degree<3||degree>16)throw runtime_error("bad degree");fmpq_mpoly_ctx_t ctx;fmpq_mpoly_ctx_init(ctx,n,ORD_LEX);flint_set_num_threads(2);
 {
 vector<Array*>W;vector<int>sizes;vector<vector<vector<int>>>nz;
 for(int d=0;d<nw;d++){int count;in>>count;sizes.push_back(count);W.push_back(new Array(coords*count,ctx));nz.emplace_back(coords);for(int i=0;i<coords;i++)for(int j=0;j<count;j++){readpoly(in,(*W[d])[i*count+j],ctx,n);if(!fmpq_mpoly_is_zero((*W[d])[i*count+j],ctx))nz[d][i].push_back(j);}}
 vector<vector<int>>table(npairs);for(int z=0;z<npairs;z++){int da=z,db=nw-1-z;table[z].resize(sizes[da]*sizes[db]);for(auto&i:table[z])in>>i;}
 vector<vector<int>>nfrows(fullcols);vector<Array*>nfcoeff(fullcols,nullptr);for(int i=0;i<fullcols;i++){int count;in>>count;nfrows[i].resize(count);if(count){nfcoeff[i]=new Array(count,ctx);for(int j=0;j<count;j++){in>>nfrows[i][j];readpoly(in,(*nfcoeff[i])[j],ctx,n);}}}
 int pairs;in>>pairs;Array result(eq*cols,ctx),z(cols,ctx),temp(2,ctx);
 auto convolution=[&](int da,int db,int i,int j,int ti){for(int a:nz[da][i])for(int b:nz[db][j]){int k=table[ti][a*sizes[db]+b];if(nfrows[k].empty())continue;fmpq_mpoly_mul(temp[0],(*W[da])[i*sizes[da]+a],(*W[db])[j*sizes[db]+b],ctx);for(size_t q=0;q<nfrows[k].size();q++){fmpq_mpoly_mul(temp[1],temp[0],(*nfcoeff[k])[q],ctx);int dest=nfrows[k][q];fmpq_mpoly_add(z[dest],z[dest],temp[1],ctx);}}};
 for(int pair=0;pair<pairs;pair++){
  int i,j,count;in>>i>>j>>count;vector<int>rows(count);Array coeff(count,ctx);for(int k=0;k<count;k++){in>>rows[k];readpoly(in,coeff[k],ctx,n);}for(int k=0;k<cols;k++)fmpq_mpoly_zero(z[k],ctx);
  for(int z=0;z<npairs;z++){int da=z,db=nw-1-z;convolution(da,db,i,j,z);if(da!=db)convolution(da,db,j,i,z);}
  vector<int>support;for(int k=0;k<cols;k++)if(!fmpq_mpoly_is_zero(z[k],ctx))support.push_back(k);
  for(int r=0;r<count;r++)for(int k:support){fmpq_mpoly_mul(temp[0],coeff[r],z[k],ctx);fmpq_mpoly_add(result[rows[r]*cols+k],result[rows[r]*cols+k],temp[0],ctx);}
  if(pair%100==0)cerr<<"pair "<<pair<<"/"<<pairs<<"\n";
 }
 if(!in)throw runtime_error("truncated input");ofstream out(argv[2]);if(!out)throw runtime_error("output open failed");out<<eq<<' '<<cols<<' '<<n<<'\n';long total=0;
 for(int r=0;r<eq;r++){vector<int>support;for(int j=0;j<cols;j++)if(!fmpq_mpoly_is_zero(result[r*cols+j],ctx))support.push_back(j);out<<support.size()<<'\n';for(int j:support){out<<j<<'\n';writepoly(out,result[r*cols+j],ctx,n);total+=fmpq_mpoly_length(result[r*cols+j],ctx);}}
 out.flush();if(!out)throw runtime_error("output write failed");cerr<<"terms "<<total<<"\n";for(auto*w:W)delete w;for(auto*x:nfcoeff)delete x;
 }
 fmpq_mpoly_ctx_clear(ctx);flint_cleanup();return 0;
}catch(const exception&e){cerr<<e.what()<<"\n";return 1;}}
