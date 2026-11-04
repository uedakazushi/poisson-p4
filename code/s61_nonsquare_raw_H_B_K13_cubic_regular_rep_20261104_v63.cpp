// One fixed-prime regular representation. One full work matrix, freed before
// allocating the ONE combined-LU minor. Every selected PLU entry is checked.
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
using namespace std;using Z=uint32_t;const Z P=32003;const int D=13,R=221,C=312;
Z mul(Z a,Z b){return uint64_t(a)*b%P;}Z sub(Z a,Z b){return(a+P-b)%P;}
Z inv(Z a){assert(a);Z r=1;for(Z n=P-2;n;n>>=1,a=mul(a,a))if(n&1)r=mul(r,a);return r;}
void arr(ostream&o,const vector<Z>&v){o<<'[';for(size_t i=0;i<v.size();i++){if(i)o<<',';o<<v[i];}o<<']';}
void iarr(ostream&o,const vector<int>&v){o<<'[';for(size_t i=0;i<v.size();i++){if(i)o<<',';o<<v[i];}o<<']';}
int main(int argc,char**argv){
 assert(argc==3);ifstream in(argv[1]);int pp,nr,nc,dd;in>>pp>>nr>>nc>>dd;assert(pp==P&&nr==17&&nc==24&&dd==D);
 vector<Z>g(D),f(17*24*D);for(auto&x:g){in>>x;assert(x<P);}for(auto&x:f){in>>x;assert(x<P);}string extra;assert(!(in>>extra));
 auto shift=[&](vector<Z>&v){Z top=v[D-1];for(int k=D-1;k>=1;k--)v[k]=sub(v[k-1],mul(top,g[k]));v[0]=sub(0,mul(top,g[0]));};
 auto column=[&](int blockrow,int fullcolumn){int b=fullcolumn/D,e=fullcolumn%D;vector<Z>v(D);for(int k=0;k<D;k++)v[k]=f[(blockrow*24+b)*D+k];for(int j=0;j<e;j++)shift(v);return v;};
 vector<Z>M(R*C);for(int i=0;i<17;i++)for(int j=0;j<C;j++){auto v=column(i,j);for(int k=0;k<D;k++)M[(i*D+k)*C+j]=v[k];}
 vector<int>cols;int rank=0;
 for(int j=0;j<C&&rank<R;j++){
  int p=rank;while(p<R&&!M[p*C+j])p++;if(p==R)continue;
  if(p!=rank)for(int k=0;k<C;k++)swap(M[p*C+k],M[rank*C+k]);Z pinv=inv(M[rank*C+j]);
  for(int i=rank+1;i<R;i++)if(M[i*C+j]){Z q=mul(M[i*C+j],pinv);M[i*C+j]=0;for(int k=j+1;k<C;k++)M[i*C+k]=sub(M[i*C+k],mul(q,M[rank*C+k]));}
  cols.push_back(j);rank++;
 }
 // No second full representation or selected-minor copy exists.
 vector<Z>().swap(M);
 ofstream out(argv[2]);out<<"{\"prime\":32003,\"regular_shape\":[221,312],\"rank_modp\":"<<rank<<",\"selected_columns\":";iarr(out,cols);
 if(rank==R){
  vector<Z>A(R*R);for(int i=0;i<17;i++)for(int j=0;j<R;j++){auto v=column(i,cols[j]);for(int k=0;k<D;k++)A[(i*D+k)*R+j]=v[k];}
  vector<int>perm(R);iota(perm.begin(),perm.end(),0);int swaps=0;Z det=1;vector<Z>diag;
  for(int k=0;k<R;k++){
   int p=k;while(p<R&&!A[p*R+k])p++;assert(p<R);
   if(p!=k){for(int j=0;j<R;j++)swap(A[p*R+j],A[k*R+j]);swap(perm[p],perm[k]);swaps++;}
   Z pivot=A[k*R+k];diag.push_back(pivot);det=mul(det,pivot);Z pinv=inv(pivot);
   for(int i=k+1;i<R;i++){Z q=mul(A[i*R+k],pinv);A[i*R+k]=q;for(int j=k+1;j<R;j++)A[i*R+j]=sub(A[i*R+j],mul(q,A[k*R+j]));}
  }
  if(swaps%2)det=sub(0,det);assert(det);
  long checks=0;for(int i=0;i<R;i++)for(int j=0;j<R;j++){
   Z s=0;for(int k=0;k<=min(i,j);k++){Z l=k==i?1:A[i*R+k];s=(s+mul(l,A[k*R+j]))%P;}
   auto original=column(perm[i]/D,cols[j]);assert(s==original[perm[i]%D]);checks++;
  }
  out<<",\"actual_selected_minor_determinant\":"<<det<<",\"row_swaps\":"<<swaps<<",\"row_permutation_PA_equals_LU\":";iarr(out,perm);
  out<<",\"diagonal_U\":";arr(out,diag);out<<",\"combined_LU221x221\":[";
  for(int i=0;i<R;i++){if(i)out<<',';out<<'[';for(int j=0;j<R;j++){if(j)out<<',';out<<A[i*R+j];}out<<']';}out<<']';
  out<<",\"all_selected_PLU_entry_checks\":"<<checks<<",\"full_work_matrix_freed_before_LU\":true";
 }
 out<<",\"peak_native_scalar_slots_upper\":"<<(R*C+int(f.size())+D)<<",\"positive_only\":true}\n";assert(out.tellp()<30*1024*1024);
}
