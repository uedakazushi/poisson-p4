#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
using namespace std;
int main(int argc,char**argv){assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);int m,n,k;in>>m>>n>>k;assert(in&&out);fmpq_mat_t A,R;fmpq_mat_init(A,m,n+k);fmpq_mat_init(R,m,n+k);for(int i=0;i<m;i++)for(int j=0;j<n+k;j++){string v;in>>v;assert(fmpq_set_str(fmpq_mat_entry(A,i,j),v.c_str(),10)==0);fmpq_canonicalise(fmpq_mat_entry(A,i,j));}slong rank=fmpq_mat_rref(R,A);int consistent=1;vector<int>piv(n,-1);for(int i=0;i<rank;i++){int j=0;while(j<n+k&&fmpq_is_zero(fmpq_mat_entry(R,i,j)))j++;assert(j<n+k);if(j>=n)consistent=0;else piv[j]=i;}out<<n<<' '<<k<<' '<<rank<<' '<<consistent<<'\n';if(consistent)for(int i=0;i<n;i++)for(int j=0;j<k;j++){if(piv[i]<0)out<<"0\n";else{char*v=fmpq_get_str(nullptr,10,fmpq_mat_entry(R,piv[i],n+j));out<<v<<'\n';flint_free(v);}}assert(in&&out);}
