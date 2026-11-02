#include "s61_binary_univariate_field_20261102.hpp"
#include <chrono>
int main(int argc,char**argv){
 assert(argc==3);ifstream in(argv[1]),ki(argv[2]);int rows,n;in>>rows>>n;
 vector<URow>A(rows);for(auto&r:A){int count;in>>count;for(int k=0;k<count;k++){int j;URat a;in>>j;uread(in,a.n);r.emplace(j,std::move(a));}}
 int kn,degree,number;ki>>kn>>degree>>number;assert(kn==n);fmpq_t c;fmpq_init(c);auto start=chrono::steady_clock::now();
 for(int a=0;a<number;a++){int count;ki>>count;vector<URat>v(n);for(int k=0;k<count;k++){int j;string s;ki>>j>>s;assert(fmpq_set_str(c,s.c_str(),10)==0);fmpq_canonicalise(c);fmpq_poly_set_coeff_fmpq(v[j%n].n,j/n,c);}URat sum,tmp;for(int i=0;i<rows;i++){fmpq_poly_zero(sum.n);for(auto&[j,q]:A[i])if(!v[j].zero()){fmpq_poly_mul(tmp.n,q.n,v[j].n);fmpq_poly_add(sum.n,sum.n,tmp.n);}assert(sum.zero());}}
 cerr<<"EXACT_KERNELS "<<number<<" ROWS "<<rows<<" IDENTITIES "<<number*rows<<" SECONDS "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';fmpq_clear(c);flint_cleanup();}
