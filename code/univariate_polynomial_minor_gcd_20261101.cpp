// Exact selected maximal minors and their polynomial gcd.
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_mat.h>
using namespace std;
void write(ostream&out,fmpz_poly_t p){out<<fmpz_poly_length(p)<<'\n';for(int i=0;i<fmpz_poly_length(p);i++){char*c=fmpz_get_str(nullptr,10,fmpz_poly_get_coeff_ptr(p,i));out<<c<<'\n';flint_free(c);}}
int main(int argc,char**argv){assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);assert(in&&out);int m,n,k,rank;in>>m>>n>>k>>rank;fmpz_poly_mat_t M,N;fmpz_poly_mat_init(M,m,n);fmpz_poly_mat_init(N,rank,rank);fmpz_t c;fmpz_init(c);for(int i=0;i<m;i++)for(int j=0;j<n;j++){int t;in>>t;while(t--){int e;string v;in>>e>>v;assert(fmpz_set_str(c,v.c_str(),10)==0);fmpz_poly_set_coeff_fmpz(fmpz_poly_mat_entry(M,i,j),e,c);}}
fmpz_poly_t d,g,p;fmpz_poly_init(d);fmpz_poly_init(g);fmpz_poly_init(p);auto st=chrono::steady_clock::now();out<<k<<'\n';for(int a=0;a<k;a++){vector<int>rr(rank),cc(rank);for(auto&i:rr)in>>i;for(auto&i:cc)in>>i;for(int i=0;i<rank;i++)for(int j=0;j<rank;j++)fmpz_poly_set(fmpz_poly_mat_entry(N,i,j),fmpz_poly_mat_entry(M,rr[i],cc[j]));fmpz_poly_mat_det(d,N);assert(!fmpz_poly_is_zero(d));fmpz_poly_primitive_part(p,d);if(!a)fmpz_poly_set(g,p);else fmpz_poly_gcd(g,g,p);write(out,p);cerr<<"minor "<<a<<" degree "<<fmpz_poly_degree(p)<<" gcd degree "<<fmpz_poly_degree(g)<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-st).count()<<'\n';}write(out,g);assert(in&&out);}
