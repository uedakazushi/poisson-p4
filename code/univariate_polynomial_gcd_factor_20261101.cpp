// Exact factorization of the previously certified gcd.
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
using namespace std;
int main(int argc,char**argv){assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);int k;in>>k;fmpz_poly_t g;fmpz_poly_init(g);fmpz_t c;fmpz_init(c);for(int a=0;a<=k;a++){int n;in>>n;fmpz_poly_zero(g);for(int i=0;i<n;i++){string v;in>>v;if(a==k){assert(fmpz_set_str(c,v.c_str(),10)==0);fmpz_poly_set_coeff_fmpz(g,i,c);}}}fmpz_poly_factor_t f;fmpz_poly_factor_init(f);fmpz_poly_factor(f,g);out<<f->num<<'\n';for(int i=0;i<f->num;i++){auto p=f->p+i;out<<f->exp[i]<<' '<<fmpz_poly_length(p)<<'\n';for(int j=0;j<fmpz_poly_length(p);j++){char*v=fmpz_get_str(nullptr,10,fmpz_poly_get_coeff_ptr(p,j));out<<v<<'\n';flint_free(v);}cerr<<"factor degree "<<fmpz_poly_degree(p)<<" exponent "<<f->exp[i]<<'\n';}assert(in&&out);}
