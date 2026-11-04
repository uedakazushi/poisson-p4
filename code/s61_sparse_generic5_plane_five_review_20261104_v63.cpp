// Independent selected-five-row source substitution and polynomial-matrix det.
// No producer solver, row selection, point, prime, or field elimination.
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <gmp.h>
#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly_mat.h>
using namespace std;
struct P{fmpq_poly_t x;P(){fmpq_poly_init(x);}P(const P&a):P(){fmpq_poly_set(x,a.x);}P(P&&a):P(){fmpq_poly_swap(x,a.x);}P&operator=(P a){fmpq_poly_swap(x,a.x);return *this;}~P(){fmpq_poly_clear(x);}};
P readp(istream&in){P p;long n;in>>n;assert(n>=0);fmpq_t q;fmpq_init(q);for(long i=0;i<n;i++){long e;string c;in>>e>>c;assert(e>=0&&!fmpq_set_str(q,c.c_str(),10));fmpq_canonicalise(q);fmpq_poly_set_coeff_fmpq(p.x,e,q);}fmpq_clear(q);return p;}
P add(const P&a,const P&b){P c;fmpq_poly_add(c.x,a.x,b.x);return c;}
P mul(const P&a,const P&b){P c;fmpq_poly_mul(c.x,a.x,b.x);return c;}
void eq(const P&a,const P&b){assert(fmpq_poly_equal(a.x,b.x));}
int main(int argc,char**argv){assert(argc==2);auto start=chrono::steady_clock::now();ifstream in(argv[1]);P common=readp(in);vector<P>C(45),D(45);vector<vector<P>>W(45,vector<P>(2));for(auto&p:C)p=readp(in);for(auto&p:D)p=readp(in);for(auto&r:W)for(auto&p:r)p=readp(in);
  vector<vector<P>>matrix(5,vector<P>(5));long source_terms=0,degree=0;
  for(int row=0;row<5;row++){int n;in>>n;for(int j=0;j<n;j++){int qm;in>>qm;assert(qm>=0&&qm<=2);vector<int>m(qm);for(auto&i:m){in>>i;assert(i>=0&&i<45);}P co=readp(in);source_terms+=fmpq_poly_length(co.x);
    if(qm==0)matrix[row][0]=add(matrix[row][0],co);
    else if(qm==1){int a=m[0];matrix[row][0]=add(matrix[row][0],mul(co,C[a]));for(int k=0;k<2;k++)matrix[row][1+k]=add(matrix[row][1+k],mul(co,W[a][k]));}
    else {int a=m[0],b=m[1];matrix[row][0]=add(matrix[row][0],mul(co,mul(C[a],C[b])));for(int k=0;k<2;k++){
      matrix[row][1+k]=add(matrix[row][1+k],mul(co,add(mul(W[a][k],C[b]),mul(C[a],W[b][k]))));
      matrix[row][3+k]=add(matrix[row][3+k],mul(co,add(mul(W[a][k],D[b]),mul(D[a],W[b][k]))));
    }}
  }matrix[row][0]=mul(common,matrix[row][0]);}
  for(int i=0;i<5;i++)for(int j=0;j<5;j++){eq(matrix[i][j],readp(in));degree=max(degree,fmpq_poly_degree(matrix[i][j].x));}
  P expected=readp(in);string extra;assert(in&&!(in>>extra));assert(!fmpq_poly_is_zero(expected.x));
  // Clear one integer denominator per row. FLINT's exact polynomial matrix
  // determinant is independent of the producer's saved field/Bareiss trace.
  fmpz_poly_mat_t M;fmpz_poly_mat_init(M,5,5);fmpz_t scale,all,ratio;fmpz_init(scale);fmpz_init(all);fmpz_init(ratio);fmpz_one(all);
  for(int i=0;i<5;i++){fmpz_one(scale);for(int j=0;j<5;j++)fmpz_lcm(scale,scale,fmpq_poly_denref(matrix[i][j].x));fmpz_mul(all,all,scale);
    for(int j=0;j<5;j++){fmpz_divexact(ratio,scale,fmpq_poly_denref(matrix[i][j].x));fmpq_poly_get_numerator(fmpz_poly_mat_entry(M,i,j),matrix[i][j].x);fmpz_poly_scalar_mul_fmpz(fmpz_poly_mat_entry(M,i,j),fmpz_poly_mat_entry(M,i,j),ratio);}}
  fmpz_poly_t determinant;fmpz_poly_init(determinant);fmpz_poly_mat_det(determinant,M);assert(!fmpz_poly_is_zero(determinant));P actual;fmpq_poly_set_fmpz_poly(actual.x,determinant);fmpq_poly_scalar_div_fmpz(actual.x,actual.x,all);eq(actual,expected);
  cout<<"{\"all_passed\":true,\"selected_rows\":5,\"actual_entries_verified\":25,\"source_coefficient_terms\":"<<source_terms<<",\"max_entry_degree\":"<<degree<<",\"determinant_degree\":"<<fmpq_poly_degree(actual.x)<<",\"nonzero_QQ_polynomial_determinant\":true,\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
  fmpz_poly_clear(determinant);fmpz_poly_mat_clear(M);fmpz_clear(scale);fmpz_clear(all);fmpz_clear(ratio);
}
