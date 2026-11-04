// Independent ONE 15x15 cofactor00 Bareiss; no matrix-det/solve API.
#include <flint/fmpz_poly.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
struct P {
 fmpz_poly_t v;
 P(){fmpz_poly_init(v);} P(const P&x):P(){fmpz_poly_set(v,x.v);}
 P(P&&x):P(){fmpz_poly_swap(v,x.v);} ~P(){fmpz_poly_clear(v);}
 P&operator=(P x){fmpz_poly_swap(v,x.v);return *this;}
};
void readp(istream&in,P&p){
 long n;in>>n;assert(in&&n>=0);fmpz_t z;fmpz_init(z);
 vector<long>seen;
 for(long j=0;j<n;j++){long degree;string c;in>>degree>>c;assert(in&&degree>=0&&find(seen.begin(),seen.end(),degree)==seen.end());seen.push_back(degree);assert(!fmpz_set_str(z,c.c_str(),10)&&!fmpz_is_zero(z));fmpz_poly_set_coeff_fmpz(p.v,degree,z);}
 assert(fmpz_poly_length(p.v)==0||fmpz_poly_degree(p.v)>=0);fmpz_clear(z);
}
int main(int argc,char**argv){
 assert(argc==3);auto start=chrono::steady_clock::now();ifstream in(argv[1]);
 vector<vector<P>>N(16,vector<P>(16));for(auto&r:N)for(auto&p:r)readp(in,p);P expected;readp(in,expected);string trailing;assert(!(in>>trailing)&&!fmpz_poly_is_zero(expected.v));
 vector<vector<P>>A(15,vector<P>(15));for(int i=0;i<15;i++)for(int j=0;j<15;j++)A[i][j]=N[i+1][j+1];
 P previous;fmpz_poly_one(previous.v);int sign=1;long divisions=0,maxdegree=0,maxbits=0;vector<long>pivotrows,pivotdegrees;
 for(int k=0;k<14;k++){
  assert(chrono::duration<double>(chrono::steady_clock::now()-start).count()<52);
  // Last nonzero pivot is deliberately different from first-pivot FFLU.
  int row=14;while(row>=k&&fmpz_poly_is_zero(A[row][k].v))row--;assert(row>=k);
  if(row!=k){swap(A[row],A[k]);sign=-sign;}pivotrows.push_back(row);P pivot=A[k][k];pivotdegrees.push_back(fmpz_poly_degree(pivot.v));
  for(int i=k+1;i<15;i++)for(int j=k+1;j<15;j++){
   P numerator,other,quotient,recheck;
   fmpz_poly_mul(numerator.v,pivot.v,A[i][j].v);fmpz_poly_mul(other.v,A[i][k].v,A[k][j].v);fmpz_poly_sub(numerator.v,numerator.v,other.v);
   assert(fmpz_poly_divides(quotient.v,numerator.v,previous.v));fmpz_poly_mul(recheck.v,quotient.v,previous.v);assert(fmpz_poly_equal(recheck.v,numerator.v));divisions++;
   maxdegree=max(maxdegree,(long)fmpz_poly_degree(numerator.v));maxbits=max(maxbits,labs((long)fmpz_poly_max_bits(numerator.v)));A[i][j]=std::move(quotient);
  }
  for(int i=k+1;i<15;i++)fmpz_poly_zero(A[i][k].v);previous=std::move(pivot);
 }
 P cofactor=A[14][14];if(sign<0)fmpz_poly_neg(cofactor.v,cofactor.v);assert(fmpz_poly_equal(cofactor.v,expected.v));
 ofstream out(argv[2]);out<<"{\"cofactor00_equals_saved_V0\":true,\"all_Bareiss_divisions_reverse_multiplied\":true,\"exact_divisions\":"<<divisions<<",\"maximum_numerator_degree\":"<<maxdegree<<",\"maximum_numerator_coefficient_bits\":"<<maxbits<<",\"pivot_rows_last_nonzero\":[";
 for(size_t i=0;i<pivotrows.size();i++){if(i)out<<',';out<<pivotrows[i];}out<<"],\"pivot_degrees\":[";for(size_t i=0;i<pivotdegrees.size();i++){if(i)out<<',';out<<pivotdegrees[i];}
 out<<"],\"cofactor00\":[";bool first=true;fmpz_t z;fmpz_init(z);for(long i=fmpz_poly_degree(cofactor.v);i>=0;i--){fmpz_poly_get_coeff_fmpz(z,cofactor.v,i);if(fmpz_is_zero(z))continue;if(!first)out<<',';first=false;char*s=fmpz_get_str(nullptr,10,z);out<<'['<<i<<",\""<<s<<"\"]";flint_free(s);}fmpz_clear(z);
 out<<"],\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";assert(out.tellp()<32*1024*1024);cerr<<"ONE_COFACTOR00_BAREISS_ALL_DIVISIONS_AND_FINAL_V0_PASS\n";flint_cleanup();
}
