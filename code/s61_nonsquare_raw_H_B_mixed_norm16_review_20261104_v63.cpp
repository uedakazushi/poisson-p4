// Independentlibrarydet ONLY. Author usedscalarBareiss; no solve/rank/point/GCD here.
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_mat.h>
#include <cassert>
#include <chrono>
#include <fstream>
#include <string>
using namespace std;
void readpoly(istream&f,fmpz_poly_t p){long n;f>>n;assert(f&&n>=0);fmpz_t c;fmpz_init(c);for(long i=0;i<n;i++){long j;string s;f>>j>>s;assert(j>=0&&!fmpz_set_str(c,s.c_str(),10));fmpz_poly_set_coeff_fmpz(p,j,c);}fmpz_clear(c);}
int main(int argc,char**argv){assert(argc==3);auto start=chrono::steady_clock::now();ifstream f(argv[1]);fmpz_poly_mat_t A;fmpz_poly_mat_init(A,16,16);for(int i=0;i<16;i++)for(int j=0;j<16;j++)readpoly(f,fmpz_poly_mat_entry(A,i,j));fmpz_poly_t expected,det;fmpz_poly_init(expected);fmpz_poly_init(det);readpoly(f,expected);string extra;assert(!(f>>extra));fmpz_poly_mat_det_fflu(det,A);assert(!fmpz_poly_is_zero(det)&&fmpz_poly_equal(det,expected));long lo=0;fmpz_t c;fmpz_init(c);while(true){fmpz_poly_get_coeff_fmpz(c,det,lo);if(!fmpz_is_zero(c))break;lo++;}ofstream o(argv[2]);o<<"{\"algorithm\":\"FLINT_fmpz_poly_mat_det_fflu\",\"independent_exact_det_equals_author\":true,\"degree\":"<<fmpz_poly_degree(det)<<",\"a_valuation\":"<<lo<<",\"determinant\":[";bool first=true;for(long j=fmpz_poly_degree(det);j>=0;j--){fmpz_poly_get_coeff_fmpz(c,det,j);if(fmpz_is_zero(c))continue;if(!first)o<<',';first=false;char*s=fmpz_get_str(nullptr,10,c);o<<'['<<j<<",\""<<s<<"\"]";flint_free(s);}o<<"],\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";assert(o.tellp()<32*1024*1024);fmpz_clear(c);fmpz_poly_clear(expected);fmpz_poly_clear(det);fmpz_poly_mat_clear(A);flint_cleanup();}
