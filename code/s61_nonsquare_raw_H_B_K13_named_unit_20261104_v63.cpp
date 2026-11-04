// Fixed degree<=13 coefficient-unit witness, used inside one named-root job.
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <string>
using namespace std;
struct P{fmpz_poly_t v;P(){fmpz_poly_init(v);}~P(){fmpz_poly_clear(v);}};
void readp(istream&f,P&p){long n;f>>n;assert(n>=0&&n<=14);fmpz_t z;fmpz_init(z);for(long j=0;j<n;j++){long a;string s;f>>a>>s;assert(f&&0<=a&&a<=13&&!fmpz_set_str(z,s.c_str(),10));fmpz_poly_set_coeff_fmpz(p.v,a,z);}fmpz_clear(z);}
void jp(ostream&o,const P&p){o<<'[';bool first=true;fmpz_t z;fmpz_init(z);for(long i=fmpz_poly_degree(p.v);i>=0;i--){fmpz_poly_get_coeff_fmpz(z,p.v,i);if(fmpz_is_zero(z))continue;if(!first)o<<',';first=false;char*s=fmpz_get_str(nullptr,10,z);o<<'['<<i<<",\""<<s<<"\"]";flint_free(s);assert(o.tellp()<28*1024*1024);}o<<']';fmpz_clear(z);}
int main(int argc,char**argv){
 assert(argc==3);ifstream in(argv[1]);P g,p,d,S,T,x,y;readp(in,g);readp(in,p);assert(fmpz_poly_degree(g.v)==13&&!fmpz_poly_is_zero(p.v));fmpz_poly_gcd_subresultant(d.v,g.v,p.v);ofstream out(argv[2]);out<<"{\"gcd\":";jp(out,d);
 if(fmpz_poly_degree(d.v)==0){fmpz_t R;fmpz_init(R);fmpz_poly_xgcd_modular(R,T.v,S.v,g.v,p.v);assert(!fmpz_is_zero(R));fmpz_poly_mul(x.v,S.v,p.v);fmpz_poly_mul(y.v,T.v,g.v);fmpz_poly_add(x.v,x.v,y.v);fmpz_poly_set_fmpz(y.v,R);assert(fmpz_poly_equal(x.v,y.v));out<<",\"S\":";jp(out,S);out<<",\"T\":";jp(out,T);char*s=fmpz_get_str(nullptr,10,R);out<<",\"R\":\""<<s<<"\",\"unit\":true";flint_free(s);fmpz_clear(R);}else out<<",\"unit\":false";
 out<<"}\n";
}
