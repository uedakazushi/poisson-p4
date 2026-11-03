// Independent INTEGER verification of saved QQ cross inverse and reductions.
// CLI: cross-input saved-cross-output report-json. No inverse is recomputed.
// Denominators are cleared globally; all checks use fmpz_mat multiplication,
// whereas the producer computes its inverse/products with fmpq_mat.
#include <flint/fmpq_mat.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpq.h>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
void readq(istream&in,fmpq_t q,fmpz_t den){string s;in>>s;assert(in&&!fmpq_set_str(q,s.c_str(),10));fmpq_canonicalise(q);assert(fmpz_fdiv_ui(fmpq_denref(q),32003));fmpz_lcm(den,den,fmpq_denref(q));}
void integerize(fmpz_mat_t Z,const fmpq_mat_t Q,const fmpz_t den){fmpz_t t;fmpz_init(t);for(slong i=0;i<fmpq_mat_nrows(Q);i++)for(slong j=0;j<fmpq_mat_ncols(Q);j++){auto q=fmpq_mat_entry(Q,i,j);fmpz_divexact(t,den,fmpq_denref(q));fmpz_mul(fmpz_mat_entry(Z,i,j),t,fmpq_numref(q));}fmpz_clear(t);}
bool scalaridentity(const fmpz_mat_t A,const fmpz_t d){for(slong i=0;i<fmpz_mat_nrows(A);i++)for(slong j=0;j<fmpz_mat_ncols(A);j++)if(i==j){if(!fmpz_equal(fmpz_mat_entry(A,i,j),d))return false;}else if(!fmpz_is_zero(fmpz_mat_entry(A,i,j)))return false;return true;}
int main(int argc,char**argv){assert(argc==4);ifstream in(argv[1]),saved(argv[2]);ofstream report(argv[3]);int nr,nc,nd;in>>nr>>nc>>nd;assert(nr==129&&nc==126&&nd==945);int n,m;saved>>n>>m;assert(n==126&&m==945);vector<int>ids(126);vector<bool>seen(129);for(int&i:ids){saved>>i;assert(i>=0&&i<129&&!seen[i]);seen[i]=true;}
 fmpq_mat_t C,D,S,T,V,P;fmpq_mat_init(C,129,126);fmpq_mat_init(D,129,945);fmpq_mat_init(S,126,126);fmpq_mat_init(T,126,945);fmpq_mat_init(V,126,126);fmpq_mat_init(P,126,945);fmpz_t delta,lambda,mu,diag;fmpz_init(delta);fmpz_init(lambda);fmpz_init(mu);fmpz_init(diag);fmpz_one(delta);fmpz_one(lambda);fmpz_one(mu);
 for(int i=0;i<129;i++)for(int j=0;j<126;j++)readq(in,fmpq_mat_entry(C,i,j),delta);for(int i=0;i<129;i++)for(int j=0;j<945;j++)readq(in,fmpq_mat_entry(D,i,j),delta);
 for(int i=0;i<126;i++){for(int j=0;j<126;j++)fmpq_set(fmpq_mat_entry(S,i,j),fmpq_mat_entry(C,ids[i],j));for(int j=0;j<945;j++)fmpq_set(fmpq_mat_entry(T,i,j),fmpq_mat_entry(D,ids[i],j));}
 for(int i=0;i<126;i++)for(int j=0;j<126;j++)readq(saved,fmpq_mat_entry(V,i,j),lambda);for(int i=0;i<126;i++)for(int j=0;j<945;j++)readq(saved,fmpq_mat_entry(P,i,j),mu);string extra;assert(!(in>>extra)&&!(saved>>extra));assert(fmpz_fdiv_ui(delta,32003)&&fmpz_fdiv_ui(lambda,32003)&&fmpz_fdiv_ui(mu,32003));cerr<<"ALL_INPUT_AND_SAVED_DENOMINATORS_PUNIT\n";
 fmpz_mat_t ZS,ZT,ZV,ZP,A,B;fmpz_mat_init(ZS,126,126);fmpz_mat_init(ZT,126,945);fmpz_mat_init(ZV,126,126);fmpz_mat_init(ZP,126,945);fmpz_mat_init(A,126,126);fmpz_mat_init(B,126,945);integerize(ZS,S,delta);integerize(ZT,T,delta);integerize(ZV,V,lambda);integerize(ZP,P,mu);fmpz_mul(diag,delta,lambda);
 fmpz_mat_mul(A,ZV,ZS);assert(scalaridentity(A,diag));fmpz_mat_mul(A,ZS,ZV);assert(scalaridentity(A,diag));cerr<<"BOTH_126_INTEGER_INVERSE_IDENTITIES_PASS\n";
 fmpz_mat_mul(B,ZS,ZP);fmpz_mat_scalar_mul_fmpz(ZT,ZT,mu);fmpz_mat_add(B,B,ZT);assert(fmpz_mat_is_zero(B));cerr<<"ALL_126x945_INTEGER_CROSS_REDUCTION_IDENTITIES_PASS\n";
 report<<"{\"all_passed\":true,\"rank\":126,\"inverse_products_checked_over_Z\":2,\"reduction_coefficients_checked_over_Z\":119070,\"all_denominators_32003_units\":true,\"input_denominator_bits\":"<<fmpz_bits(delta)<<",\"inverse_denominator_bits\":"<<fmpz_bits(lambda)<<",\"reduction_denominator_bits\":"<<fmpz_bits(mu)<<"}\n";
 fmpq_mat_clear(C);fmpq_mat_clear(D);fmpq_mat_clear(S);fmpq_mat_clear(T);fmpq_mat_clear(V);fmpq_mat_clear(P);fmpz_mat_clear(ZS);fmpz_mat_clear(ZT);fmpz_mat_clear(ZV);fmpz_mat_clear(ZP);fmpz_mat_clear(A);fmpz_mat_clear(B);fmpz_clear(delta);fmpz_clear(lambda);fmpz_clear(mu);fmpz_clear(diag);flint_cleanup();}
