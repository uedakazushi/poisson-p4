// ONE exact integer multivariate s-resultant; bounded by the external120s/2GB job.
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <openssl/evp.h>
using namespace std;
fmpz_mpoly_ctx_t ctx;
string zstr(const fmpz_t z){char*p=fmpz_get_str(nullptr,10,z);string s(p);flint_free(p);return s;}
string termline(const fmpz_mpoly_t a,long k){fmpz_t z;fmpz_init(z);ulong e[6];fmpz_mpoly_get_term_coeff_fmpz(z,a,k,ctx);fmpz_mpoly_get_term_exp_ui(e,a,k,ctx);ostringstream o;o<<zstr(z);for(auto x:e)o<<' '<<x;o<<'\n';fmpz_clear(z);return o.str();}
void readpoly(istream&f,fmpz_mpoly_t a){long n;f>>n;assert(n>0);fmpz_t z;fmpz_init(z);ulong e[6];for(long k=0;k<n;k++){string c;f>>c;for(auto&x:e)f>>x;assert(f&&fmpz_set_str(z,c.c_str(),10)==0);fmpz_mpoly_push_term_fmpz_ui(a,z,e,ctx);}fmpz_mpoly_sort_terms(a,ctx);fmpz_mpoly_combine_like_terms(a,ctx);fmpz_clear(z);assert(fmpz_mpoly_length(a,ctx)==n);}
void front(const string&dir,const string&stage,long n=-1){ofstream f(dir+"/native-frontier.json");f<<"{\"stage\":\""<<stage<<"\",\"resultant_call_count\":"<<(stage=="input_loaded"?0:1)<<",\"resultant_terms_if_returned\":"<<n<<"}\n";}
int main(int argc,char**argv){assert(argc==4);fmpz_mpoly_ctx_init(ctx,6,ORD_DEGLEX);auto start=chrono::steady_clock::now();int result=0;{
 fmpz_mpoly_t a,b,res;fmpz_mpoly_init(a,ctx);fmpz_mpoly_init(b,ctx);fmpz_mpoly_init(res,ctx);
 try{
  ifstream in(argv[1]);int vars,elim;in>>vars>>elim;assert(vars==6&&elim==0);readpoly(in,a);readpoly(in,b);string extra;assert(!(in>>extra));
  long budget=stol(argv[3]);assert(budget>0&&budget<=32L*1024*1024);string dir=argv[2];front(dir,"input_loaded");
  long da=fmpz_mpoly_degree_si(a,0,ctx),db=fmpz_mpoly_degree_si(b,0,ctx);assert(da>0&&db>0);cerr<<"INPUT_S_DEGREES "<<da<<' '<<db<<" TERMS "<<fmpz_mpoly_length(a,ctx)<<' '<<fmpz_mpoly_length(b,ctx)<<endl;
  front(dir,"calling_single_fmpz_mpoly_resultant");int success=fmpz_mpoly_resultant(res,a,b,0,ctx);
  long n=fmpz_mpoly_length(res,ctx);front(dir,"library_resultant_returned",n);bool zero=fmpz_mpoly_is_zero(res,ctx);
  slong degrees[6];fmpz_mpoly_degrees_si(degrees,res,ctx);if(!zero)assert(degrees[0]==0);
  EVP_MD_CTX*digest=EVP_MD_CTX_new();assert(digest&&EVP_DigestInit_ex(digest,EVP_sha256(),nullptr));uint64_t bytes=0;long bits=0;fmpz_t z;fmpz_init(z);
  for(long k=0;k<n;k++){string line=termline(res,k);assert(EVP_DigestUpdate(digest,line.data(),line.size()));bytes+=line.size();fmpz_mpoly_get_term_coeff_fmpz(z,res,k,ctx);bits=max(bits,long(fmpz_bits(z)));}
  unsigned char hash[EVP_MAX_MD_SIZE];unsigned int len;assert(EVP_DigestFinal_ex(digest,hash,&len));EVP_MD_CTX_free(digest);assert(len==32);ostringstream hex;for(unsigned int i=0;i<len;i++)hex<<std::hex<<setw(2)<<setfill('0')<<unsigned(hash[i]);fmpz_clear(z);
  string header="6 "+to_string(n)+"\n";bool written=success&&!zero&&bytes+header.size()<=uint64_t(budget);
  if(written){ofstream out(dir+"/resultant.txt");out<<header;for(long k=0;k<n;k++)out<<termline(res,k);out.close();}
  string status=!success?"FLINT_resultant_returned_failure":zero?"zero_resultant_structural_frontier":!written?"resultant_exceeds_artifact_cap_hash_only":"nonzero_resultant_saved_pending_independent_identity_audit";
  ofstream out(dir+"/native-report.json");out<<"{\"all_passed\":"<<(success?"true":"false")<<",\"status\":\""<<status<<"\",\"FLINT_resultant_return_code\":"<<success<<",\"resultant_call_count\":1,\"zero\":"<<(zero?"true":"false")<<",\"terms\":"<<n<<",\"multidegree\":[";for(int i=0;i<6;i++){if(i)out<<',';out<<degrees[i];}out<<"],\"max_coefficient_bits\":"<<bits<<",\"canonical_term_encoding_bytes\":"<<bytes<<",\"canonical_term_SHA256\":\""<<hex.str()<<"\",\"hash_encoding\":\"decimal coefficient, six decimal exponents separated byspaces, newline; noheader\",\"full_coefficients_saved\":"<<(written?"true":"false")<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";out.close();front(dir,status,n);cerr<<status<<" TERMS "<<n<<" BYTES "<<bytes<<endl;
 }catch(const exception&e){cerr<<"BOUNDED_RESULTANT_STOP "<<e.what()<<endl;result=2;}
 fmpz_mpoly_clear(a,ctx);fmpz_mpoly_clear(b,ctx);fmpz_mpoly_clear(res,ctx);
 }fmpz_mpoly_ctx_clear(ctx);flint_cleanup();return result;
}
