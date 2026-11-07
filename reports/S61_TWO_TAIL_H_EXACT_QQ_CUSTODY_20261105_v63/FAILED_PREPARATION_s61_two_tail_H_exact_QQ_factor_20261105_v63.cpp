#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_factor.h>
#include <flint/mpoly.h>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
using namespace std;
namespace fs=std::filesystem;
static fmpz_mpoly_ctx_t ctx;
static chrono::steady_clock::time_point begun;
static string directory,phase="not_started";
static long base_RAW=0,peak=0,status=-1;
static double soft=45;
static const long RAW=4194304,OWN=50000;
static void require(bool b,const string&s){if(!b)throw runtime_error(s);}
static double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-begun).count();}
static long files(){long z=0;for(auto&p:fs::directory_iterator(directory))if(p.is_regular_file())z+=p.file_size();return z;}
static void guard(long owned=12000){peak=max(peak,owned);require(owned<=OWN,"accessible ownership guard");require(base_RAW+files()+65536<RAW,"RAW guard");require(elapsed()<soft,"soft clock guard outside opaque library");}
static void journal(const string&s){phase=s;ofstream f(directory+"/native-phase.jsonl",ios::app);f<<"{\"phase\":\""<<s<<"\",\"seconds\":"<<elapsed()<<"}\n";f.flush();}
static long cells(const fmpz_mpoly_struct*p){return p->alloc*(1+mpoly_words_per_exp(p->bits,ctx->minfo));}
static void checkedWrite(ostream&o,const string&s){require(base_RAW+files()+long(s.size())+65536<RAW,"output RAW guard");o<<s;require(bool(o),"output I/O");}
int main(int argc,char**argv){
 if(argc==2&&string(argv[1])=="--help"){cout<<"--run FIXED_H94 OUTPUTDIR BASE_RAW SOFT_SECONDS; ONE3var fmpz_mpoly_factor only\n";return 0;}
 if(argc!=6||string(argv[1])!="--run")return 2;
 directory=argv[3];base_RAW=stol(argv[4]);soft=stod(argv[5]);begun=chrono::steady_clock::now();
 rlimit ram{1900000000,1900000000};setrlimit(RLIMIT_AS,&ram);rlimit cpu{52,52};setrlimit(RLIMIT_CPU,&cpu);
 fmpz_mpoly_ctx_init(ctx,3,ORD_LEX);fmpz_mpoly_t H,A,B;fmpz_mpoly_init(H,ctx);fmpz_mpoly_init(A,ctx);fmpz_mpoly_init(B,ctx);
 fmpz_mpoly_factor_t fac;fmpz_mpoly_factor_init(fac,ctx);fmpz_t c,g;fmpz_init(c);fmpz_init(g);
 try{
  journal("fixedH94 input before ONE opaque library call");ifstream in(argv[2]);long n;in>>n;require(n==94,"literal94");
  for(long i=0;i<n;i++){ulong e[3];string v;in>>e[0]>>e[1]>>e[2]>>v;require(bool(in)&&e[0]<=6&&e[1]<=9&&e[2]<=4,"input box");require(fmpz_set_str(c,v.c_str(),10)==0&&fmpz_bits(c)<=40,"input40bit");fmpz_mpoly_push_term_fmpz_ui(H,c,e,ctx);}
  string extra;require(!(in>>extra),"input trailing bytes");fmpz_mpoly_sort_terms(H,ctx);fmpz_mpoly_combine_like_terms(H,ctx);
  require(fmpz_mpoly_length(H,ctx)==94&&fmpz_mpoly_total_degree_si(H,ctx)==11,"literalH shape");guard(12000+cells(H));
  journal("ONE fmpz_mpoly_factor opaque3var call");status=fmpz_mpoly_factor(fac,H,ctx);
  journal("opaque factor API returned");guard();require(status==1,"factor API STATUS failure");
  long nf=fmpz_mpoly_factor_length(fac,ctx);require(nf>=1&&nf<=11,"factor count");
  long returned_visible=12000+cells(H)+fac->alloc;
  for(long j=0;j<fac->alloc;j++)returned_visible+=cells(fac->poly+j);
  guard(returned_visible);
  ofstream out(directory+"/canonical-factor-candidates.json");checkedWrite(out,"{\"factor_API_status\":1,\"canonical_candidates_only\":true,\"opaque_internal_growth_not_384_certified\":true,\"factors\":[");
  long aggregate=0,total=0;slong degree_sum[3]={0,0,0};long maxbits=0;
  for(long i=0;i<nf;i++){
   long exp=fmpz_mpoly_factor_get_exp_si(fac,i,ctx);require(exp>=1&&exp<=11,"multiplicity");
   guard(returned_visible+2*cells(fac->poly+i));
   journal("canonical primitive factor "+to_string(i));fmpz_mpoly_factor_get_base(A,fac,i,ctx);
   // This library normalization remains inside the disclosed opaque phase.
   fmpz_mpoly_primitive_part(B,A,ctx);require(fmpz_mpoly_is_canonical(B,ctx),"canonical base");
   long nt=fmpz_mpoly_length(B,ctx);require(nt>=1&&nt<=350,"factor350 box");
   fmpz_mpoly_get_term_coeff_fmpz(c,B,0,ctx);if(fmpz_sgn(c)<0)fmpz_mpoly_neg(B,B,ctx);
   slong deg[3];fmpz_mpoly_degrees_si(deg,B,ctx);long td=fmpz_mpoly_total_degree_si(B,ctx);
   require(td>=1&&td<=11&&deg[0]<=6&&deg[1]<=9&&deg[2]<=4,"factor degrees");
   total+=exp*td;for(int j=0;j<3;j++)degree_sum[j]+=exp*deg[j];aggregate+=nt;require(aggregate<=3850,"factor aggregate");
   long visible=12000+cells(H)+cells(A)+cells(B)+fac->alloc;
   for(long j=0;j<fac->alloc;j++)visible+=cells(fac->poly+j);guard(visible);
   fmpz_zero(g);long fb=0;
   for(long j=0;j<nt;j++){fmpz_mpoly_get_term_coeff_fmpz(c,B,j,ctx);fb=max(fb,long(fmpz_bits(c)));require(fb<=357,"Mignotte357 output guard");fmpz_gcd(g,g,c);}
   require(fmpz_is_one(g),"primitive coefficient content");maxbits=max(maxbits,fb);
   if(i)checkedWrite(out,",");checkedWrite(out,"{\"multiplicity\":"+to_string(exp)+",\"total_degree\":"+to_string(td)+",\"degrees\":["+to_string(deg[0])+","+to_string(deg[1])+","+to_string(deg[2])+"],\"max_bits\":"+to_string(fb)+",\"terms\":[");
   for(long j=0;j<nt;j++){
    guard(visible);ulong e[3];require(fmpz_mpoly_term_exp_fits_ui(B,j,ctx),"exp ui");fmpz_mpoly_get_term_exp_ui(e,B,j,ctx);fmpz_mpoly_get_term_coeff_fmpz(c,B,j,ctx);require(fmpz_bits(c)<=357,"emission357");
    char*v=fmpz_get_str(nullptr,10,c);string line=(j?",":"")+string("[[")+to_string(e[0])+","+to_string(e[1])+","+to_string(e[2])+"],\""+string(v)+"\"]";flint_free(v);checkedWrite(out,line);
   }
   checkedWrite(out,"]}");
  }
  require(total==11&&degree_sum[0]==6&&degree_sum[1]==9&&degree_sum[2]==4,"all factor degree sums");
  checkedWrite(out,"],\"factor_count\":"+to_string(nf)+",\"aggregate_terms\":"+to_string(aggregate)+",\"max_canonical_bits\":"+to_string(maxbits)+",\"native_accessible_peak\":"+to_string(peak)+",\"unit_restoration_not_assumed_before_full_product\":true}\n");out.flush();require(bool(out),"candidate close");
  journal("all primitive positive candidates emitted; identity UNVERIFIED");ofstream done(directory+"/native-status.json");done<<"{\"API_completed\":true,\"candidate_list_completed\":true,\"factor_API_status\":1,\"opaque_internal_growth_not_384_certified\":true,\"accessible_peak\":"<<peak<<",\"seconds\":"<<elapsed()<<"}\n";done.flush();
 }catch(const exception&e){ofstream f(directory+"/native-status.json");f<<"{\"API_completed\":"<<(status==1?"true":"false")<<",\"candidate_list_completed\":false,\"factor_API_status\":"<<status<<",\"guard_or_failure\":true,\"opaque_internal_state_not_serialized\":true,\"phase\":\""<<phase<<"\",\"reason\":\""<<e.what()<<"\",\"seconds\":"<<elapsed()<<"}\n";cerr<<e.what()<<"\n";return 3;}
 fmpz_clear(c);fmpz_clear(g);fmpz_mpoly_factor_clear(fac,ctx);fmpz_mpoly_clear(B,ctx);fmpz_mpoly_clear(A,ctx);fmpz_mpoly_clear(H,ctx);fmpz_mpoly_ctx_clear(ctx);flint_cleanup();return 0;
}
