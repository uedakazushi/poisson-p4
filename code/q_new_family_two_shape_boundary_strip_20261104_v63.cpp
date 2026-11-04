// ONE exact named-boundary stripping job. NO resultant or general factorization.
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
struct P{fmpz_mpoly_t v;P(){fmpz_mpoly_init(v,ctx);}P(const P&a):P(){fmpz_mpoly_set(v,a.v,ctx);}P(P&&a):P(){fmpz_mpoly_swap(v,a.v,ctx);}P&operator=(P a){fmpz_mpoly_swap(v,a.v,ctx);return*this;}~P(){fmpz_mpoly_clear(v,ctx);}};
string zstr(const fmpz_t z){char*p=fmpz_get_str(nullptr,10,z);string s(p);flint_free(p);return s;}
string line(const P&a,long k){fmpz_t z;fmpz_init(z);ulong e[6];fmpz_mpoly_get_term_coeff_fmpz(z,a.v,k,ctx);fmpz_mpoly_get_term_exp_ui(e,a.v,k,ctx);ostringstream o;o<<zstr(z);for(auto x:e)o<<' '<<x;o<<'\n';fmpz_clear(z);return o.str();}
struct Info{long terms=0,bits=0;uint64_t bytes=0;string digest;slong degree[6];};
Info info(const P&a){Info out;out.terms=fmpz_mpoly_length(a.v,ctx);fmpz_mpoly_degrees_si(out.degree,a.v,ctx);EVP_MD_CTX*d=EVP_MD_CTX_new();assert(d&&EVP_DigestInit_ex(d,EVP_sha256(),nullptr));fmpz_t z;fmpz_init(z);for(long k=0;k<out.terms;k++){string l=line(a,k);out.bytes+=l.size();assert(EVP_DigestUpdate(d,l.data(),l.size()));fmpz_mpoly_get_term_coeff_fmpz(z,a.v,k,ctx);out.bits=max(out.bits,long(fmpz_bits(z)));}unsigned char hash[EVP_MAX_MD_SIZE];unsigned int len;assert(EVP_DigestFinal_ex(d,hash,&len));EVP_MD_CTX_free(d);fmpz_clear(z);ostringstream h;for(unsigned i=0;i<len;i++)h<<hex<<setw(2)<<setfill('0')<<unsigned(hash[i]);out.digest=h.str();return out;}
string jsoninfo(const Info&i){ostringstream o;o<<"{\"terms\":"<<i.terms<<",\"max_coefficient_bits\":"<<i.bits<<",\"canonical_bytes\":"<<i.bytes<<",\"canonical_SHA256\":\""<<i.digest<<"\",\"multidegree\":[";for(int j=0;j<6;j++){if(j)o<<',';o<<i.degree[j];}o<<"]}";return o.str();}
void readpoly(istream&f,P&a){int vars;long n;f>>vars>>n;assert(vars==6&&n>0);fmpz_t z;fmpz_init(z);ulong e[6];for(long k=0;k<n;k++){string c;f>>c;for(auto&v:e)f>>v;assert(f&&fmpz_set_str(z,c.c_str(),10)==0);fmpz_mpoly_push_term_fmpz_ui(a.v,z,e,ctx);}fmpz_mpoly_sort_terms(a.v,ctx);fmpz_mpoly_combine_like_terms(a.v,ctx);fmpz_clear(z);assert(fmpz_mpoly_length(a.v,ctx)==n);}
P mul(const P&a,const P&b){P c;fmpz_mpoly_mul(c.v,a.v,b.v,ctx);return c;}
P power(const P&a,long n){P b;assert(fmpz_mpoly_pow_ui(b.v,a.v,n,ctx));return b;}
P term(long coeff,const vector<ulong>&e){P p;fmpz_mpoly_set_coeff_si_ui(p.v,coeff,e.data(),ctx);return p;}
void divrem(P&q,P&rem,const P&a,const P&b){fmpz_mpoly_divrem(q.v,rem.v,a.v,b.v,ctx);P check=mul(q,b);fmpz_mpoly_add(check.v,check.v,rem.v,ctx);assert(fmpz_mpoly_equal(check.v,a.v,ctx));}
void frontier(const string&dir,const string&stage,const vector<string>&steps){ofstream f(dir+"/native-frontier.json");f<<"{\"stage\":\""<<stage<<"\",\"resultant_call_count\":0,\"verified_steps\":[";for(int i=0;i<int(steps.size());i++){if(i)f<<',';f<<steps[i];}f<<"]}\n";}
long totalbytes=0,budget=0;
void save(const string&dir,const string&name,const P&a,const Info&i){string header="6 "+to_string(i.terms)+"\n";if(totalbytes+header.size()+i.bytes>ulong(budget))throw runtime_error("32MiB total coefficient artifact cap");ofstream f(dir+"/"+name+".txt");f<<header;for(long k=0;k<i.terms;k++)f<<line(a,k);f.close();totalbytes+=header.size()+i.bytes;}
int main(int argc,char**argv){assert(argc==4);fmpz_mpoly_ctx_init(ctx,6,ORD_DEGLEX);int code=0;auto started=chrono::steady_clock::now();{
 try{
  string dir=argv[2];budget=stol(argv[3]);assert(budget>0&&budget<=32L*1024*1024);ifstream in(argv[1]);P raw;readpoly(in,raw);string extra;assert(!(in>>extra));Info rawinfo=info(raw);vector<string>steps;
  ulong minima[6];for(auto&v:minima)v=~ulong(0);ulong e[6];for(long k=0;k<rawinfo.terms;k++){fmpz_mpoly_get_term_exp_ui(e,raw.v,k,ctx);for(int j=0;j<6;j++)minima[j]=min(minima[j],e[j]);}assert(minima[1]==22&&minima[3]==1);P mon=term(1,{0,22,0,1,0,0}),cur,rem;divrem(cur,rem,raw,mon);assert(fmpz_mpoly_is_zero(rem.v,ctx));steps.push_back("{\"factor\":\"r^22*Z0\",\"exact_identity_verified\":true,\"quotient\":"+jsoninfo(info(cur))+"}");frontier(dir,"visible_monomial_removed",steps);
  fmpz_t content,z;fmpz_init(content);fmpz_init(z);fmpz_zero(content);for(long k=0;k<fmpz_mpoly_length(cur.v,ctx);k++){fmpz_mpoly_get_term_coeff_fmpz(z,cur.v,k,ctx);fmpz_gcd(content,content,z);}assert(fmpz_sgn(content)>0);string contentstr=zstr(content);P before(cur);fmpz_mpoly_scalar_divexact_fmpz(cur.v,cur.v,content,ctx);P check;fmpz_mpoly_scalar_mul_fmpz(check.v,cur.v,content,ctx);assert(fmpz_mpoly_equal(check.v,before.v,ctx));steps.push_back("{\"factor\":\"integer_content\",\"value\":\""+contentstr+"\",\"exact_identity_verified\":true,\"quotient\":"+jsoninfo(info(cur))+"}");frontier(dir,"integer_content_removed",steps);
  P rm1=term(1,{0,1,0,0,0,0});P one=term(1,{0,0,0,0,0,0});fmpz_mpoly_sub(rm1.v,rm1.v,one.v,ctx);P T=term(1,{0,1,1,0,0,0});P minus2r=term(-2,{0,1,0,0,0,0}),minusR=term(-1,{0,0,1,0,0,0});fmpz_mpoly_add(T.v,T.v,minus2r.v,ctx);fmpz_mpoly_add(T.v,T.v,minusR.v,ctx);
  long counts[2]={0,0};vector<P>factors={rm1,T};vector<string>names={"r-1","T=(R-2)r-R"};
  for(int f=0;f<2;f++){
   while(true){P q,r;divrem(q,r,cur,factors[f]);if(!fmpz_mpoly_is_zero(r.v,ctx))break;cur=move(q);counts[f]++;Info qi=info(cur);steps.push_back("{\"factor\":\""+names[f]+"\",\"factor_power_so_far\":"+to_string(counts[f])+",\"exact_identity_verified\":true,\"quotient\":"+jsoninfo(qi)+"}");frontier(dir,"named_factor_removed",steps);}
  }
  P factor=mul(mon,mul(power(rm1,counts[0]),power(T,counts[1])));fmpz_mpoly_scalar_mul_fmpz(factor.v,factor.v,content,ctx);P all=mul(factor,cur);assert(fmpz_mpoly_equal(all.v,raw.v,ctx));
  fmpz_zero(z);for(long k=0;k<fmpz_mpoly_length(cur.v,ctx);k++){fmpz_mpoly_get_term_coeff_fmpz(content,cur.v,k,ctx);fmpz_gcd(z,z,content);}assert(fmpz_is_one(z));fmpz_clear(content);fmpz_clear(z);
  Info finalinfo=info(cur);save(dir,"quotient",cur,finalinfo);vector<string>remainders;
  for(int f=0;f<2;f++){P q,r;divrem(q,r,cur,factors[f]);assert(!fmpz_mpoly_is_zero(r.v,ctx));Info qi=info(q),ri=info(r);string label=f==0?"rminus1":"T";save(dir,"nondiv_quotient_"+label,q,qi);save(dir,"nondiv_remainder_"+label,r,ri);remainders.push_back("{\"factor\":\""+names[f]+"\",\"q_times_factor_plus_remainder_equals_final_quotient\":true,\"remainder_nonzero\":true,\"division_quotient\":"+jsoninfo(qi)+",\"division_remainder\":"+jsoninfo(ri)+"}");}
  ofstream out(dir+"/native-report.json");out<<"{\"all_passed\":true,\"resultant_calls\":0,\"general_factorization_calls\":0,\"input\":"<<jsoninfo(rawinfo)<<",\"visible_r_power\":22,\"visible_Z0_power\":1,\"integer_content\":\""<<contentstr<<"\",\"rminus1_power\":"<<counts[0]<<",\"T_power\":"<<counts[1]<<",\"global_factorization_identity_verified\":true,\"final_integer_content_one\":true,\"final_quotient\":"<<jsoninfo(finalinfo)<<",\"all_successful_quotient_step_identities\":[";for(int i=0;i<int(steps.size());i++){if(i)out<<',';out<<steps[i];}out<<"],\"final_maximality_remainder_identities\":[";for(int i=0;i<int(remainders.size());i++){if(i)out<<',';out<<remainders[i];}out<<"],\"all_intermediate_quotients_reconstructible_from_final_and_factor_chain\":true,\"coefficient_artifact_bytes\":"<<totalbytes<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"}\n";out.close();frontier(dir,"all_named_boundary_identities_saved",steps);cerr<<"BOUNDARY_STRIP_PASS rminus1 "<<counts[0]<<" T "<<counts[1]<<" FINAL_TERMS "<<finalinfo.terms<<endl;
 }catch(const exception&e){cerr<<"BOUNDED_NAMED_STRIP_STOP "<<e.what()<<endl;code=2;}
 }fmpz_mpoly_ctx_clear(ctx);flint_cleanup();return code;
}
