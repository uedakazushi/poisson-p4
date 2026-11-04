// ONE point-free coefficient-content job: exact3-variable gcd/cofactors only.
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
string zstr(const fmpz_t z){char*p=fmpz_get_str(nullptr,10,z);string r(p);flint_free(p);return r;}
string termline(const fmpz_mpoly_t a,long k,const fmpz_mpoly_ctx_t c,int vars){fmpz_t z;fmpz_init(z);ulong e[6];fmpz_mpoly_get_term_coeff_fmpz(z,a,k,c);fmpz_mpoly_get_term_exp_ui(e,a,k,c);ostringstream o;o<<zstr(z);for(int j=0;j<vars;j++)o<<' '<<e[j];o<<'\n';fmpz_clear(z);return o.str();}
string info(const fmpz_mpoly_t a,const fmpz_mpoly_ctx_t c,int vars){long n=fmpz_mpoly_length(a,c);slong deg[6];fmpz_mpoly_degrees_si(deg,a,c);EVP_MD_CTX*d=EVP_MD_CTX_new();assert(d&&EVP_DigestInit_ex(d,EVP_sha256(),nullptr));long bytes=0,bits=0;fmpz_t z;fmpz_init(z);for(long k=0;k<n;k++){string line=termline(a,k,c,vars);bytes+=line.size();assert(EVP_DigestUpdate(d,line.data(),line.size()));fmpz_mpoly_get_term_coeff_fmpz(z,a,k,c);bits=max(bits,long(fmpz_bits(z)));}unsigned char hash[EVP_MAX_MD_SIZE];unsigned int len;assert(EVP_DigestFinal_ex(d,hash,&len));EVP_MD_CTX_free(d);fmpz_clear(z);ostringstream h;for(unsigned int i=0;i<len;i++)h<<hex<<setw(2)<<setfill('0')<<unsigned(hash[i]);ostringstream o;o<<"{\"terms\":"<<n<<",\"max_coefficient_bits\":"<<bits<<",\"canonical_bytes\":"<<bytes<<",\"canonical_SHA256\":\""<<h.str()<<"\",\"multidegree\":[";for(int i=0;i<vars;i++){if(i)o<<',';o<<deg[i];}o<<"]}";return o.str();}
void readterms(istream&f,fmpz_mpoly_t a,const fmpz_mpoly_ctx_t c,int vars,long n){fmpz_t z;fmpz_init(z);ulong e[6];for(long k=0;k<n;k++){string word;f>>word;for(int j=0;j<vars;j++)f>>e[j];assert(f&&fmpz_set_str(z,word.c_str(),10)==0);fmpz_mpoly_push_term_fmpz_ui(a,z,e,c);}fmpz_mpoly_sort_terms(a,c);fmpz_mpoly_combine_like_terms(a,c);fmpz_clear(z);assert(fmpz_mpoly_length(a,c)==n);}
void savepoly(const string&path,const fmpz_mpoly_t a,const fmpz_mpoly_ctx_t c,int vars){ofstream f(path);long bytes=0,n=fmpz_mpoly_length(a,c);f<<vars<<' '<<n<<'\n';for(long k=0;k<n;k++){string s=termline(a,k,c,vars);bytes+=s.size();if(bytes>16L*1024*1024)throw runtime_error("single coefficient file exceeds16MiB");f<<s;}}
P mul(const P&a,const P&b){P c;fmpz_mpoly_mul(c.v,a.v,b.v,ctx);return c;}
void exact(const P&a,const P&g,const P&q){P z=mul(g,q);assert(fmpz_mpoly_equal(a.v,z.v,ctx));}
long gcdcalls=0;
P gcd_verified(const P&a,const P&b,vector<string>&log,const string&stage,int index){P g,qa,qb;int status=fmpz_mpoly_gcd_cofactors(g.v,qa.v,qb.v,a.v,b.v,ctx);gcdcalls++;log.push_back("{\"stage\":\""+stage+"\",\"group_index\":"+to_string(index)+",\"gcd_API_return_status\":"+to_string(status)+"}");if(!status)throw runtime_error("FLINT gcd_cofactors returnedstatus0");assert(!fmpz_mpoly_is_zero(g.v,ctx));exact(a,g,qa);exact(b,g,qb);return g;}
void front(const string&dir,const string&stage,int processed,const P&g,const vector<string>&log){savepoly(dir+"/current-content.txt",g.v,ctx,3);ofstream f(dir+"/native-frontier.json");f<<"{\"stage\":\""<<stage<<"\",\"processed_groups\":"<<processed<<",\"gcd_API_calls\":"<<gcdcalls<<",\"current_gcd\":"<<info(g.v,ctx,3)<<",\"API_status_log\":[";for(int i=0;i<int(log.size());i++){if(i)f<<',';f<<log[i];}f<<"]}\n";}
int main(int argc,char**argv){assert(argc==4);fmpz_mpoly_ctx_init(ctx,3,ORD_DEGLEX);fmpz_mpoly_ctx_t six;fmpz_mpoly_ctx_init(six,6,ORD_DEGLEX);int result=0;auto start=chrono::steady_clock::now();{
 try{
  string dir=argv[3];ifstream in(argv[1]);int vars,ng;in>>vars>>ng;assert(vars==3&&ng==15);vector<P>groups(ng);vector<pair<int,int>>labels(ng);vector<long>counts;
  for(int i=0;i<ng;i++){long n;in>>labels[i].first>>labels[i].second>>n;assert(n>0);readterms(in,groups[i].v,ctx,3,n);counts.push_back(n);}string extra;assert(!(in>>extra));assert(labels[0]==make_pair(0,3)&&labels[1]==make_pair(5,0));assert(fmpz_mpoly_degree_si(groups[0].v,2,ctx)==0&&fmpz_mpoly_degree_si(groups[1].v,2,ctx)==0);
  vector<string>api;P g=gcd_verified(groups[0],groups[1],api,"initial_two_pure_groups",1);front(dir,"initial_pure_coefficient_gcd",2,g,api);
  for(int i=2;i<ng;i++){P q;if(fmpz_mpoly_divides(q.v,groups[i].v,g.v,ctx))exact(groups[i],g,q);else g=gcd_verified(g,groups[i],api,"remaining_group_update",i);front(dir,"content_gcd_all_coefficient_groups",i+1,g,api);}
  assert(fmpz_mpoly_degree_si(g.v,2,ctx)==0);savepoly(dir+"/content.txt",g.v,ctx,3);vector<P>cofactors(ng);for(int i=0;i<ng;i++){assert(fmpz_mpoly_divides(cofactors[i].v,groups[i].v,g.v,ctx));exact(groups[i],g,cofactors[i]);}
  {ofstream f(dir+"/cofactors.txt");f<<"3 "<<ng<<'\n';long bytes=0;for(int i=0;i<ng;i++){long n=fmpz_mpoly_length(cofactors[i].v,ctx);f<<labels[i].first<<' '<<labels[i].second<<' '<<n<<'\n';for(long k=0;k<n;k++){string line=termline(cofactors[i].v,k,ctx,3);bytes+=line.size();if(bytes>16L*1024*1024)throw runtime_error("cofactors artifact16MiBcap");f<<line;}}}
  // Verify literal six-variable Q = content * G, independently of the grouping writer.
  fmpz_mpoly_t raw,regroup,g6,G6,product;for(auto*a:{raw,regroup,g6,G6,product})fmpz_mpoly_init(a,six);
  {ifstream f(argv[2]);int nv;long n;f>>nv>>n;assert(nv==6);readterms(f,raw,six,6,n);assert(!(f>>extra));}
  fmpz_t z;fmpz_init(z);ulong e3[3],e6[6];
  auto embed=[&](const P&p,fmpz_mpoly_t dest,int zi,int zj){for(long k=0;k<fmpz_mpoly_length(p.v,ctx);k++){fmpz_mpoly_get_term_coeff_fmpz(z,p.v,k,ctx);fmpz_mpoly_get_term_exp_ui(e3,p.v,k,ctx);e6[0]=0;e6[1]=e3[0];e6[2]=e3[1];e6[3]=e3[2];e6[4]=zi;e6[5]=zj;fmpz_mpoly_push_term_fmpz_ui(dest,z,e6,six);}};
  for(int i=0;i<ng;i++){embed(groups[i],regroup,labels[i].first,labels[i].second);embed(cofactors[i],G6,labels[i].first,labels[i].second);}embed(g,g6,0,0);
  for(auto*a:{regroup,g6,G6}){fmpz_mpoly_sort_terms(a,six);fmpz_mpoly_combine_like_terms(a,six);}assert(fmpz_mpoly_equal(raw,regroup,six));fmpz_mpoly_mul(product,g6,G6,six);assert(fmpz_mpoly_equal(product,raw,six));savepoly(dir+"/primitive.txt",G6,six,6);fmpz_clear(z);
  // Exact native maximality check on ALL quotient coefficient polynomials.
  P qg=cofactors[0];for(int i=1;i<ng&&!fmpz_mpoly_is_one(qg.v,ctx);i++){P q;if(fmpz_mpoly_divides(q.v,cofactors[i].v,qg.v,ctx))exact(cofactors[i],qg,q);else qg=gcd_verified(qg,cofactors[i],api,"quotient_coefficient_gcd_check",i);}
  bool unit=fmpz_mpoly_is_one(qg.v,ctx);savepoly(dir+"/quotient-coefficient-gcd.txt",qg.v,ctx,3);front(dir,"all_content_identities_and_native_maximality_finished",ng,g,api);
  ofstream out(dir+"/native-report.json");out<<"{\"all_passed\":"<<(unit?"true":"false")<<",\"content\":"<<info(g.v,ctx,3)<<",\"primitive\":"<<info(G6,six,6)<<",\"all15_exact_cofactor_identities_verified\":true,\"source_grouping_all6919_terms_verified\":true,\"full_Q_equals_content_times_primitive_verified\":true,\"content_independent_of_Z0\":true,\"quotient_coefficient_gcd_is_one\":"<<(unit?"true":"false")<<",\"quotient_coefficient_gcd\":"<<info(qg.v,ctx,3)<<",\"gcd_API_call_count\":"<<gcdcalls<<",\"API_status_log\":[";for(int i=0;i<int(api.size());i++){if(i)out<<',';out<<api[i];}out<<"],\"group_labels\":[";for(int i=0;i<ng;i++){if(i)out<<',';out<<'['<<labels[i].first<<','<<labels[i].second<<']';}out<<"],\"cofactor_infos\":[";for(int i=0;i<ng;i++){if(i)out<<',';out<<info(cofactors[i].v,ctx,3);}out<<"],\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";out.close();
  for(auto*a:{raw,regroup,g6,G6,product})fmpz_mpoly_clear(a,six);cerr<<"CONTENT_GCD_PASS unit "<<unit<<" calls "<<gcdcalls<<endl;if(!unit)result=2;
 }catch(const exception&e){cerr<<"BOUNDED_CONTENT_STOP "<<e.what()<<endl;result=2;}
 }fmpz_mpoly_ctx_clear(six);fmpz_mpoly_ctx_clear(ctx);flint_cleanup();return result;
}
