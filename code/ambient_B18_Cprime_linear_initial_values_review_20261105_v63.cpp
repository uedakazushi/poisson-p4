// Independent PROVIDED rational-matrix products and streamed span hashes.
// No rank search, inverse solve or determinant algorithm.
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <openssl/sha.h>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <sys/resource.h>
#include <unistd.h>
using namespace std;
void scalar(istream&in,fmpq_t x){string s;in>>s;assert(in&&!fmpq_set_str(x,s.c_str(),10));fmpq_canonicalise(x);}
void column(istream&in,fmpq_mat_t x){fmpq_mat_zero(x);long n;in>>n;assert(in&&n>=0&&n<=fmpq_mat_nrows(x));vector<bool>seen(fmpq_mat_nrows(x));for(long i=0;i<n;i++){long r;in>>r;assert(r>=0&&r<fmpq_mat_nrows(x)&&!seen[r]);seen[r]=true;scalar(in,fmpq_mat_entry(x,r,0));assert(!fmpq_is_zero(fmpq_mat_entry(x,r,0)));}}
string hexsum(unsigned char*d){ostringstream s;s<<hex<<setfill('0');for(int i=0;i<SHA256_DIGEST_LENGTH;i++)s<<setw(2)<<unsigned(d[i]);return s.str();}
int main(int argc,char**argv){assert(argc==3);alarm(60);rlimit mem{2000000000,2000000000};setrlimit(RLIMIT_AS,&mem);auto start=chrono::steady_clock::now();ifstream in(argv[1]);ofstream out(argv[2]);int blocks;in>>blocks;assert(blocks==2);out<<"{\"all_passed\":true,\"blocks\":[";long peak=0,total=0;
 for(int b=0;b<blocks;b++){
  long n,m,r;in>>n>>m>>r;assert((b==0&&n==27&&m==44&&r==24)||(b==1&&n==135&&m==396&&r==131));vector<long>rr(r),cc(r);for(auto&x:rr)in>>x;for(auto&x:cc)in>>x;
  fmpq_mat_t U,A,iv,prod,v,sel,fac,back;fmpq_mat_init(U,n,r);fmpq_mat_init(A,r,r);fmpq_mat_init(iv,r,r);fmpq_mat_init(prod,r,r);fmpq_mat_init(v,n,1);fmpq_mat_init(sel,r,1);fmpq_mat_init(fac,r,1);fmpq_mat_init(back,n,1);
  peak=max(peak,n*r+3*r*r+2*n+2*r+12000);
  for(long j=0;j<r;j++){column(in,v);for(long i=0;i<n;i++)fmpq_set(fmpq_mat_entry(U,i,j),fmpq_mat_entry(v,i,0));}
  for(long i=0;i<r;i++)for(long j=0;j<r;j++){assert(rr[i]>=0&&rr[i]<n);fmpq_set(fmpq_mat_entry(A,i,j),fmpq_mat_entry(U,rr[i],j));scalar(in,fmpq_mat_entry(iv,i,j));}
  fmpq_mat_mul(prod,A,iv);for(long i=0;i<r;i++)for(long j=0;j<r;j++)assert(i==j?fmpq_is_one(fmpq_mat_entry(prod,i,j)):fmpq_is_zero(fmpq_mat_entry(prod,i,j)));
  fmpq_mat_mul(prod,iv,A);for(long i=0;i<r;i++)for(long j=0;j<r;j++)assert(i==j?fmpq_is_one(fmpq_mat_entry(prod,i,j)):fmpq_is_zero(fmpq_mat_entry(prod,i,j)));fmpq_mat_clear(A);fmpq_mat_clear(prod);
  SHA256_CTX state;SHA256_Init(&state);
  for(long j=0;j<m;j++){
   column(in,v);for(long k=0;k<r;k++)if(cc[k]==j)for(long i=0;i<n;i++)assert(fmpq_equal(fmpq_mat_entry(v,i,0),fmpq_mat_entry(U,i,k)));
   for(long i=0;i<r;i++)fmpq_set(fmpq_mat_entry(sel,i,0),fmpq_mat_entry(v,rr[i],0));fmpq_mat_mul(fac,iv,sel);fmpq_mat_mul(back,U,fac);assert(fmpq_mat_equal(back,v));total++;
   string line="[";for(long i=0;i<r;i++){if(i)line+=",";char*s=fmpq_get_str(nullptr,10,fmpq_mat_entry(fac,i,0));string value=s;flint_free(s);if(value.size()>2&&value.substr(value.size()-2)=="/1")value.erase(value.size()-2);line+='"';line+=value;line+='"';}line+="]\n";SHA256_Update(&state,line.data(),line.size());assert(chrono::duration<double>(chrono::steady_clock::now()-start).count()<55);
  }
  unsigned char digest[SHA256_DIGEST_LENGTH];SHA256_Final(digest,&state);if(b)out<<',';out<<"{\"rows\":"<<n<<",\"columns\":"<<m<<",\"fixed_rank\":"<<r<<",\"both_inverse_products\":true,\"all_span_columns\":"<<m<<",\"factorization_sha256\":\""<<hexsum(digest)<<"\"}";
  fmpq_mat_clear(U);fmpq_mat_clear(iv);fmpq_mat_clear(v);fmpq_mat_clear(sel);fmpq_mat_clear(fac);fmpq_mat_clear(back);
 }
 string x;assert(!(in>>x));assert(peak<85555);out<<"],\"all_streamed_columns\":"<<total<<",\"owned_matrix_slots_including12000_metadata\":"<<peak<<",\"no_rank_inverse_or_det_algorithm\":true,\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";flint_cleanup();}
