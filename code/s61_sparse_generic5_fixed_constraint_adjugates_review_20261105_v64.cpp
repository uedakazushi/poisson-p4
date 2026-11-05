// Independent integer-polynomial inverse witness checker, no determinant,
// cofactor, rank, GCD, point evaluation or division algorithm.
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <sys/resource.h>
#include <unistd.h>
using namespace std;
fmpz_mpoly_ctx_t C;
struct Z{fmpz_mpoly_t p;Z(){fmpz_mpoly_init(p,C);}~Z(){fmpz_mpoly_clear(p,C);}};
void readpoly(istream&in,Z&v){
 long n;in>>n;assert(in&&n>=0&&n<150000);fmpz_t c;fmpz_init(c);ulong e[5];
 for(long i=0;i<n;i++){string s;in>>s;assert(!fmpz_set_str(c,s.c_str(),10)&&!fmpz_is_zero(c));for(auto&x:e)in>>x;assert(in);fmpz_mpoly_set_coeff_fmpz_ui(v.p,c,e,C);}
 fmpz_clear(c);assert(fmpz_mpoly_length(v.p,C)==n);
}
int main(int argc,char**argv){
 assert(argc==3);alarm(60);rlimit mem{2000000000,2000000000};setrlimit(RLIMIT_AS,&mem);rlimit cpu{55,55};setrlimit(RLIMIT_CPU,&cpu);
 fmpz_mpoly_ctx_init(C,5,ORD_DEGREVLEX);auto began=chrono::steady_clock::now();long checks=0,muls=0,peak=0,nonzero=0,loaded=0;
 {ifstream in(argv[1]);int blocks;in>>blocks;assert(blocks==2);
 for(int b=0;b<blocks;b++){
  string name;in>>name;assert(name==(b==0?"alpha":"beta"));vector<unique_ptr<Z>>H,A;long held=0;
  for(int k=0;k<25;k++){H.emplace_back(make_unique<Z>());readpoly(in,*H.back());held+=fmpz_mpoly_length(H.back()->p,C);loaded++;}
  for(int k=0;k<25;k++){A.emplace_back(make_unique<Z>());readpoly(in,*A.back());held+=fmpz_mpoly_length(A.back()->p,C);loaded++;}
  Z d,tmp,sum;readpoly(in,d);held+=fmpz_mpoly_length(d.p,C);loaded++;assert(!fmpz_mpoly_is_zero(d.p,C));nonzero++;
  for(int order=0;order<2;order++)for(int i=0;i<5;i++)for(int j=0;j<5;j++){
   fmpz_mpoly_zero(sum.p,C);
   for(int k=0;k<5;k++){
    auto&u=order==0?H[5*i+k]:A[5*i+k];auto&v=order==0?A[5*k+j]:H[5*k+j];fmpz_mpoly_mul(tmp.p,u->p,v->p,C);muls++;
    fmpz_mpoly_add(sum.p,sum.p,tmp.p,C);long owned=25000+held+2*fmpz_mpoly_length(tmp.p,C)+3*fmpz_mpoly_length(sum.p,C);peak=max(peak,owned);assert(owned<150000);assert(chrono::duration<double>(chrono::steady_clock::now()-began).count()<55);
   }
   if(i==j)fmpz_mpoly_sub(sum.p,sum.p,d.p,C);assert(fmpz_mpoly_is_zero(sum.p,C));checks++;
  }
 }
 string extra;assert(!(in>>extra));}
 assert(checks==100&&muls==500&&nonzero==2&&loaded==102);fmpz_mpoly_ctx_clear(C);flint_cleanup();
 ofstream out(argv[2]);out<<"{\"all_passed\":true,\"literal_integer_polynomials\":"<<loaded<<",\"both_nonzero_denominators\":"<<nonzero<<",\"all_left_right_identities\":"<<checks<<",\"fixed_polynomial_products\":"<<muls<<",\"owned_slots_including25000_metadata_upper\":"<<peak<<",\"no_determinant_cofactor_rank_or_GCD_algorithm\":true,\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";
}
