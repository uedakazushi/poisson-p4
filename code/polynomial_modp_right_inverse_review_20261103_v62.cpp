// Independent outer-product proof of all100 identities, no Singular/GB calls.
#include <flint/nmod_mpoly.h>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
nmod_mpoly_ctx_t context;
constexpr ulong review_prime=32003;
struct Poly{
 nmod_mpoly_t a;
 Poly(){nmod_mpoly_init(a,context);}
 Poly(const Poly&b):Poly(){nmod_mpoly_set(a,b.a,context);}
 Poly(Poly&&b):Poly(){nmod_mpoly_swap(a,b.a,context);}
 Poly&operator=(Poly b){nmod_mpoly_swap(a,b.a,context);return *this;}
 ~Poly(){nmod_mpoly_clear(a,context);}
};
vector<vector<Poly>> input(const char*path,int n,int m){
 ifstream file(path);long nr,nc,unused;file>>nr>>nc>>unused;assert(nr==n&&nc==m&&unused==0);
 vector<vector<Poly>> result(n,vector<Poly>(m));
 for(auto&row:result)for(auto&poly:row){
  long count;file>>count;assert(file&&count>=0);
  for(long k=0;k<count;k++){
   long long c,ex[3];file>>c>>ex[0]>>ex[1]>>ex[2];assert(file);
   assert(ex[0]>=0&&ex[1]>=0&&ex[2]>=0);ulong power[3]={(ulong)ex[0],(ulong)ex[1],(ulong)ex[2]};
   c%=(long long)review_prime;if(c<0)c+=review_prime;assert(c);
   nmod_mpoly_push_term_ui_ui(poly.a,c,power,context);
  }
  nmod_mpoly_sort_terms(poly.a,context);nmod_mpoly_combine_like_terms(poly.a,context);
  assert(nmod_mpoly_length(poly.a,context)==count);
 }
 string tail;assert(!(file>>tail));return result;
}
void term(Poly&poly,long coefficient,ulong a,ulong e,ulong p){
 coefficient%=(long)review_prime;if(coefficient<0)coefficient+=review_prime;ulong ex[3]={a,e,p};
 nmod_mpoly_set_coeff_ui_ui(poly.a,coefficient,ex,context);
}
int main(int argc,char**argv){
 assert(argc==3);nmod_mpoly_ctx_init(context,3,ORD_DEGREVLEX,review_prime);
 {
  auto M=input(argv[1],10,53),L=input(argv[2],53,10);
  vector<vector<Poly>> answer(10,vector<Poly>(10));Poly product;
  // Column outer products, unlike the producer checker's separate entry sums.
  for(int k=0;k<53;k++)for(int i=0;i<10;i++)if(!nmod_mpoly_is_zero(M[i][k].a,context))
   for(int j=0;j<10;j++)if(!nmod_mpoly_is_zero(L[k][j].a,context)){
    nmod_mpoly_mul(product.a,M[i][k].a,L[k][j].a,context);
    nmod_mpoly_add(answer[i][j].a,answer[i][j].a,product.a,context);
   }
  Poly C,J,J2,aP3,temp,g;
  term(C,100,1,0,0);term(C,18,0,2,0);term(C,-9,0,1,0);term(C,-2,0,0,0);
  term(J,5,0,0,1);term(J,-6,0,1,0);term(J,-4,0,0,0);term(aP3,1,1,0,3);
  nmod_mpoly_mul(J2.a,J.a,J.a,context);nmod_mpoly_mul(temp.a,aP3.a,C.a,context);
  nmod_mpoly_mul(g.a,temp.a,J2.a,context);assert(!nmod_mpoly_is_zero(g.a,context));
  int checks=0;
  for(int i=0;i<10;i++)for(int j=0;j<10;j++){
   assert(i==j?nmod_mpoly_equal(answer[i][j].a,g.a,context):nmod_mpoly_is_zero(answer[i][j].a,context));checks++;
  }
  assert(checks==100);cout<<"INDEPENDENT_OUTER_PRODUCT_ALL100_IDENTITIES_PASS "<<checks<<'\n';
 }
 nmod_mpoly_ctx_clear(context);flint_cleanup();
}
