// Independent generic composition audit. Reversed variable order and DEGREVLEX.
#include <flint/fmpq_mpoly.h>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
void load(istream& in, fmpq_mpoly_t p, const fmpq_mpoly_ctx_t ctx, int n){
 long count;in>>count; assert(in && count>=0);vector<ulong> exp(n);fmpq_t c;fmpq_init(c);
 for(long i=0;i<count;i++){string s;in>>s;assert(fmpq_set_str(c,s.c_str(),10)==0);fmpq_canonicalise(c);
  for(auto& x:exp)in>>x;assert(in);fmpq_mpoly_set_coeff_fmpq_ui(p,c,exp.data(),ctx);}
 assert(fmpq_mpoly_length(p,ctx)==count); fmpq_clear(c);
}
int main(int argc,char**argv){assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);
 fmpq_mpoly_ctx_t source,target;fmpq_mpoly_ctx_init(source,3,ORD_DEGREVLEX);fmpq_mpoly_ctx_init(target,2,ORD_DEGREVLEX);
 fmpq_mpoly_t images[3];fmpq_mpoly_struct* ptrs[3];
 for(int i=0;i<3;i++){fmpq_mpoly_init(images[i],target);load(in,images[i],target,2);ptrs[i]=images[i];}
 long n;in>>n;out<<"{\"all_passed\":true,\"composition_checks\":[";
 for(long i=0;i<n;i++){
  string name;in>>name;fmpq_mpoly_t q,want,got;fmpq_mpoly_init(q,source);fmpq_mpoly_init(want,target);fmpq_mpoly_init(got,target);
  load(in,q,source,3);load(in,want,target,2);
  int ok=fmpq_mpoly_compose_fmpq_mpoly(got,q,ptrs,source,target);
  assert(ok && fmpq_mpoly_equal(got,want,target));
  if(i)out<<',';out<<"{\"name\":\""<<name<<"\",\"source_terms\":"<<fmpq_mpoly_length(q,source)<<",\"image_terms\":"<<fmpq_mpoly_length(got,target)<<"}";
  fmpq_mpoly_clear(q,source);fmpq_mpoly_clear(want,target);fmpq_mpoly_clear(got,target);
 }
 string extra;assert(!(in>>extra));out<<"]}\n";
 for(int i=0;i<3;i++)fmpq_mpoly_clear(images[i],target);fmpq_mpoly_ctx_clear(source);fmpq_mpoly_ctx_clear(target);flint_cleanup();
}
