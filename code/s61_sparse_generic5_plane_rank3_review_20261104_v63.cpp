// Independent certificate interpreter. No pivot selection, GCD search, or solver.
#include <cassert>
#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <gmp.h>
#include <flint/fmpz_mpoly_q.h>
using namespace std;
fmpz_mpoly_ctx_t ctx;
struct Q { fmpz_mpoly_q_t x; Q(){fmpz_mpoly_q_init(x,ctx);} ~Q(){fmpz_mpoly_q_clear(x,ctx);} };
void poly(istream& in,fmpz_mpoly_t p){
  fmpz_mpoly_zero(p,ctx);long n; in>>n; assert(in && n>=0); fmpz_t c;fmpz_init(c);ulong e[5];
  for(long k=0;k<n;k++){string s;in>>s;assert(!fmpz_set_str(c,s.c_str(),10));for(auto&v:e)in>>v;fmpz_mpoly_set_coeff_fmpz_ui(p,c,e,ctx);}
  fmpz_clear(c);
}
int main(){
  fmpz_mpoly_ctx_init(ctx,5,ORD_DEGLEX);auto start=chrono::steady_clock::now();
  long checks=0,nonzero=0,arithmetic=0,loads=0,maxnum=0,maxden=0;vector<unique_ptr<Q>> r;
  string op;
  while(cin>>op){
    if(op=="L") {long id;cin>>id;assert(id==(long)r.size());r.emplace_back(make_unique<Q>());poly(cin,fmpz_mpoly_q_numref(r.back()->x));poly(cin,fmpz_mpoly_q_denref(r.back()->x));assert(!fmpz_mpoly_is_zero(fmpz_mpoly_q_denref(r.back()->x),ctx));fmpz_mpoly_q_canonicalise(r.back()->x,ctx);loads++;}
    else if(op=="E") {long a,b;cin>>a>>b;assert(r.at(a)&&r.at(b));if(!fmpz_mpoly_q_equal(r[a]->x,r[b]->x,ctx)){cerr<<"FAILED_EQUALITY "<<checks<<" REG "<<a<<' '<<b<<endl;const char*names[]={"H","L","s","q","t"};fmpz_mpoly_q_print_pretty(r[a]->x,names,ctx);puts("");fmpz_mpoly_q_print_pretty(r[b]->x,names,ctx);puts("");return 2;}checks++;}
    else if(op=="N") {long a;cin>>a;assert(r.at(a)&&!fmpz_mpoly_q_is_zero(r[a]->x,ctx));nonzero++;}
    else if(op=="X") {long a;cin>>a;r.at(a).reset();}
    else {long id,a,b;cin>>id>>a>>b;assert(id==(long)r.size()&&r.at(a)&&r.at(b));r.emplace_back(make_unique<Q>());
      if(op=="A")fmpz_mpoly_q_add(r[id]->x,r[a]->x,r[b]->x,ctx);
      else if(op=="S")fmpz_mpoly_q_sub(r[id]->x,r[a]->x,r[b]->x,ctx);
      else if(op=="M")fmpz_mpoly_q_mul(r[id]->x,r[a]->x,r[b]->x,ctx);
      else {assert(op=="D"&&!fmpz_mpoly_q_is_zero(r[b]->x,ctx));fmpz_mpoly_q_div(r[id]->x,r[a]->x,r[b]->x,ctx);}arithmetic++;
    }
    if(op!="E"&&op!="N"&&op!="X") {maxnum=max(maxnum,fmpz_mpoly_length(fmpz_mpoly_q_numref(r.back()->x),ctx));maxden=max(maxden,fmpz_mpoly_length(fmpz_mpoly_q_denref(r.back()->x),ctx));assert(maxnum<150000&&maxden<150000);}
  }
  assert(cin.eof());cout<<"{\"all_passed\":true,\"equalities\":"<<checks<<",\"nonzero_checks\":"<<nonzero<<",\"rational_operations\":"<<arithmetic<<",\"loaded_rationals\":"<<loads<<",\"max_num_terms\":"<<maxnum<<",\"max_den_terms\":"<<maxden<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
  r.clear();fmpz_mpoly_ctx_clear(ctx);flint_cleanup();
}
