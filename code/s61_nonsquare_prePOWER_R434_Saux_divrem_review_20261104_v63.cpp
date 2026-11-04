// Fixed multiply-back/normal-form certificate only: NO division or factor API.
#include <flint/fmpq.h>
#include <flint/nmod_mpoly.h>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
using namespace std;constexpr ulong prime=32003;nmod_mpoly_ctx_t ctx;long denchecks=0;
ulong power(ulong a,ulong n){ulong v=1;for(;n;n>>=1,a=a*a%prime)if(n&1)v=v*a%prime;return v;}
ulong qq(const string&s){fmpq_t q;fmpq_init(q);assert(!fmpq_set_str(q,s.c_str(),10));fmpq_canonicalise(q);ulong a=fmpz_fdiv_ui(fmpq_numref(q),prime),b=fmpz_fdiv_ui(fmpq_denref(q),prime);assert(b);denchecks++;fmpq_clear(q);return a*power(b,prime-2)%prime;}
void get(nmod_mpoly_t p,istream&in,bool rational,long expected){long n;in>>n;assert(n==expected);set<string>seen;for(long i=0;i<n;i++){ulong e[3];string c;in>>e[0]>>e[1]>>e[2]>>c;ulong v=rational?qq(c):stoul(c);assert(v<prime);if(!rational)assert(v);string k=to_string(e[0])+","+to_string(e[1])+","+to_string(e[2]);assert(seen.insert(k).second);nmod_mpoly_set_coeff_ui_ui(p,v,e,ctx);}}
int main(int argc,char**argv){assert(argc==3);auto start=chrono::steady_clock::now();nmod_mpoly_ctx_init(ctx,3,ORD_DEGREVLEX,prime);{
 nmod_mpoly_t R,S,q,rem,product,total;nmod_mpoly_init(R,ctx);nmod_mpoly_init(S,ctx);nmod_mpoly_init(q,ctx);nmod_mpoly_init(rem,ctx);nmod_mpoly_init(product,ctx);nmod_mpoly_init(total,ctx);
 ifstream in(argv[1]);get(R,in,true,434);get(S,in,true,7022);get(q,in,false,3212);get(rem,in,false,6729);
 ulong expectedLM[3],lc;in>>expectedLM[0]>>expectedLM[1]>>expectedLM[2]>>lc;assert(lc==27306);string extra;assert(!(in>>extra));
 ulong LM[3];nmod_mpoly_get_term_exp_ui(LM,R,0,ctx);assert(LM[0]==2&&LM[1]==15&&LM[2]==3);for(int i=0;i<3;i++)assert(LM[i]==expectedLM[i]);assert(nmod_mpoly_get_term_coeff_ui(R,0,ctx)==lc);
 assert(nmod_mpoly_length(q,ctx)==3212&&nmod_mpoly_length(rem,ctx)==6729&&!nmod_mpoly_is_zero(rem,ctx));
 ulong m[3];for(long i=0;i<nmod_mpoly_length(rem,ctx);i++){nmod_mpoly_get_term_exp_ui(m,rem,i,ctx);assert(!(m[0]>=LM[0]&&m[1]>=LM[1]&&m[2]>=LM[2]));}
 nmod_mpoly_mul(product,q,R,ctx);nmod_mpoly_add(total,product,rem,ctx);assert(nmod_mpoly_equal(total,S,ctx));
 long visible=15000;for(auto*p:{R,S,q,rem,product,total})visible+=nmod_mpoly_length(p,ctx);assert(visible<150000);
 ofstream out(argv[2]);out<<"{\"all_passed\":true,\"prime\":32003,\"order\":\"DEGREVLEX(a,e,d)\",\"division_calls\":0,\"factor_calls\":0,\"source_R_LM\":[2,15,3],\"source_R_LC\":27306,\"quotient_terms\":3212,\"remainder_terms\":6729,\"fixed_qR_plus_remainder_equals_full_Saux\":true,\"all6729_remainder_terms_normal\":true,\"remainder_nonzero\":true,\"QQ_denominator_checks\":"<<denchecks<<",\"visible_coefficients_with_reserve\":"<<visible<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
 nmod_mpoly_clear(R,ctx);nmod_mpoly_clear(S,ctx);nmod_mpoly_clear(q,ctx);nmod_mpoly_clear(rem,ctx);nmod_mpoly_clear(product,ctx);nmod_mpoly_clear(total,ctx);
 }nmod_mpoly_ctx_clear(ctx);flint_cleanup();}
