// One complete symbolic F32003 multivariate factorization, no evaluation.
#include <flint/fmpq.h>
#include <flint/nmod_mpoly.h>
#include <flint/nmod_mpoly_factor.h>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;constexpr ulong prime=32003;nmod_mpoly_ctx_t ctx;
struct P;set<const P*>live;long peak=0;auto start=chrono::steady_clock::now();
struct P{nmod_mpoly_t p;P(){nmod_mpoly_init(p,ctx);live.insert(this);}P(const P&a):P(){nmod_mpoly_set(p,a.p,ctx);}P(P&&a):P(){nmod_mpoly_swap(p,a.p,ctx);}P&operator=(P a){nmod_mpoly_swap(p,a.p,ctx);return *this;}~P(){live.erase(this);nmod_mpoly_clear(p,ctx);}};
void guard(long factor_terms=0){long n=15000+factor_terms;for(auto*p:live)n+=nmod_mpoly_length(p->p,ctx);peak=max(peak,n);if(n>150000)throw runtime_error("visible coefficient guard");if(chrono::duration<double>(chrono::steady_clock::now()-start).count()>48)throw runtime_error("48s guard");}
ulong power(ulong a,ulong n){ulong v=1;for(;n;n>>=1,a=a*a%prime)if(n&1)v=v*a%prime;return v;}
ulong qq(const string&s){fmpq_t q;fmpq_init(q);assert(!fmpq_set_str(q,s.c_str(),10));fmpq_canonicalise(q);ulong n=fmpz_fdiv_ui(fmpq_numref(q),prime),d=fmpz_fdiv_ui(fmpq_denref(q),prime);assert(d);fmpq_clear(q);return n*power(d,prime-2)%prime;}
P readpoly(istream&in,bool integer){long n;in>>n;assert(n==434);P p;for(long i=0;i<n;i++){ulong ex[3];string c;in>>ex[0]>>ex[1]>>ex[2]>>c;ulong x;if(integer){fmpz_t z;fmpz_init(z);assert(!fmpz_set_str(z,c.c_str(),10));x=fmpz_fdiv_ui(z,prime);fmpz_clear(z);}else x=qq(c);nmod_mpoly_set_coeff_ui_ui(p.p,x,ex,ctx);}guard();return p;}
long degree(const P&p){long d=-1;ulong e[3];for(long i=0;i<nmod_mpoly_length(p.p,ctx);i++){nmod_mpoly_get_term_exp_ui(e,p.p,i,ctx);d=max(d,long(e[0]+e[1]+e[2]));}return d;}
void emit(ostream&out,const P&p){out<<'[';ulong e[3];for(long i=0;i<nmod_mpoly_length(p.p,ctx);i++){if(i)out<<',';nmod_mpoly_get_term_exp_ui(e,p.p,i,ctx);out<<"[["<<e[0]<<','<<e[1]<<','<<e[2]<<"],"<<nmod_mpoly_get_term_coeff_ui(p.p,i,ctx)<<']';}out<<']';}
int main(int argc,char**argv){assert(argc==4);nmod_mpoly_ctx_init(ctx,3,ORD_LEX,prime);ofstream stage(argv[3]);try{{
 ifstream in(argv[1]);string cs;in>>cs;ulong content=qq(cs);assert(content);P original=readpoly(in,false),primitive=readpoly(in,true),scaled;nmod_mpoly_scalar_mul_ui(scaled.p,primitive.p,content,ctx);assert(nmod_mpoly_equal(original.p,scaled.p,ctx));assert(degree(original)==20&&degree(primitive)==20);
 string extra;assert(!(in>>extra));stage<<"all434_QQ_and_primitive_coefficients_bound;content_unit;degree20_preserved\n"<<flush;
 nmod_mpoly_factor_t fac;nmod_mpoly_factor_init(fac,ctx);guard();int ok=nmod_mpoly_factor(fac,primitive.p,ctx);stage<<"factor_API_status="<<ok<<'\n'<<flush;assert(ok);nmod_mpoly_factor_sort(fac,ctx);
 long count=nmod_mpoly_factor_length(fac,ctx);assert(count>=1&&count<=20);vector<P>bases;vector<long>exponents;long factor_terms=0,total_degree=0;
 for(long i=0;i<count;i++){P b;nmod_mpoly_factor_get_base(b.p,fac,i,ctx);long e=nmod_mpoly_factor_get_exp_si(fac,i,ctx);assert(e>0&&degree(b)>0);factor_terms+=nmod_mpoly_length(b.p,ctx);total_degree+=degree(b)*e;bases.push_back(move(b));exponents.push_back(e);guard(factor_terms);}
 assert(total_degree==20);ulong fc=nmod_mpoly_factor_get_constant_ui(fac,ctx);assert(fc);P product;nmod_mpoly_set_ui(product.p,fc,ctx);
 for(long i=0;i<count;i++)for(long j=0;j<exponents[i];j++){P temp;nmod_mpoly_mul(temp.p,product.p,bases[i].p,ctx);product=move(temp);assert(degree(product)<=20);guard(factor_terms);}
 assert(nmod_mpoly_equal(product.p,primitive.p,ctx));bool irreducible=count==1&&exponents[0]==1&&degree(bases[0])==20;
 ofstream out(argv[2]);out<<"{\"all_factorization_checks_passed\":true,\"prime\":32003,\"factor_API_status\":"<<ok<<",\"complete_factor_calls\":1,\"degree_preserved\":20,\"content_mod32003\":"<<content<<",\"factor_constant_mod32003\":"<<fc<<",\"irreducible_mod32003\":"<<(irreducible?"true":"false")<<",\"factors\":[";
 for(long i=0;i<count;i++){if(i)out<<',';out<<"{\"multiplicity\":"<<exponents[i]<<",\"degree\":"<<degree(bases[i])<<",\"terms\":"<<nmod_mpoly_length(bases[i].p,ctx)<<",\"polynomial\":";emit(out,bases[i]);out<<'}';}
 out<<"],\"full_factor_product_matches_primitive\":true,\"max_visible_coefficients_with_reserve\":"<<peak<<",\"internal_factor_API_support_not_instrumented\":true,\"no_parameter_evaluation_divisibility_Saux_or_new_prime\":true,\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
 nmod_mpoly_factor_clear(fac,ctx);stage<<"FINAL_product_equal;irreducible="<<irreducible<<'\n'<<flush;
 }}catch(const exception&e){stage<<"STOP "<<e.what()<<'\n'<<flush;cerr<<e.what()<<'\n';nmod_mpoly_ctx_clear(ctx);return 2;}nmod_mpoly_ctx_clear(ctx);flint_cleanup();}
