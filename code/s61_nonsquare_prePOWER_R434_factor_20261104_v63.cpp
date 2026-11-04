#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <flint/fmpq_mpoly_factor.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sys/resource.h>
using namespace std;
fmpq_mpoly_ctx_t ctx;chrono::steady_clock::time_point started;long peak=0,factor_terms=0;
struct P{fmpq_mpoly_t p;P(){fmpq_mpoly_init(p,ctx);}P(const P&x):P(){fmpq_mpoly_set(p,x.p,ctx);}P(P&&x):P(){fmpq_mpoly_swap(p,x.p,ctx);}P&operator=(P x){fmpq_mpoly_swap(p,x.p,ctx);return *this;}~P(){fmpq_mpoly_clear(p,ctx);}};
long len(const P&p){return fmpq_mpoly_length(p.p,ctx);}
void guard(long work=0){long live=434+434+21+2*factor_terms+15000+work;peak=max(peak,live);assert(live<=150000);assert(chrono::duration<double>(chrono::steady_clock::now()-started).count()<55);}
void readp(istream&f,P&p){long n;f>>n;assert(f&&n>=0);fmpq_t c;fmpq_init(c);ulong e[3];for(long i=0;i<n;i++){string s;f>>s>>e[0]>>e[1]>>e[2];assert(f&&!fmpq_set_str(c,s.c_str(),10));fmpq_canonicalise(c);fmpq_mpoly_push_term_fmpq_ui(p.p,c,e,ctx);}fmpq_mpoly_sort_terms(p.p,ctx);fmpq_mpoly_combine_like_terms(p.p,ctx);assert(len(p)==n);fmpq_clear(c);}
bool eq(const P&a,const P&b){return fmpq_mpoly_equal(a.p,b.p,ctx);}
P mul(const P&a,const P&b){guard(len(a)+len(b));P z;fmpq_mpoly_mul(z.p,a.p,b.p,ctx);guard(len(a)+len(b)+2*len(z));return z;}
P add(const P&a,const P&b){P z;fmpq_mpoly_add(z.p,a.p,b.p,ctx);guard(len(a)+len(b)+2*len(z));return z;}
void jp(ostream&o,const P&p){o<<'[';fmpq_t c;fmpq_init(c);ulong e[3];for(long i=0;i<len(p);i++){if(i)o<<',';fmpq_mpoly_get_term_exp_ui(e,p.p,i,ctx);fmpq_mpoly_get_term_coeff_fmpq(c,p.p,i,ctx);char*s=fmpq_get_str(nullptr,10,c);o<<"[["<<e[0]<<','<<e[1]<<','<<e[2]<<"],\""<<s<<"\"]";flint_free(s);assert(o.tellp()<30*1024*1024);}o<<']';fmpq_clear(c);}
long degree(const P&p){return fmpq_mpoly_total_degree_si(p.p,ctx);}
int main(int argc,char**argv){assert(argc==3);started=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);{
 ifstream in(argv[1]);P original,primitive,p;readp(in,original);readp(in,primitive);readp(in,p);string sc;in>>sc;assert(in);fmpq_t content;fmpq_init(content);assert(!fmpq_set_str(content,sc.c_str(),10));fmpq_canonicalise(content);P check;fmpq_mpoly_scalar_mul_fmpq(check.p,primitive.p,content,ctx);assert(eq(check,original));guard();
 fmpq_mpoly_factor_t fac;fmpq_mpoly_factor_init(fac,ctx);int status=fmpq_mpoly_factor(fac,primitive.p,ctx);assert(status);assert(fmpq_mpoly_factor_make_monic(fac,ctx));fmpq_mpoly_factor_sort(fac,ctx);guard();
 vector<P>factors;vector<long>exponents;for(long i=0;i<fmpq_mpoly_factor_length(fac,ctx);i++){P f;fmpq_mpoly_factor_get_base(f.p,fac,i,ctx);long e=fmpq_mpoly_factor_get_exp_si(fac,i,ctx);assert(e>0&&degree(f)>0);factor_terms+=len(f);guard();factors.push_back(std::move(f));exponents.push_back(e);}
 fmpq_t constant;fmpq_init(constant);fmpq_mpoly_factor_get_constant_fmpq(constant,fac,ctx);fmpq_mul(constant,constant,content);P prod;fmpq_mpoly_set_fmpq(prod.p,constant,ctx);vector<long>product_terms;for(long i=0;i<long(factors.size());i++)for(long j=0;j<exponents[i];j++){prod=mul(prod,factors[i]);product_terms.push_back(len(prod));}assert(eq(prod,original));
 ofstream out(argv[2]);char*con=fmpq_get_str(nullptr,10,constant);out<<"{\"all_passed\":true,\"factor_API_status\":"<<status<<",\"single_factor_call\":true,\"factor_constant_for_original_R\":\""<<con<<"\",\"factors\":[";flint_free(con);
 for(long i=0;i<long(factors.size());i++){if(i)out<<',';out<<"{\"multiplicity\":"<<exponents[i]<<",\"total_degree\":"<<degree(factors[i])<<",\"terms\":"<<len(factors[i])<<",\"factor\":";jp(out,factors[i]);out<<'}';}out<<"],\"sequential_product_terms\":[";for(long i=0;i<long(product_terms.size());i++){if(i)out<<',';out<<product_terms[i];}out<<"],\"all434_original_coefficients_reproduct_verified\":true,\"named_p_factor_checks\":[";
 P gcd,remaining=p;fmpq_mpoly_one(gcd.p,ctx);for(long i=0;i<long(factors.size());i++){if(i)out<<',';out<<"{\"factor_index\":"<<i<<",\"degree_excludes_dividing_p\":"<<(degree(factors[i])>degree(p)?"true":"false")<<",\"division_steps\":[";long removed=0;
  if(degree(factors[i])<=degree(p))for(long j=0;j<exponents[i];j++){P q,r;fmpq_mpoly_divrem(q.p,r.p,remaining.p,factors[i].p,ctx);assert(eq(add(mul(q,factors[i]),r),remaining));guard(len(q)+len(r));if(j)out<<',';out<<"{\"before\":";jp(out,remaining);out<<",\"quotient\":";jp(out,q);out<<",\"remainder\":";jp(out,r);out<<'}';if(!fmpq_mpoly_is_zero(r.p,ctx))break;gcd=mul(gcd,factors[i]);remaining=q;removed++;}
  out<<"],\"removed_p_multiplicity\":"<<removed<<'}';}
 P rq;assert(fmpq_mpoly_divides(rq.p,original.p,gcd.p,ctx));assert(eq(mul(rq,gcd),original));assert(eq(mul(remaining,gcd),p));out<<"],\"factor_list_gcd_with_p\":";jp(out,gcd);out<<",\"R_quotient_by_factor_list_gcd\":";jp(out,rq);out<<",\"p_quotient_by_factor_list_gcd\":";jp(out,remaining);out<<",\"all_gcd_quotient_products_verified\":true,\"max_accessible_terms_including_source_factors_scratch_reserve\":"<<peak<<",\"factor_internal_support_not_instrumented\":true,\"no_point_prime_other_minor_or_Saux_factor\":true,\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"}\n";out.flush();assert(out.good());fmpq_clear(content);fmpq_clear(constant);fmpq_mpoly_factor_clear(fac,ctx);
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
