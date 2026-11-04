#include <flint/nmod_mpoly.h>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;nmod_mpoly_ctx_t ctx;chrono::steady_clock::time_point started;long peak=0;
struct P{nmod_mpoly_t p;P(){nmod_mpoly_init(p,ctx);}~P(){nmod_mpoly_clear(p,ctx);}void release(){nmod_mpoly_clear(p,ctx);nmod_mpoly_init(p,ctx);}};
long len(const P&p){return nmod_mpoly_length(p.p,ctx);}
void guard(initializer_list<const P*>ps,long scratch=0){long s=10000+scratch;for(auto p:ps)s+=len(*p);peak=max(peak,s);assert(s<=150000);assert(chrono::duration<double>(chrono::steady_clock::now()-started).count()<55);}
void readp(istream&f,P&p){long n;f>>n;assert(f&&n>=0);ulong e[3],c;for(long i=0;i<n;i++){f>>c>>e[0]>>e[1]>>e[2];assert(f&&c>0&&c<32003);nmod_mpoly_push_term_ui_ui(p.p,c,e,ctx);}nmod_mpoly_sort_terms(p.p,ctx);nmod_mpoly_combine_like_terms(p.p,ctx);assert(len(p)==n);}
void jp(ostream&o,const P&p){o<<'[';ulong e[3];for(long i=0;i<len(p);i++){if(i)o<<',';nmod_mpoly_get_term_exp_ui(e,p.p,i,ctx);o<<"[["<<e[0]<<','<<e[1]<<','<<e[2]<<"],"<<nmod_mpoly_get_term_coeff_ui(p.p,i,ctx)<<']';assert(o.tellp()<29*1024*1024);}o<<']';}
int main(int argc,char**argv){assert(argc==3);started=chrono::steady_clock::now();nmod_mpoly_ctx_init(ctx,3,ORD_DEGREVLEX,32003);{
 ifstream in(argv[1]);P f,R,q,r,t,product;readp(in,f);readp(in,R);assert(nmod_mpoly_total_degree_si(f.p,ctx)==61&&nmod_mpoly_total_degree_si(R.p,ctx)==20);ulong lm[3];nmod_mpoly_get_term_exp_ui(lm,R.p,0,ctx);assert(lm[0]+lm[1]+lm[2]==20);guard({&f,&R});
 nmod_mpoly_divrem(q.p,r.p,f.p,R.p,ctx);guard({&f,&R,&q,&r});assert(len(q)<=13244&&len(r)<=28420);assert(nmod_mpoly_total_degree_si(q.p,ctx)<=41&&nmod_mpoly_total_degree_si(r.p,ctx)<=61);
 ulong e[3];for(long i=0;i<len(r);i++){nmod_mpoly_get_term_exp_ui(e,r.p,i,ctx);assert(e[0]<lm[0]||e[1]<lm[1]||e[2]<lm[2]);}
 ofstream out(argv[2]);out<<"{\"prime\":32003,\"order\":\"DEGREVLEX(a,e,d)\",\"single_divrem_call\":true,\"R_leading_monomial\":["<<lm[0]<<','<<lm[1]<<','<<lm[2]<<"],\"R_leading_coefficient\":"<<nmod_mpoly_get_term_coeff_ui(R.p,0,ctx)<<",\"quotient\":";jp(out,q);out<<",\"remainder\":";jp(out,r);bool remzero=nmod_mpoly_is_zero(r.p,ctx);long qn=len(q),rn=len(r);out.flush();
 nmod_mpoly_sub(t.p,f.p,r.p,ctx);guard({&f,&R,&q,&r,&t},len(t));assert(len(t)<=35442);r.release();guard({&f,&R,&q,&t});
 nmod_mpoly_mul(product.p,R.p,q.p,ctx);guard({&f,&R,&q,&t,&product},len(product));assert(len(product)<=41664);assert(nmod_mpoly_equal(t.p,product.p,ctx));
 out<<",\"all_passed\":true,\"remainder_is_zero\":"<<(remzero?"true":"false")<<",\"Rbar_does_not_divide_Sauxbar\":"<<(remzero?"false":"true")<<",\"whole_fixed_polynomial_identity_verified\":true,\"every_remainder_monomial_normal\":true,\"quotient_terms\":"<<qn<<",\"remainder_terms\":"<<rn<<",\"max_accessible_coefficient_slots_with_scratch_reserve\":"<<peak<<",\"division_internal_support_not_instrumented\":true,\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"}\n";out.flush();assert(out.good());
 }nmod_mpoly_ctx_clear(ctx);flint_cleanup();}
