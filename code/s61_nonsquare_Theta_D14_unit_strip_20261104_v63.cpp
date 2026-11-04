// Fixed gcd/division ONLY against six already certified unit localizers.
#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;fmpq_mpoly_ctx_t ctx;long peak=0,stored=0;auto start=chrono::steady_clock::now();
struct P{fmpq_mpoly_t v;P(){fmpq_mpoly_init(v,ctx);}P(const P&a):P(){fmpq_mpoly_set(v,a.v,ctx);}P(P&&a):P(){fmpq_mpoly_swap(v,a.v,ctx);}P&operator=(P a){fmpq_mpoly_swap(v,a.v,ctx);return*this;}~P(){fmpq_mpoly_clear(v,ctx);}};
void guard(const P&p){assert(fmpq_mpoly_is_canonical(p.v,ctx));peak=max(peak,fmpq_mpoly_length(p.v,ctx));assert(peak<=150000);assert(chrono::duration<double>(chrono::steady_clock::now()-start).count()<50);}
P mul(const P&a,const P&b){P p;fmpq_mpoly_mul(p.v,a.v,b.v,ctx);guard(p);return p;}
P quotient(const P&a,const P&b){P p;assert(fmpq_mpoly_divides(p.v,a.v,b.v,ctx));assert(fmpq_mpoly_equal(mul(p,b).v,a.v,ctx));guard(p);return p;}
void monic(P&p){assert(!fmpq_mpoly_is_zero(p.v,ctx));fmpq_t c;fmpq_init(c);fmpq_mpoly_get_term_coeff_fmpq(c,p.v,0,ctx);fmpq_mpoly_scalar_div_fmpq(p.v,p.v,c,ctx);fmpq_clear(c);guard(p);}
struct L{P p;long off;};
L read(istream&f){long n;f>>n;struct T{string c;long z,d;};vector<T>v;long off=0;for(long i=0;i<n;i++){T t;f>>t.c>>t.z>>t.d;assert(t.z>=0);if(i==0)off=t.d;else off=min(off,t.d);v.push_back(t);}P p;fmpq_t c;fmpq_init(c);ulong ex[2];for(auto&t:v){assert(!fmpq_set_str(c,t.c.c_str(),10));fmpq_canonicalise(c);ex[0]=t.z;ex[1]=t.d-off;fmpq_mpoly_push_term_fmpq_ui(p.v,c,ex,ctx);}fmpq_mpoly_sort_terms(p.v,ctx);fmpq_mpoly_combine_like_terms(p.v,ctx);fmpq_clear(c);assert(fmpq_mpoly_length(p.v,ctx)==n);return{p,off};}
void jp(ostream&o,const P&p,long off=0){stored+=fmpq_mpoly_length(p.v,ctx);assert(stored<=150000);o<<'[';fmpq_t c;fmpq_init(c);ulong ex[2];for(long i=0;i<fmpq_mpoly_length(p.v,ctx);i++){if(i)o<<',';fmpq_mpoly_get_term_coeff_fmpq(c,p.v,i,ctx);fmpq_mpoly_get_term_exp_ui(ex,p.v,i,ctx);char*s=fmpq_get_str(nullptr,10,c);o<<"[["<<ex[0]<<','<<long(ex[1])+off<<"],\""<<s<<"\"]";flint_free(s);}o<<']';fmpq_clear(c);}
int main(int argc,char**argv){assert(argc==3);start=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,2,ORD_LEX);{
 ifstream f(argv[1]);long n;f>>n;assert(n==2);vector<string>names;vector<L>inputs;for(int i=0;i<n;i++){string s;f>>s;names.push_back(s);inputs.push_back(read(f));}long k;f>>k;assert(k==6);vector<string>fn;vector<P>units;for(int i=0;i<k;i++){string s;f>>s;L p=read(f);assert(p.off==0);fn.push_back(s);units.push_back(p.p);}string extra;assert(!(f>>extra));ofstream o(argv[2]);o<<"{\"results\":[";
 for(int i=0;i<n;i++){if(i)o<<',';P r=inputs[i].p,original=r;assert(!fmpq_mpoly_is_zero(r.v,ctx));fmpq_t content;fmpq_init(content);fmpq_mpoly_content(content,r.v,ctx);fmpq_mpoly_scalar_div_fmpq(r.v,r.v,content,ctx);P primitive=r;vector<P>gs;char*s=fmpq_get_str(nullptr,10,content);o<<"{\"name\":\""<<names[i]<<"\",\"d_unit_exponent\":"<<inputs[i].off<<",\"QQ_content\":\""<<s<<"\",\"source\":";flint_free(s);jp(o,original,inputs[i].off);o<<",\"steps\":[";bool first=true;
  for(int j=0;j<k;j++){for(int rep=0;;rep++){assert(rep<100);P g;assert(fmpq_mpoly_gcd(g.v,r.v,units[j].v,ctx));monic(g);if(fmpq_mpoly_is_one(g.v,ctx))break;P q=quotient(r,g),u=quotient(units[j],g);gs.push_back(g);if(!first)o<<',';first=false;o<<"{\"localizer\":\""<<fn[j]<<"\",\"gcd\":";jp(o,g);o<<",\"source_quotient\":";jp(o,q);o<<",\"localizer_quotient\":";jp(o,u);o<<'}';r=q;}}
  P reconstructed=r;for(auto&g:gs)reconstructed=mul(reconstructed,g);assert(fmpq_mpoly_equal(reconstructed.v,primitive.v,ctx));fmpq_mpoly_scalar_mul_fmpq(reconstructed.v,reconstructed.v,content,ctx);assert(fmpq_mpoly_equal(reconstructed.v,original.v,ctx));o<<"],\"remaining_unknown\":";jp(o,r);o<<",\"remaining_is_nonzero_QQ_constant\":"<<(fmpq_mpoly_is_fmpq(r.v,ctx)?"true":"false")<<",\"all_source_and_localizer_quotient_identities_verified\":true}";fmpq_clear(content);
 }o<<"],\"stored_output_terms\":"<<stored<<",\"peak_polynomial_terms\":"<<peak<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
