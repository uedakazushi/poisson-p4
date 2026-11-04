// Fixed S9/S14/S13 plus old V, no solver or minor test here.
#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
fmpq_mpoly_ctx_t ctx;long peak=0,stored=0;auto started=chrono::steady_clock::now();
struct P{fmpq_mpoly_t v;P(){fmpq_mpoly_init(v,ctx);}P(const P&x):P(){fmpq_mpoly_set(v,x.v,ctx);}P(P&&x):P(){fmpq_mpoly_swap(v,x.v,ctx);}P&operator=(P x){fmpq_mpoly_swap(v,x.v,ctx);return *this;}~P(){fmpq_mpoly_clear(v,ctx);}};
void guard(const P&p){assert(fmpq_mpoly_is_canonical(p.v,ctx));peak=max(peak,fmpq_mpoly_length(p.v,ctx));assert(peak<=150000);assert(chrono::duration<double>(chrono::steady_clock::now()-started).count()<50);}
P one(){P p;fmpq_mpoly_one(p.v,ctx);return p;}
P constant(const string&s){P p;fmpq_t c;fmpq_init(c);assert(!fmpq_set_str(c,s.c_str(),10));fmpq_canonicalise(c);fmpq_mpoly_set_fmpq(p.v,c,ctx);fmpq_clear(c);return p;}
P add(const P&a,const P&b){P p;fmpq_mpoly_add(p.v,a.v,b.v,ctx);guard(p);return p;}
P mul(const P&a,const P&b){P p;fmpq_mpoly_mul(p.v,a.v,b.v,ctx);guard(p);return p;}
P scalar(const P&a,const string&s){P p;fmpq_t c;fmpq_init(c);assert(!fmpq_set_str(c,s.c_str(),10));fmpq_canonicalise(c);fmpq_mpoly_scalar_mul_fmpq(p.v,a.v,c,ctx);fmpq_clear(c);guard(p);return p;}
P shift(const P&a,long dz,long dd){assert(dz>=0&&dd>=0);P p;ulong ex[2];fmpq_t c;fmpq_init(c);for(long i=0;i<fmpq_mpoly_length(a.v,ctx);i++){fmpq_mpoly_get_term_exp_ui(ex,a.v,i,ctx);ex[0]+=dz;ex[1]+=dd;fmpq_mpoly_get_term_coeff_fmpq(c,a.v,i,ctx);fmpq_mpoly_push_term_fmpq_ui(p.v,c,ex,ctx);}fmpq_mpoly_sort_terms(p.v,ctx);fmpq_mpoly_combine_like_terms(p.v,ctx);fmpq_clear(c);guard(p);return p;}
P readtarget(istream&f){long n;f>>n;P p;ulong ex[2];fmpq_t c;fmpq_init(c);for(long i=0;i<n;i++){string s;f>>s>>ex[0]>>ex[1];assert(!fmpq_set_str(c,s.c_str(),10));fmpq_canonicalise(c);fmpq_mpoly_push_term_fmpq_ui(p.v,c,ex,ctx);}fmpq_mpoly_sort_terms(p.v,ctx);fmpq_mpoly_combine_like_terms(p.v,ctx);assert(fmpq_mpoly_length(p.v,ctx)==n);fmpq_clear(c);return p;}
void writepoly(ostream&o,const P&p,long dshift){o<<'[';fmpq_t c;fmpq_init(c);ulong ex[2];for(long i=0;i<fmpq_mpoly_length(p.v,ctx);i++){if(i)o<<',';fmpq_mpoly_get_term_exp_ui(ex,p.v,i,ctx);fmpq_mpoly_get_term_coeff_fmpq(c,p.v,i,ctx);char*s=fmpq_get_str(nullptr,10,c);o<<"[["<<ex[0]<<','<<(long(ex[1])-dshift)<<"],\""<<s<<"\"]";flint_free(s);}fmpq_clear(c);o<<']';}
struct Term{string c;long a,e,d;};
int main(int argc,char**argv){assert(argc==3);started=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,2,ORD_LEX);{
 ifstream f(argv[1]);P Anum=readtarget(f),Enum=readtarget(f);P z=shift(one(),1,0),d=shift(one(),0,1);
 P zp20=add(z,constant("20"));P K=add(mul(zp20,add(z,constant("-10"))),scalar(d,"50"));
 assert(fmpq_mpoly_equal(Anum.v,mul(zp20,K).v,ctx));assert(fmpq_mpoly_equal(Enum.v,add(add(mul(z,z),constant("-400")),scalar(d,"40")).v,ctx));
 vector<P>ap{one()},ep{one()};P A=scalar(Anum,"1/500"),E=scalar(Enum,"1/60");long n;f>>n;assert(n==4);
 vector<string>names;vector<vector<Term>>polys;for(long i=0;i<n;i++){string name;long k;f>>name>>k;names.push_back(name);polys.emplace_back();for(long j=0;j<k;j++){Term t;f>>t.c>>t.a>>t.e>>t.d;assert(t.a>=0&&t.e>=0);polys.back().push_back(t);}}
 string extra;assert(!(f>>extra));ofstream o(argv[2]);o<<"{\"images\":[";bool stopped=false;
 for(long i=0;i<n;i++){
  long clearD=0;for(auto&t:polys[i])clearD=max(clearD,2*t.a+t.e-t.d);P r;
  for(auto&t:polys[i]){while(long(ap.size())<=t.a)ap.push_back(mul(ap.back(),A));while(long(ep.size())<=t.e)ep.push_back(mul(ep.back(),E));P term=scalar(shift(mul(ap[t.a],ep[t.e]),0,t.d+clearD-2*t.a-t.e),t.c);r=add(r,term);}
  stored+=fmpq_mpoly_length(r.v,ctx);assert(stored<=150000);if(i)o<<',';o<<"{\"name\":\""<<names[i]<<"\",\"known_d_denominator_exponent\":"<<clearD<<",\"polynomial\":";writepoly(o,r,clearD);o<<'}';
  // All four fixed source entries are evaluated, even if one is zero.
 }
 o<<"],\"unused_stop_flag\":"<<(stopped?"true":"false")<<",\"stored_image_terms\":"<<stored<<",\"peak_single_polynomial_terms\":"<<peak<<",\"source_substitution_only\":true,\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"}\n";
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
