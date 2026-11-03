// Verify one exact QQ polynomial row/cofactor identity; classify its image
// using ONLY the supplied localization factors. CLI row cofactor units.
#include <flint/fmpq_mpoly.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cassert>
using namespace std;fmpq_mpoly_ctx_t ctx;
struct Poly{fmpq_mpoly_t p;Poly(){fmpq_mpoly_init(p,ctx);}Poly(const Poly&q):Poly(){fmpq_mpoly_set(p,q.p,ctx);}Poly(Poly&&q):Poly(){fmpq_mpoly_swap(p,q.p,ctx);}Poly&operator=(Poly q){fmpq_mpoly_swap(p,q.p,ctx);return *this;}~Poly(){fmpq_mpoly_clear(p,ctx);}};
void rd(istream&f,Poly&p){long nt;f>>nt;assert(f&&nt>=0);fmpq_t c;fmpq_init(c);ulong ex[3];for(long j=0;j<nt;j++){string s;f>>s>>ex[0]>>ex[1]>>ex[2];assert(f&&!fmpq_set_str(c,s.c_str(),10));fmpq_canonicalise(c);fmpq_mpoly_push_term_fmpq_ui(p.p,c,ex,ctx);}fmpq_mpoly_sort_terms(p.p,ctx);fmpq_mpoly_combine_like_terms(p.p,ctx);assert(fmpq_mpoly_length(p.p,ctx)==nt);fmpq_clear(c);}
int main(int argc,char**argv){assert(argc==4);fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);{ifstream f(argv[1]),w(argv[2]),u(argv[3]);long n,m,k,wm;f>>n>>m>>k;w>>wm;assert(n==1&&k==0&&m==wm);vector<Poly>A(m),L(m);for(auto&p:A)rd(f,p);for(auto&p:L)rd(w,p);string extra;assert(!(f>>extra)&&!(w>>extra));Poly sum,tmp;for(long j=0;j<m;j++){fmpq_mpoly_mul(tmp.p,A[j].p,L[j].p,ctx);fmpq_mpoly_add(sum.p,sum.p,tmp.p,ctx);}assert(!fmpq_mpoly_is_zero(sum.p,ctx));const char*vars[]={"a","e","d"};char*target=fmpq_mpoly_get_str_pretty(sum.p,vars,ctx);cout<<"EXACT_SCALAR_IMAGE "<<target<<'\n';flint_free(target);vector<int>values;string text;while(getline(u,text)){if(text.empty())continue;Poly q;assert(!fmpq_mpoly_set_str_pretty(q.p,text.c_str(),vars,ctx));assert(!fmpq_mpoly_is_fmpq(q.p,ctx));int e=0;while(fmpq_mpoly_divides(tmp.p,sum.p,q.p,ctx)){sum=tmp;e++;}values.push_back(e);}assert(fmpq_mpoly_is_fmpq(sum.p,ctx));fmpq_t scalar;fmpq_init(scalar);fmpq_mpoly_get_fmpq(scalar,sum.p,ctx);assert(!fmpq_is_zero(scalar));char*s=fmpq_get_str(nullptr,10,scalar);cout<<"EXACT_KNOWN_UNIT_SCALAR "<<s<<" EXPONENTS";flint_free(s);for(int e:values)cout<<' '<<e;cout<<"\nEXACT_FULL_QQ_SCALAR_UNIT_IDENTITY_PASS\n";fmpq_clear(scalar);}fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
