// Remove ONLY declared common unit factors from complete polynomial columns.
// CLI: input output known-units trace. Metadata is preserved. Every accepted
// division is reversed against the complete previous polynomial, including all
// coefficients. No unknown factor or function-field denominator is inverted.
#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;fmpq_mpoly_ctx_t ctx;
struct Poly{fmpq_mpoly_t p;Poly(){fmpq_mpoly_init(p,ctx);}Poly(const Poly&a):Poly(){fmpq_mpoly_set(p,a.p,ctx);}Poly(Poly&&a):Poly(){fmpq_mpoly_swap(p,a.p,ctx);}Poly&operator=(Poly a){fmpq_mpoly_swap(p,a.p,ctx);return *this;}~Poly(){fmpq_mpoly_clear(p,ctx);}};
void rd(istream&in,Poly&p){long n;in>>n;assert(in&&n>=0);fmpq_t q;fmpq_init(q);ulong ex[3];for(long t=0;t<n;t++){string s;in>>s>>ex[0]>>ex[1]>>ex[2];assert(in&&!fmpq_set_str(q,s.c_str(),10));fmpq_canonicalise(q);fmpq_mpoly_push_term_fmpq_ui(p.p,q,ex,ctx);}fmpq_mpoly_sort_terms(p.p,ctx);fmpq_mpoly_combine_like_terms(p.p,ctx);fmpq_clear(q);}
void wr(ostream&out,const Poly&p){fmpq_t q;fmpq_init(q);ulong ex[3];long n=fmpq_mpoly_length(p.p,ctx);out<<n;for(long t=0;t<n;t++){fmpq_mpoly_get_term_coeff_fmpq(q,p.p,t,ctx);fmpq_mpoly_get_term_exp_ui(ex,p.p,t,ctx);char*c=fmpq_get_str(nullptr,10,q);out<<' '<<c<<' '<<ex[0]<<' '<<ex[1]<<' '<<ex[2];flint_free(c);}out<<'\n';fmpq_clear(q);}
int main(int argc,char**argv){assert(argc==5);fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);{
 ifstream in(argv[1]),uf(argv[3]);ofstream out(argv[2]),trace(argv[4]);int n,m,k;in>>n>>m>>k;vector<int>used(k),ids(n);for(int&x:used)in>>x;for(int&x:ids)in>>x;vector<vector<Poly>>A(n,vector<Poly>(m));long before=0;for(auto&row:A)for(auto&p:row){rd(in,p);before+=fmpq_mpoly_length(p.p,ctx);}
 vector<Poly>units;string s;const char*names[]={"a","e","d"};while(getline(uf,s)){if(s.empty())continue;Poly q;assert(!fmpq_mpoly_set_str_pretty(q.p,s.c_str(),names,ctx));assert(!fmpq_mpoly_is_fmpq(q.p,ctx));units.push_back(q);}assert(!units.empty());trace<<m<<' '<<units.size()<<'\n';long divisions=0;
 for(int j=0;j<m;j++){
  bool nonzero=false;for(int i=0;i<n;i++)nonzero|=!fmpq_mpoly_is_zero(A[i][j].p,ctx);vector<long>powers(units.size());
  if(nonzero)for(int u=0;u<(int)units.size();u++)while(true){
   vector<Poly>quot(n);bool all=true;for(int i=0;i<n;i++)if(!fmpq_mpoly_divides(quot[i].p,A[i][j].p,units[u].p,ctx)){all=false;break;}if(!all)break;
   Poly reverse;for(int i=0;i<n;i++){fmpq_mpoly_mul(reverse.p,quot[i].p,units[u].p,ctx);assert(fmpq_mpoly_equal(reverse.p,A[i][j].p,ctx));A[i][j]=quot[i];}powers[u]++;divisions++;
  }
  trace<<j;for(long p:powers)trace<<' '<<p;trace<<'\n';cerr<<"EXACT_DECLARED_UNIT_COLUMN "<<j;for(long p:powers)cerr<<' '<<p;cerr<<endl;
 }
 out<<n<<' '<<m<<' '<<k<<'\n';for(int x:used)out<<x<<' ';out<<'\n';for(int x:ids)out<<x<<' ';out<<'\n';long after=0;for(auto&row:A)for(auto&p:row){wr(out,p);after+=fmpq_mpoly_length(p.p,ctx);}cerr<<"EXACT_COLUMN_UNIT_DIAGONAL_PASS DIVISIONS "<<divisions<<" TERMS_BEFORE "<<before<<" TERMS_AFTER "<<after<<endl;
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
