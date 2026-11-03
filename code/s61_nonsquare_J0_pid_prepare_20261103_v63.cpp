// Exact J=0 substitution on the J-safe36x70 cubic cokernel matrix.
// d=(6e+4)/(5a); each source row is multiplied only by a^M to clear a-poles.
#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <cassert>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;fmpq_mpoly_ctx_t ctx;
struct Poly{fmpq_mpoly_t p;Poly(){fmpq_mpoly_init(p,ctx);}Poly(const Poly&a):Poly(){fmpq_mpoly_set(p,a.p,ctx);}Poly(Poly&&a):Poly(){fmpq_mpoly_swap(p,a.p,ctx);}Poly&operator=(Poly a){fmpq_mpoly_swap(p,a.p,ctx);return *this;}~Poly(){fmpq_mpoly_clear(p,ctx);}};
void read(istream&in,Poly&p){int n;in>>n;fmpq_t q;fmpq_init(q);ulong ex[3];for(int i=0;i<n;i++){string c;in>>c>>ex[0]>>ex[1]>>ex[2];assert(!fmpq_set_str(q,c.c_str(),10));fmpq_canonicalise(q);fmpq_mpoly_push_term_fmpq_ui(p.p,q,ex,ctx);}fmpq_mpoly_sort_terms(p.p,ctx);fmpq_mpoly_combine_like_terms(p.p,ctx);fmpq_clear(q);}
void write(ostream&out,const Poly&p){fmpq_t q;fmpq_init(q);ulong ex[3];int n=fmpq_mpoly_length(p.p,ctx);out<<n;for(int i=0;i<n;i++){fmpq_mpoly_get_term_coeff_fmpq(q,p.p,i,ctx);fmpq_mpoly_get_term_exp_ui(ex,p.p,i,ctx);char*c=fmpq_get_str(nullptr,10,q);out<<' '<<c<<' '<<ex[0]<<' '<<ex[1]<<' '<<ex[2];flint_free(c);}out<<'\n';fmpq_clear(q);}
int main(int argc,char**argv){assert(argc==4);ifstream in(argv[1]);ofstream out(argv[2]),meta(argv[3]);fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);{
 int n,m,k;in>>n>>m>>k;assert(n==36&&m==70&&k==21);int label;for(int i=0;i<k+n;i++)in>>label;
 out<<n<<' '<<m<<" 0\n";const char*names[]={"a","e","d"};Poly D;assert(!fmpq_mpoly_set_str_pretty(D.p,"(6*e+4)/5",names,ctx));vector<Poly>powers(1);fmpq_mpoly_one(powers[0].p,ctx);
 fmpq_t c;fmpq_init(c);ulong ex[3];long allterms=0;
 for(int row=0;row<n;row++){
  vector<Poly>source(m);long M=0,maxd=0;for(auto&p:source){read(in,p);for(int j=0;j<fmpq_mpoly_length(p.p,ctx);j++){fmpq_mpoly_get_term_exp_ui(ex,p.p,j,ctx);M=max(M,(long)ex[2]-(long)ex[0]);maxd=max(maxd,(long)ex[2]);}}
  while((long)powers.size()<=maxd){Poly p;fmpq_mpoly_mul(p.p,powers.back().p,D.p,ctx);powers.push_back(move(p));}
  long count=0;
  for(auto&p:source){Poly result;
   for(int j=0;j<fmpq_mpoly_length(p.p,ctx);j++){
    fmpq_mpoly_get_term_exp_ui(ex,p.p,j,ctx);fmpq_mpoly_get_term_coeff_fmpq(c,p.p,j,ctx);long d=ex[2];ex[0]=ex[0]+M-d;ex[2]=0;
    Poly mon,term;fmpq_mpoly_push_term_fmpq_ui(mon.p,c,ex,ctx);fmpq_mpoly_mul(term.p,mon.p,powers[d].p,ctx);fmpq_mpoly_add(result.p,result.p,term.p,ctx);
   }count+=fmpq_mpoly_length(result.p,ctx);write(out,result);
  }meta<<row<<' '<<M<<' '<<maxd<<' '<<count<<'\n';meta.flush();allterms+=count;cerr<<"J0_ROW "<<row<<" a_clear "<<M<<" maxd "<<maxd<<" terms "<<count<<endl;
 }cerr<<"EXACT_J0_SUBSTITUTION_COMPLETE TERMS "<<allterms<<endl;fmpq_clear(c);
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();
}
