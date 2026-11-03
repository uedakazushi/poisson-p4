// Primitive integer rows, then columns, by NONZERO QQ SCALARS ONLY.
// CLI input output factors metadata(0/1). Every scalar operation is reversed
// and compared with the full previous polynomial; no core polynomial is divided.
#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;fmpq_mpoly_ctx_t ctx;
struct Poly{fmpq_mpoly_t p;Poly(){fmpq_mpoly_init(p,ctx);}Poly(const Poly&a):Poly(){fmpq_mpoly_set(p,a.p,ctx);}Poly(Poly&&a):Poly(){fmpq_mpoly_swap(p,a.p,ctx);}Poly&operator=(Poly a){fmpq_mpoly_swap(p,a.p,ctx);return *this;}~Poly(){fmpq_mpoly_clear(p,ctx);}};
void rd(istream&in,Poly&p){long n;in>>n;assert(in&&n>=0);fmpq_t q;fmpq_init(q);ulong ex[3];for(long t=0;t<n;t++){string c;in>>c>>ex[0]>>ex[1]>>ex[2];assert(in&&!fmpq_set_str(q,c.c_str(),10));fmpq_canonicalise(q);fmpq_mpoly_push_term_fmpq_ui(p.p,q,ex,ctx);}fmpq_mpoly_sort_terms(p.p,ctx);fmpq_mpoly_combine_like_terms(p.p,ctx);fmpq_clear(q);}
void wr(ostream&out,const Poly&p){fmpq_t q;fmpq_init(q);ulong ex[3];long n=fmpq_mpoly_length(p.p,ctx);out<<n;for(long t=0;t<n;t++){fmpq_mpoly_get_term_coeff_fmpq(q,p.p,t,ctx);fmpq_mpoly_get_term_exp_ui(ex,p.p,t,ctx);char*c=fmpq_get_str(nullptr,10,q);out<<' '<<c<<' '<<ex[0]<<' '<<ex[1]<<' '<<ex[2];flint_free(c);}out<<'\n';fmpq_clear(q);}
void primitive(vector<Poly*>&ps,fmpq_t factor){
 fmpz_t den,content,coefficient,quot;fmpz_init(den);fmpz_init(content);fmpz_init(coefficient);fmpz_init(quot);fmpz_one(den);fmpq_t q;fmpq_init(q);bool nonzero=false;
 for(auto*p:ps)for(long t=0;t<fmpq_mpoly_length(p->p,ctx);t++){nonzero=true;fmpq_mpoly_get_term_coeff_fmpq(q,p->p,t,ctx);fmpz_lcm(den,den,fmpq_denref(q));}
 if(!nonzero)fmpq_one(factor);else{
  for(auto*p:ps)for(long t=0;t<fmpq_mpoly_length(p->p,ctx);t++){fmpq_mpoly_get_term_coeff_fmpq(q,p->p,t,ctx);fmpz_divexact(quot,den,fmpq_denref(q));fmpz_mul(coefficient,fmpq_numref(q),quot);fmpz_gcd(content,content,coefficient);}
  assert(!fmpz_is_zero(content));fmpq_set_fmpz_frac(factor,den,content);assert(!fmpq_is_zero(factor));
 }
 for(auto*p:ps){Poly before=*p,reverse;fmpq_mpoly_scalar_mul_fmpq(p->p,p->p,factor,ctx);fmpq_mpoly_scalar_div_fmpq(reverse.p,p->p,factor,ctx);assert(fmpq_mpoly_equal(reverse.p,before.p,ctx));for(long t=0;t<fmpq_mpoly_length(p->p,ctx);t++){fmpq_mpoly_get_term_coeff_fmpq(q,p->p,t,ctx);assert(fmpz_is_one(fmpq_denref(q)));}}
 fmpq_clear(q);fmpz_clear(den);fmpz_clear(content);fmpz_clear(coefficient);fmpz_clear(quot);
}
int main(int argc,char**argv){assert(argc==5);fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);{
 ifstream in(argv[1]);ofstream out(argv[2]),factors(argv[3]);int n,m,k;in>>n>>m>>k;bool metadata=stoi(argv[4]);vector<int>used(k),ids(n);if(metadata){for(int&v:used)in>>v;for(int&v:ids)in>>v;}
 vector<vector<Poly>>A(n,vector<Poly>(m));long terms=0;for(auto&row:A)for(auto&p:row){rd(in,p);terms+=fmpq_mpoly_length(p.p,ctx);}factors<<n<<' '<<m<<'\n';fmpq_t q;fmpq_init(q);
 for(int i=0;i<n;i++){vector<Poly*>ps;for(int j=0;j<m;j++)ps.push_back(&A[i][j]);primitive(ps,q);char*c=fmpq_get_str(nullptr,10,q);factors<<"ROW "<<i<<' '<<c<<'\n';flint_free(c);cerr<<"EXACT_QQ_ROW_SCALAR "<<i<<endl;}
 for(int j=0;j<m;j++){vector<Poly*>ps;for(int i=0;i<n;i++)ps.push_back(&A[i][j]);primitive(ps,q);char*c=fmpq_get_str(nullptr,10,q);factors<<"COLUMN "<<j<<' '<<c<<'\n';flint_free(c);}
 out<<n<<' '<<m<<' '<<k<<'\n';if(metadata){for(int v:used)out<<v<<' ';out<<'\n';for(int v:ids)out<<v<<' ';out<<'\n';}for(auto&row:A)for(auto&p:row)wr(out,p);
 cerr<<"EXACT_QQ_DIAGONAL_SCALING_PASS TERMS "<<terms<<" SHAPE "<<n<<' '<<m<<endl;fmpq_clear(q);
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();
}
