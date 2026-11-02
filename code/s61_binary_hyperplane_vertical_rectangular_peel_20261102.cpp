#include <flint/fmpq_mpoly.h>
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <chrono>
using namespace std;
fmpq_mpoly_ctx_t ctx;
struct Poly {
 fmpq_mpoly_t p;
 Poly(){fmpq_mpoly_init(p,ctx);}
 Poly(const Poly&a):Poly(){fmpq_mpoly_set(p,a.p,ctx);}
 Poly(Poly&&a):Poly(){fmpq_mpoly_swap(p,a.p,ctx);}
 Poly&operator=(Poly a){fmpq_mpoly_swap(p,a.p,ctx);return *this;}
 ~Poly(){fmpq_mpoly_clear(p,ctx);}
};
void read(istream&in,Poly&p,int dim){
 int n;in>>n;assert(n>=0);ulong ex[3]={};string s;fmpq_t q;fmpq_init(q);
 for(int i=0;i<n;i++){
  in>>s;assert(!fmpq_set_str(q,s.c_str(),10));fmpq_canonicalise(q);
  assert(fmpz_fdiv_ui(fmpq_denref(q),32003));
  for(int j=0;j<dim;j++)in>>ex[j];assert(ex[0]==0);
  fmpq_mpoly_push_term_fmpq_ui(p.p,q,ex,ctx);
 }
 assert(in);fmpq_mpoly_sort_terms(p.p,ctx);fmpq_mpoly_combine_like_terms(p.p,ctx);fmpq_clear(q);
}
void write(ostream&out,const Poly&p,int dim){
 ulong ex[3];fmpq_t c;fmpq_init(c);int n=fmpq_mpoly_length(p.p,ctx);out<<n;
 for(int i=0;i<n;i++){
  fmpq_mpoly_get_term_coeff_fmpq(c,p.p,i,ctx);fmpq_mpoly_get_term_exp_ui(ex,p.p,i,ctx);
  char*s=fmpq_get_str(nullptr,10,c);out<<' '<<s;flint_free(s);
  for(int j=0;j<dim;j++)out<<' '<<ex[j];
 }
 out<<'\n';fmpq_clear(c);
}
Poly truncate(const Poly&p,int degree,bool exact=false){
 Poly r;fmpq_t c;fmpq_init(c);ulong ex[3];
 for(int i=0;i<fmpq_mpoly_length(p.p,ctx);i++){
  fmpq_mpoly_get_term_exp_ui(ex,p.p,i,ctx);
  if(exact)assert(ex[2]<=(ulong)degree);
  if(ex[2]>(ulong)degree)continue;
  fmpq_mpoly_get_term_coeff_fmpq(c,p.p,i,ctx);fmpq_mpoly_push_term_fmpq_ui(r.p,c,ex,ctx);
 }
 fmpq_mpoly_sort_terms(r.p,ctx);fmpq_clear(c);return r;
}
int main(int argc,char**argv){
 assert(argc==7);auto start=chrono::steady_clock::now();
 auto sec=[&](){return chrono::duration<double>(chrono::steady_clock::now()-start).count();};
 int precision=stoi(argv[6]);assert(precision==1||precision==2);
 ifstream in(argv[1]),lead(argv[2]),ker(argv[3]);ofstream out(argv[4]),o0(argv[5]);assert(in&&lead&&ker&&out&&o0);
 fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);
 {
  int nr,nc,lr,lc,r,k;in>>nr>>nc;lead>>lr>>lc>>r>>k;
  assert(nr==lr&&nc==lc&&nr>=nc&&r+k==nc);
  vector<int>rows(r),cols(r),free(k);for(int&i:rows)lead>>i;for(int&i:cols)lead>>i;for(int&i:free)lead>>i;
  vector<vector<Poly>>A(nr,vector<Poly>(nc)),B(nr,vector<Poly>(nc));
  for(auto&row:A)for(auto&p:row){read(in,p,3);truncate(p,precision,true);}
  string extra;assert(!(in>>extra));
  for(int i=0;i<nr;i++)for(int j=0;j<r;j++)B[i][j]=truncate(A[i][cols[j]],precision-1);
  Poly delta;fmpq_mpoly_gen(delta.p,2,ctx);long terms=0;
  for(int c=0;c<k;c++){
   string s;int t;ker>>s>>t;assert(s=="KERNEL"&&t==c);Poly den;read(ker,den,2);
   vector<Poly>K(nc);for(auto&p:K)read(ker,p,2);
   for(int j=0;j<k;j++)assert(j==c?fmpq_mpoly_equal(K[free[j]].p,den.p,ctx):fmpq_mpoly_is_zero(K[free[j]].p,ctx));
   assert(!fmpq_mpoly_is_zero(den.p,ctx));
   for(int i=0;i<nr;i++){
    Poly z,tmp,div;
    for(int j=0;j<nc;j++)if(!fmpq_mpoly_is_zero(K[j].p,ctx)&&!fmpq_mpoly_is_zero(A[i][j].p,ctx)){
     fmpq_mpoly_mul(tmp.p,A[i][j].p,K[j].p,ctx);fmpq_mpoly_add(z.p,z.p,tmp.p,ctx);
    }
    assert(fmpq_mpoly_divides(div.p,z.p,delta.p,ctx));
    B[i][r+c]=truncate(div,precision-1,true);terms+=fmpq_mpoly_length(div.p,ctx);
   }
   cerr<<"COLUMN "<<c<<" KERNEL_TERMS "<<terms<<" SEC "<<sec()<<'\n';
   assert(sec()<60&&terms<2000000);
  }
  assert(!(ker>>extra));
  out<<nr<<' '<<nc<<'\n';o0<<nr<<' '<<nc<<'\n';long allterms=0,leadingterms=0;
  for(auto&row:B)for(auto&p:row){
   allterms+=fmpq_mpoly_length(p.p,ctx);write(out,p,3);
   Poly c=truncate(p,0);leadingterms+=fmpq_mpoly_length(c.p,ctx);write(o0,c,2);
  }
  assert(allterms<2000000);out.flush();o0.flush();assert(out.good()&&o0.good());
  cerr<<"PEEL_PASS "<<k<<" PRECISION "<<precision<<" TO "<<precision-1<<" TERMS "<<allterms<<" LEADING "<<leadingterms<<" SEC "<<sec()<<'\n';
 }
 fmpq_mpoly_ctx_clear(ctx);flint_cleanup();
}
