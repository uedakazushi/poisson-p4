// Independent streaming QQ[H,L,s,q,t] source coefficient checker.
// No field arithmetic, pivots, row primitive division, or full output matrix.
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <gmp.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
using namespace std;
fmpq_mpoly_ctx_t ctx;long maxinter=0,products=0;
struct P{fmpq_mpoly_t x;P(){fmpq_mpoly_init(x,ctx);}P(const P&a):P(){fmpq_mpoly_set(x,a.x,ctx);}P(P&&a):P(){fmpq_mpoly_swap(x,a.x,ctx);}P&operator=(P a){fmpq_mpoly_swap(x,a.x,ctx);return *this;}~P(){fmpq_mpoly_clear(x,ctx);}};
void guard(const P&p){maxinter=max(maxinter,fmpq_mpoly_length(p.x,ctx));assert(maxinter<150000);}
P readsource(istream&in){P p;long n;in>>n;assert(in&&n>=0);fmpq_t q;fmpq_init(q);ulong e[5];for(long i=0;i<n;i++){string c;in>>c;for(auto&v:e)in>>v;assert(!fmpq_set_str(q,c.c_str(),10));fmpq_canonicalise(q);fmpq_mpoly_set_coeff_fmpq_ui(p.x,q,e,ctx);}fmpq_clear(q);guard(p);return p;}
P expected(){P p;string content;long n;cin>>content>>n;assert(cin&&n>=0);assert(!fmpq_set_str(p.x->content,content.c_str(),10));fmpq_canonicalise(p.x->content);fmpz_t c;fmpz_init(c);ulong e[5];for(long i=0;i<n;i++){string z;cin>>z;for(auto&v:e)cin>>v;assert(!fmpz_set_str(c,z.c_str(),10));assert(!fmpz_is_zero(c));fmpz_mpoly_set_coeff_fmpz_ui(p.x->zpoly,c,e,ctx->zctx);}fmpz_clear(c);fmpq_mpoly_reduce(p.x,ctx);assert(fmpq_mpoly_length(p.x,ctx)==n);return p;}
P add(const P&a,const P&b){P p;fmpq_mpoly_add(p.x,a.x,b.x,ctx);guard(p);return p;}
P mul(const P&a,const P&b){P p;fmpq_mpoly_mul(p.x,a.x,b.x,ctx);guard(p);products++;return p;}
int main(int argc,char**argv){assert(argc==2);auto start=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,5,ORD_DEGLEX);long comparisons=0,total=0,maxentry=0,maxrow=0,terms=0;{
  ifstream in(argv[1]);int nr,nq;in>>nr>>nq;assert(nr==134&&nq==45);P den=readsource(in);vector<P>C(45),d(45);vector<vector<P>>W(45,vector<P>(2));for(auto&p:C)p=readsource(in);for(auto&p:d)p=readsource(in);for(auto&row:W)for(auto&p:row)p=readsource(in);
  for(int row=0;row<134;row++){
    vector<P>A(5);int groups;in>>groups;assert(in&&groups>=0);for(int k=0;k<groups;k++){int degree;in>>degree;assert(degree>=0&&degree<=2);vector<int>m(degree);for(auto&i:m){in>>i;assert(i>=0&&i<45);}P co=readsource(in);terms+=fmpq_mpoly_length(co.x,ctx);
      if(degree==0)A[0]=add(A[0],co);
      else if(degree==1){int a=m[0];A[0]=add(A[0],mul(co,C[a]));for(int h=0;h<2;h++)A[1+h]=add(A[1+h],mul(co,W[a][h]));}
      else{int a=m[0],b=m[1];A[0]=add(A[0],mul(co,mul(C[a],C[b])));for(int h=0;h<2;h++){
        A[1+h]=add(A[1+h],mul(co,add(mul(W[a][h],C[b]),mul(C[a],W[b][h]))));
        A[3+h]=add(A[3+h],mul(co,add(mul(W[a][h],d[b]),mul(d[a],W[b][h]))));
      }}
    }A[0]=mul(den,A[0]);long rowterms=0;
    for(int j=0;j<5;j++){P E=expected();if(!fmpq_mpoly_equal(A[j].x,E.x,ctx)){cerr<<"SOURCE_MISMATCH row="<<row<<" col="<<j<<endl;return 2;}long n=fmpq_mpoly_length(A[j].x,ctx);rowterms+=n;maxentry=max(maxentry,n);total+=n;comparisons++;}
    assert(rowterms<=10415);maxrow=max(maxrow,rowterms);
  }
  string extra;assert(in&&!(in>>extra));assert(cin&&!(cin>>extra));assert(comparisons==670&&total==692953&&maxentry==1685);
}fmpq_mpoly_ctx_clear(ctx);flint_cleanup();cout<<"{\"all_passed\":true,\"source_rows\":134,\"literal_QQ_polynomial_equalities\":"<<comparisons<<",\"all_expected_terms_compared\":"<<total<<",\"max_entry_terms\":"<<maxentry<<",\"max_row_terms\":"<<maxrow<<",\"source_coefficient_terms\":"<<terms<<",\"polynomial_products\":"<<products<<",\"max_intermediate_terms\":"<<maxinter<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";}
