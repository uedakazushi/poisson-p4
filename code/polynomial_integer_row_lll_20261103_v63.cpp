// Exact constant integer row operations for an integer sparse polynomial matrix.
// CLI: input output transform. Metadata is preserved; a saved U with det +/-1
// certifies cokernel equivalence. LLL is only a basis-choice heuristic: full
// Gram congruence and exact polynomial row construction are checked afterwards.
#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_lll.h>
#include <cassert>
#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std; fmpq_mpoly_ctx_t ctx;
struct Poly{fmpq_mpoly_t p;Poly(){fmpq_mpoly_init(p,ctx);}Poly(const Poly&a):Poly(){fmpq_mpoly_set(p,a.p,ctx);}Poly(Poly&&a):Poly(){fmpq_mpoly_swap(p,a.p,ctx);}Poly&operator=(Poly a){fmpq_mpoly_swap(p,a.p,ctx);return *this;}~Poly(){fmpq_mpoly_clear(p,ctx);}};
void rd(istream&in,Poly&p){long n;in>>n;assert(in&&n>=0);fmpq_t q;fmpq_init(q);ulong ex[3];for(long t=0;t<n;t++){string s;in>>s>>ex[0]>>ex[1]>>ex[2];assert(in&&!fmpq_set_str(q,s.c_str(),10));fmpq_canonicalise(q);assert(fmpz_is_one(fmpq_denref(q)));fmpq_mpoly_push_term_fmpq_ui(p.p,q,ex,ctx);}fmpq_mpoly_sort_terms(p.p,ctx);fmpq_mpoly_combine_like_terms(p.p,ctx);fmpq_clear(q);}
void wr(ostream&out,const Poly&p){fmpq_t q;fmpq_init(q);ulong ex[3];long n=fmpq_mpoly_length(p.p,ctx);out<<n;for(long t=0;t<n;t++){fmpq_mpoly_get_term_coeff_fmpq(q,p.p,t,ctx);fmpq_mpoly_get_term_exp_ui(ex,p.p,t,ctx);char*c=fmpq_get_str(nullptr,10,q);out<<' '<<c<<' '<<ex[0]<<' '<<ex[1]<<' '<<ex[2];flint_free(c);}out<<'\n';fmpq_clear(q);}
void dot(fmpz_t out,const Poly&a,const Poly&b){
 fmpz_zero(out);long i=0,j=0,na=fmpq_mpoly_length(a.p,ctx),nb=fmpq_mpoly_length(b.p,ctx);array<ulong,3>x,y;fmpq_t q,r;fmpq_init(q);fmpq_init(r);
 while(i<na&&j<nb){fmpq_mpoly_get_term_exp_ui(x.data(),a.p,i,ctx);fmpq_mpoly_get_term_exp_ui(y.data(),b.p,j,ctx);if(x>y)i++;else if(x<y)j++;else{fmpq_mpoly_get_term_coeff_fmpq(q,a.p,i++,ctx);fmpq_mpoly_get_term_coeff_fmpq(r,b.p,j++,ctx);assert(fmpz_is_one(fmpq_denref(q))&&fmpz_is_one(fmpq_denref(r)));fmpz_addmul(out,fmpq_numref(q),fmpq_numref(r));}}
 fmpq_clear(q);fmpq_clear(r);
}
pair<long,long long>bits(const vector<vector<Poly>>&A){long maxbits=0;long long s=0;fmpq_t q;fmpq_init(q);for(auto&row:A)for(auto&p:row)for(long k=0;k<fmpq_mpoly_length(p.p,ctx);k++){fmpq_mpoly_get_term_coeff_fmpq(q,p.p,k,ctx);long b=fmpz_bits(fmpq_numref(q));maxbits=max(maxbits,b);s+=b;}fmpq_clear(q);return{maxbits,s};}
int main(int argc,char**argv){assert(argc==4);fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);{
 ifstream in(argv[1]);ofstream out(argv[2]),trace(argv[3]);int n,m,k;in>>n>>m>>k;assert(n>0&&n<=57);vector<int>used(k),ids(n);for(int&x:used)in>>x;for(int&x:ids)in>>x;vector<vector<Poly>>A(n,vector<Poly>(m));for(auto&row:A)for(auto&p:row)rd(in,p);auto before=bits(A);cerr<<"READ_INTEGER_MATRIX "<<n<<' '<<m<<" MAXBITS "<<before.first<<" SUMBITS "<<before.second<<endl;
 fmpz_mat_t G,H,U,T,UT,C;fmpz_mat_init(G,n,n);fmpz_mat_init(H,n,n);fmpz_mat_init(U,n,n);fmpz_mat_init(T,n,n);fmpz_mat_init(UT,n,n);fmpz_mat_init(C,n,n);fmpz_t z,det;fmpz_init(z);fmpz_init(det);
 for(int i=0;i<n;i++)for(int j=0;j<=i;j++){for(int col=0;col<m;col++){dot(z,A[i][col],A[j][col]);fmpz_add(fmpz_mat_entry(G,i,j),fmpz_mat_entry(G,i,j),z);}fmpz_set(fmpz_mat_entry(G,j,i),fmpz_mat_entry(G,i,j));}cerr<<"FULL_COEFFICIENT_GRAM_READY\n";
 fmpz_mat_set(H,G);fmpz_mat_one(U);fmpz_lll_t fl;fmpz_lll_context_init(fl,0.75,0.51,GRAM,EXACT);fmpz_lll(H,U,fl);fmpz_mat_det(det,U);assert(fmpz_is_pm1(det));fmpz_mat_transpose(UT,U);fmpz_mat_mul(T,U,G);fmpz_mat_mul(C,T,UT);assert(fmpz_mat_equal(C,H));cerr<<"EXACT_UNIMODULAR_AND_GRAM_PASS\n";
 vector<vector<Poly>>B(n,vector<Poly>(m));fmpq_t scalar;fmpq_init(scalar);Poly tmp;
 for(int i=0;i<n;i++){for(int j=0;j<m;j++)for(int h=0;h<n;h++){if(fmpz_is_zero(fmpz_mat_entry(U,i,h)))continue;fmpq_set_fmpz(scalar,fmpz_mat_entry(U,i,h));fmpq_mpoly_scalar_mul_fmpq(tmp.p,A[h][j].p,scalar,ctx);fmpq_mpoly_add(B[i][j].p,B[i][j].p,tmp.p,ctx);}cerr<<"EXACT_ROW_EXPORTED "<<i<<endl;}
 // Independently recompute the transformed full Gram, not merely its 20x20
 // integer congruence. This also detects coefficient-row mapping mistakes.
 for(int i=0;i<n;i++)for(int j=0;j<=i;j++){fmpz_zero(det);for(int col=0;col<m;col++){dot(z,B[i][col],B[j][col]);fmpz_add(det,det,z);}assert(fmpz_equal(det,fmpz_mat_entry(H,i,j)));}
 trace<<n<<'\n';for(int i=0;i<n;i++){for(int j=0;j<n;j++){char*s=fmpz_get_str(nullptr,10,fmpz_mat_entry(U,i,j));trace<<s<<' ';flint_free(s);}trace<<'\n';}
 out<<n<<' '<<m<<' '<<k<<'\n';for(int v:used)out<<v<<' ';out<<'\n';for(int v:ids)out<<v<<' ';out<<'\n';for(auto&row:B)for(auto&p:row)wr(out,p);auto after=bits(B);cerr<<"EXACT_FULL_POLYNOMIAL_UNIMODULAR_PASS BEFORE_MAXBITS "<<before.first<<" AFTER_MAXBITS "<<after.first<<" BEFORE_SUMBITS "<<before.second<<" AFTER_SUMBITS "<<after.second<<endl;
 fmpq_clear(scalar);fmpz_clear(z);fmpz_clear(det);fmpz_mat_clear(G);fmpz_mat_clear(H);fmpz_mat_clear(U);fmpz_mat_clear(T);fmpz_mat_clear(UT);fmpz_mat_clear(C);
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
