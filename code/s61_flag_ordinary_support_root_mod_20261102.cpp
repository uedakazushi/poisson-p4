// Exact finite polynomial-matrix certificate; no formal-tail extrapolation.
#include <flint/fmpq.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_mat.h>
#include <flint/ulong_extras.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <cassert>
#include <chrono>
using namespace std;const ulong P=32003;
struct Poly{nmod_poly_t p;Poly(){nmod_poly_init(p,P);}Poly(const Poly&a):Poly(){nmod_poly_set(p,a.p);}Poly(Poly&&a):Poly(){nmod_poly_swap(p,a.p);}Poly&operator=(Poly a){nmod_poly_swap(p,a.p);return *this;}~Poly(){nmod_poly_clear(p);}};
int main(int argc,char**argv){assert(argc==4);ulong root=stoul(argv[3]);ifstream in(argv[1]);ofstream out(argv[2]);int nr,n;in>>nr>>n;assert(nr==n);vector<vector<Poly>>A(n,vector<Poly>(n));vector<int>rowdegrees(n);fmpq_t q;fmpq_init(q);long bound=0;for(int i=0;i<n;i++){int nz;in>>nz;for(int k=0;k<nz;k++){int j,len;in>>j>>len;rowdegrees[i]=max(rowdegrees[i],len-1);for(int e=0;e<len;e++){string s;in>>s;assert(!fmpq_set_str(q,s.c_str(),10));fmpq_canonicalise(q);ulong den=fmpz_fdiv_ui(fmpq_denref(q),P);assert(den);ulong v=fmpz_fdiv_ui(fmpq_numref(q),P)*n_invmod(den,P)%P;nmod_poly_set_coeff_ui(A[i][j].p,e,v);}}for(int j=0;j<n;j++){Poly b;nmod_poly_taylor_shift(b.p,A[i][j].p,root);A[i][j]=std::move(b);}bound+=rowdegrees[i];}fmpq_clear(q);nmod_mat_t E,R;nmod_mat_init(E,n,n,P);nmod_mat_init(R,n,n,P);long valuation=0;auto st=chrono::steady_clock::now();for(int step=0;step<1000;step++){for(int i=0;i<n;i++)for(int j=0;j<n;j++)nmod_mat_entry(E,i,j)=nmod_poly_get_coeff_ui(A[i][j].p,0);nmod_mat_set(R,E);int rank=nmod_mat_rref(R);double sec=chrono::duration<double>(chrono::steady_clock::now()-st).count();out<<"STEP "<<step<<" RANK "<<rank<<" VALUATION "<<valuation<<" SEC "<<sec<<'\n';out.flush();cerr<<"STEP "<<step<<" RANK "<<rank<<" VAL "<<valuation<<" SEC "<<sec<<'\n';if(rank==n){out<<"DEGREE "<<bound-valuation<<" LEADING "<<nmod_mat_det(E)<<" BOUND "<<bound<<" VALUATION "<<valuation<<'\n';break;}vector<int>piv;vector<bool>used(n);for(int i=0;i<rank;i++){int j=0;while(!nmod_mat_entry(R,i,j))j++;piv.push_back(j);used[j]=true;}for(int j=0;j<n;j++)if(!used[j])for(int i=0;i<n;i++){Poly b=A[i][j],z;for(int k=0;k<rank;k++){ulong c=nmod_mat_entry(R,k,j);if(c){nmod_poly_scalar_mul_nmod(z.p,A[i][piv[k]].p,c);nmod_poly_sub(b.p,b.p,z.p);}}assert(!nmod_poly_get_coeff_ui(b.p,0));nmod_poly_shift_right(A[i][j].p,b.p,1);}valuation+=n-rank;if(sec>55){out<<"BOUNDED_STOP\n";return 2;}}
nmod_mat_clear(E);nmod_mat_clear(R);}
