#include <cassert>
#include <fstream>
#include <iostream>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
using namespace std;
int main(int argc,char**argv){assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);int m,n,k;in>>m>>n>>k;assert(m==142&&n==111&&k==223);fmpq_mat_t A,R,H,T,P,V,W,Q,rankmat;fmpq_mat_init(A,m,n+k);fmpq_mat_init(R,m,n+k);
for(int i=0;i<m;i++)for(int j=0;j<n+k;j++){string z;in>>z;assert(fmpq_set_str(fmpq_mat_entry(A,i,j),z.c_str(),10)==0);fmpq_canonicalise(fmpq_mat_entry(A,i,j));}
assert(in);assert(fmpq_mat_rref(R,A)==n);for(int i=0;i<m;i++)for(int j=0;j<n;j++){if(i==j)assert(fmpq_is_one(fmpq_mat_entry(R,i,j)));else assert(fmpq_is_zero(fmpq_mat_entry(R,i,j)));}
for(int i=n;i<m;i++)for(int j=n;j<n+k;j++)assert(fmpq_is_zero(fmpq_mat_entry(R,i,j)));
fmpq_mat_init(H,n,n);fmpq_mat_init(T,n,n);fmpq_mat_init(P,n,n);fmpq_mat_init(Q,n,n);fmpq_mat_init(rankmat,n,n);fmpq_mat_init(V,n,1);fmpq_mat_init(W,n,1);
for(int i=0;i<n;i++){for(int j=0;j<n;j++){fmpq_set(fmpq_mat_entry(H,i,j),fmpq_mat_entry(R,i,n+j));fmpq_set(fmpq_mat_entry(T,i,j),fmpq_mat_entry(R,i,2*n+j));}fmpq_set(fmpq_mat_entry(V,i,0),fmpq_mat_entry(R,i,3*n));}
assert(fmpq_mat_is_zero(T));assert(!fmpq_mat_is_zero(V));fmpq_mat_mul(W,H,V);assert(!fmpq_mat_is_zero(W));fmpq_mat_mul(V,H,W);assert(fmpq_mat_is_zero(V));
fmpq_mat_one(P);cout<<"H ranks";for(int r=0;r<7;r++){cout<<" "<<fmpq_mat_rref(rankmat,P);fmpq_mat_mul(Q,P,H);fmpq_mat_set(P,Q);if(fmpq_mat_is_zero(P)){cout<<" 0";break;}}cout<<"\nT zero; t nonzero; Ht nonzero; H2t zero\n";
out<<n<<" "<<k<<"\n";for(int i=0;i<n;i++)for(int j=n;j<n+k;j++){char*z=fmpq_get_str(nullptr,10,fmpq_mat_entry(R,i,j));out<<z<<"\n";flint_free(z);}assert(out);}
