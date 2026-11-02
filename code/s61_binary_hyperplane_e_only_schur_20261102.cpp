#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/nmod_mat.h>
#include <flint/ulong_extras.h>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
using namespace std;
struct RM{fmpz_mat_t A;fmpz_t D;RM(int r=0,int c=0){fmpz_mat_init(A,r,c);fmpz_init_set_ui(D,1);}RM(const RM&o){fmpz_mat_init_set(A,o.A);fmpz_init_set(D,o.D);}RM(RM&&o):RM(){fmpz_mat_swap(A,o.A);fmpz_swap(D,o.D);}RM&operator=(RM o){fmpz_mat_swap(A,o.A);fmpz_swap(D,o.D);return *this;}~RM(){fmpz_clear(D);fmpz_mat_clear(A);}int r()const{return fmpz_mat_nrows(A);}int c()const{return fmpz_mat_ncols(A);}void norm(){fmpz_t g;fmpz_init(g);fmpz_mat_content(g,A);fmpz_gcd(g,g,D);if(!fmpz_is_one(g)){fmpz_mat_scalar_divexact_fmpz(A,A,g);fmpz_divexact(D,D,g);}if(fmpz_sgn(D)<0){fmpz_neg(D,D);fmpz_mat_neg(A,A);}fmpz_clear(g);assert(fmpz_fdiv_ui(D,32003));}};
RM block(const RM&a,int r0,int r1,int c0,int c1){RM b(r1-r0,c1-c0);fmpz_mat_t W;fmpz_mat_window_init(W,a.A,r0,c0,r1,c1);fmpz_mat_set(b.A,W);fmpz_mat_window_clear(W);fmpz_set(b.D,a.D);b.norm();return b;}
RM mul(const RM&a,const RM&b){assert(a.c()==b.r());RM z(a.r(),b.c());fmpz_mat_mul(z.A,a.A,b.A);fmpz_mul(z.D,a.D,b.D);z.norm();return z;}
RM sub(const RM&a,const RM&b){assert(a.r()==b.r()&&a.c()==b.c());RM x(a),y(b);fmpz_t g,f,t;for(auto z:{g,f,t})fmpz_init(z);fmpz_gcd(g,a.D,b.D);fmpz_divexact(f,b.D,g);fmpz_divexact(t,a.D,g);fmpz_mat_scalar_mul_fmpz(x.A,x.A,f);fmpz_mat_scalar_mul_fmpz(y.A,y.A,t);fmpz_mat_sub(x.A,x.A,y.A);fmpz_mul(x.D,a.D,f);x.norm();for(auto z:{g,f,t})fmpz_clear(z);return x;}
RM inv(const RM&a){assert(a.r()==a.c());RM x(a.r(),a.c());assert(fmpz_mat_inv(x.A,x.D,a.A));fmpz_mat_scalar_mul_fmpz(x.A,x.A,a.D);x.norm();return x;}
unsigned long detmod(const RM&a){const unsigned long p=32003;nmod_mat_t x;nmod_mat_init(x,a.r(),a.c(),p);for(int i=0;i<a.r();i++)for(int j=0;j<a.c();j++)nmod_mat_entry(x,i,j)=fmpz_fdiv_ui(fmpz_mat_entry(a.A,i,j),p);unsigned long v=nmod_mat_det(x);nmod_mat_clear(x);unsigned long u=n_invmod(fmpz_fdiv_ui(a.D,p),p);for(int i=0;i<a.r();i++)v=v*u%p;return v;}
int main(int argc,char**argv){assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);int n,prec,k;in>>n>>prec>>k;assert(k==prec);vector<int>ranks(k);for(int&i:ranks)in>>i;vector<RM>M;for(int j=0;j<=prec;j++)M.emplace_back(n,n);for(int i=0;i<n;i++)for(int j=0;j<n;j++){int ct;in>>ct;for(int l=0;l<ct;l++){int d;string s;in>>d>>s;assert(d>=0&&d<=prec&&!fmpz_set_str(fmpz_mat_entry(M[d].A,i,j),s.c_str(),10));}}string extra;assert(!(in>>extra));auto start=chrono::steady_clock::now();int val=0;unsigned long lc=1;
for(auto r:ranks){assert(0<r&&r<n);auto stamp=[&](){return chrono::duration<double>(chrono::steady_clock::now()-start).count();};cerr<<"BEGIN "<<n<<" PIV "<<r<<" PREC "<<prec<<" SEC "<<stamp()<<"\n";vector<RM>A,B,C,E,Y,N;for(auto&v:M){A.push_back(block(v,0,r,0,r));B.push_back(block(v,0,r,r,n));C.push_back(block(v,r,n,0,r));E.push_back(block(v,r,n,r,n));}M.clear();auto u=detmod(A[0]);assert(u);lc=lc*u%32003;RM I=inv(A[0]);cerr<<"INVERSE "<<r<<" DENBITS "<<fmpz_bits(I.D)<<" SEC "<<stamp()<<"\n";
for(int j=0;j<=prec;j++){RM rhs=B[j];for(int i=1;i<=j;i++)rhs=sub(rhs,mul(A[i],Y[j-i]));Y.push_back(mul(I,rhs));RM h=E[j];for(int i=0;i<=j;i++)h=sub(h,mul(C[i],Y[j-i]));if(!j)assert(fmpz_mat_is_zero(h.A));else N.push_back(std::move(h));}
val+=n-r;out<<"STAGE "<<n<<" RANK "<<r<<" VAL "<<val<<" PIVOT_UNIT "<<u<<"\n";out.flush();cerr<<"DONE "<<r<<" VAL "<<val<<" SEC "<<stamp()<<"\n";M=std::move(N);n-=r;prec--;}
assert(n==1&&prec==0&&M.size()==1);auto v=detmod(M[0]);assert(v);lc=lc*v%32003;out<<"TOTAL "<<val<<" LAST "<<v<<" LC "<<lc<<"\n";out.flush();assert(out.good());}
