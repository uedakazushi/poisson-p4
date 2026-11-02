#include "s61_binary_deeper_quartic_field_20261101.hpp"
#include <flint/nmod_poly_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_mat.h>
#include <flint/nmod.h>
#include <flint/ulong_extras.h>
#include <array>
#include <chrono>
#include <climits>
using Mon=array<int,9>;int nw[9]={1,1,1,1,1,2,1,2,1};
void mons(vector<Mon>&a,int d,int pos=0,Mon m={}){if(pos==9){if(!d)a.push_back(m);return;}for(int j=0;j*nw[pos]<=d;j++){m[pos]=j;mons(a,d-j*nw[pos],pos+1,m);}}
struct Uni{nmod_poly_struct p[1];Uni(){nmod_poly_init(p,32003);}Uni(const Uni&a):Uni(){nmod_poly_set(p,a.p);}Uni(Uni&&a):Uni(){nmod_poly_swap(p,a.p);}Uni&operator=(Uni a){nmod_poly_swap(p,a.p);return *this;}~Uni(){nmod_poly_clear(p);}};
Uni specialize(const fmpq_mpoly_t f,ulong e){Uni u;fmpq_t c;fmpq_init(c);ulong exp[2];nmod_t mod;nmod_init(&mod,32003);for(int i=0;i<fmpq_mpoly_length(f,ctx);i++){fmpq_mpoly_get_term_coeff_fmpq(c,f,i,ctx);fmpq_mpoly_get_term_exp_ui(exp,f,i,ctx);ulong den=fmpz_fdiv_ui(fmpq_denref(c),mod.n);assert(den);ulong v=nmod_mul(fmpz_fdiv_ui(fmpq_numref(c),mod.n),n_invmod(den,mod.n),mod);v=nmod_mul(v,nmod_pow_ui(e,exp[0],mod),mod);v=nmod_add(v,nmod_poly_get_coeff_ui(u.p,exp[1]),mod);nmod_poly_set_coeff_ui(u.p,exp[1],v);}fmpq_clear(c);return u;}
long max_assignment(const vector<vector<int>>&a){int n=a.size();vector<long>u(n+1),v(n+1);vector<int>p(n+1),way(n+1);for(int i=1;i<=n;i++){p[0]=i;int j0=0;vector<long>minv(n+1,LONG_MAX/4);vector<char>used(n+1);do{used[j0]=true;int i0=p[j0],j1=0;long delta=LONG_MAX/4;for(int j=1;j<=n;j++)if(!used[j]){long cur=-a[i0-1][j-1]-u[i0]-v[j];if(cur<minv[j]){minv[j]=cur;way[j]=j0;}if(minv[j]<delta){delta=minv[j];j1=j;}}for(int j=0;j<=n;j++)if(used[j]){u[p[j]]+=delta;v[j]-=delta;}else minv[j]-=delta;j0=j1;}while(p[j0]);do{int j1=way[j0];p[j0]=p[j1];j0=j1;}while(j0);}long total=0;for(int j=1;j<=n;j++){assert(a[p[j]-1][j-1]>-10000);total+=a[p[j]-1][j-1];}return total;}
int main(int argc,char**argv){assert(argc==4);fmpq_mpoly_ctx_init(ctx,2,ORD_LEX);{ifstream in(argv[1]);ofstream out(argv[3]);ulong ee=stoul(argv[2]);vector<Mon>mm;mons(mm,4);map<Mon,int>ix;for(int i=0;i<269;i++)ix[mm[i]]=i;vector<map<int,Uni>>A;vector<map<int,int>>Deg;int nr;in>>nr;for(int i=0;i<nr;i++){int id,w,n;in>>id>>w>>n;map<Mon,Uni>q;map<Mon,int>deg;for(int k=0;k<n;k++){Mon m;for(auto&a:m)in>>a;Rat f;readpoly(in,f.n);q.emplace(m,specialize(f.n,ee));deg[m]=fmpq_mpoly_degree_si(f.n,1,ctx);}vector<Mon>sh;mons(sh,4-w);for(auto&s:sh){map<int,Uni>row;map<int,int>drow;for(auto&[m,c]:q){Mon a=m;for(int j=0;j<9;j++)a[j]+=s[j];int j=ix.at(a);row.emplace(j,c);drow[j]=deg.at(m);}A.push_back(std::move(row));Deg.push_back(std::move(drow));}}assert(A.size()==390);auto start=chrono::steady_clock::now();
vector<vector<Uni>> basis(269);vector<int>degrees(269,-1);long steps=0;int rank=0;
auto leading=[&](const vector<Uni>&r){int deg=-1,pos=-1;for(int j=0;j<269;j++){int d=nmod_poly_degree(r[j].p);if(d>=0&&(d>deg||(d==deg&&j>pos))){deg=d;pos=j;}}return pair<int,int>(deg,pos);};
Uni tmp,shift;
for(int i=0;i<390;i++){
 vector<Uni>r(269);for(auto&[j,v]:A[i])r[j]=std::move(v);
 while(true){auto[deg,pos]=leading(r);if(pos<0)break;
  if(basis[pos].empty()){basis[pos]=std::move(r);degrees[pos]=deg;rank++;break;}
  if(deg<degrees[pos]){swap(r,basis[pos]);swap(deg,degrees[pos]);}
  ulong c=nmod_poly_get_coeff_ui(r[pos].p,deg)*n_invmod(nmod_poly_get_coeff_ui(basis[pos][pos].p,degrees[pos]),32003)%32003;
  for(int j=0;j<269;j++){nmod_poly_scalar_mul_nmod(tmp.p,basis[pos][j].p,c);nmod_poly_shift_left(shift.p,tmp.p,deg-degrees[pos]);nmod_poly_sub(r[j].p,r[j].p,shift.p);}
  steps++;if(steps%10000==0)cerr<<"STEPS "<<steps<<" ROW "<<i<<" RANK "<<rank<<" SEC "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';
  assert(chrono::duration<double>(chrono::steady_clock::now()-start).count()<150);
 }
}
assert(rank==269);long sd=0;for(int d:degrees)sd+=d;cerr<<"WEAK_POPOV ROWS "<<rank<<" DEGREE "<<sd<<" STEPS "<<steps<<" SEC "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';
nmod_mat_t special;nmod_mat_init(special,269,269,32003);for(int i=0;i<269;i++)for(int j=0;j<269;j++)nmod_mat_entry(special,i,j)=nmod_poly_evaluate_nmod(basis[i][j].p,32003-2162);cerr<<"AT_KNOWN_ROOT_RANK "<<nmod_mat_rank(special)<<" MIN_ROW_DEGREE "<<*min_element(degrees.begin(),degrees.end())<<" MAX_ROW_DEGREE "<<*max_element(degrees.begin(),degrees.end())<<'\n';nmod_mat_clear(special);
nmod_poly_mat_t M;nmod_poly_mat_init(M,269,269,32003);for(int i=0;i<269;i++)for(int j=0;j<269;j++)nmod_poly_set(nmod_poly_mat_entry(M,i,j),basis[i][j].p);Uni det;nmod_poly_mat_det_interpolate(det.p,M);assert(nmod_poly_degree(det.p)==sd);out<<ee<<" 32003 "<<rank<<' '<<sd<<' '<<steps<<'\n';out<<nmod_poly_length(det.p);for(int k=0;k<nmod_poly_length(det.p);k++)out<<' '<<nmod_poly_get_coeff_ui(det.p,k);out<<'\n';cerr<<"DONE SEC "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';nmod_poly_mat_clear(M);}fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
