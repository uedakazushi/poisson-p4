// Exact finite polynomial-matrix certificate; no formal-tail extrapolation.
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly_factor.h>
#include <flint/ulong_extras.h>
#include <fstream>
#include <vector>
#include <chrono>
#include <iostream>
#include <cassert>
using namespace std;ulong prime;
struct P{nmod_poly_t p;P(){nmod_poly_init(p,prime);}P(const P&o):P(){nmod_poly_set(p,o.p);}P(P&&o):P(){nmod_poly_swap(p,o.p);}P&operator=(P o){nmod_poly_swap(p,o.p);return *this;}~P(){nmod_poly_clear(p);}};
using Row=vector<P>;vector<ulong>roots;
pair<int,int>lead(const Row&r){int d=-1,k=-1;for(int j=0;j<(int)r.size();j++){int e=nmod_poly_degree(r[j].p);if(e>=0&&e>=d){d=e;k=j;}}return{d,k};}
int strip(Row&r){int steps=0;for(auto a:roots)while(true){bool all=true,nonzero=false;for(auto&v:r){nonzero|=!nmod_poly_is_zero(v.p);if(nmod_poly_evaluate_nmod(v.p,a)){all=false;break;}}if(!all||!nonzero)break;for(auto&v:r){P q;assert(nmod_poly_div_root(q.p,v.p,a)==0);v=std::move(q);}steps++;}return steps;}
int main(int argc,char**argv){assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);int nr,n;in>>prime>>nr>>n;roots={(prime-6)%prime,(prime-12)*n_invmod(5,prime)%prime,(prime-119)*n_invmod(20,prime)%prime,(prime-4)*n_invmod(25,prime)%prime};vector<Row>basis(n);vector<int>degs(n,-1);int rank=0;long operations=0,stripped=0;P tmp,sh;auto st=chrono::steady_clock::now();auto sec=[&](){return chrono::duration<double>(chrono::steady_clock::now()-st).count();};
for(int i=0;i<nr;i++){Row row(n);int nz;in>>nz;for(int z=0;z<nz;z++){int j,len;in>>j>>len;for(int k=0;k<len;k++){ulong c;in>>c;nmod_poly_set_coeff_ui(row[j].p,k,c);}}assert(in);stripped+=strip(row);while(true){auto[d,k]=lead(row);if(k<0)break;if(basis[k].empty()){basis[k]=std::move(row);degs[k]=d;rank++;break;}if(d<degs[k]){swap(row,basis[k]);swap(d,degs[k]);}ulong lc=nmod_poly_get_coeff_ui(row[k].p,d),b=nmod_poly_get_coeff_ui(basis[k][k].p,degs[k]);ulong c=lc*n_invmod(b,prime)%prime;for(int j=0;j<n;j++){nmod_poly_scalar_mul_nmod(tmp.p,basis[k][j].p,c);nmod_poly_shift_left(sh.p,tmp.p,d-degs[k]);nmod_poly_sub(row[j].p,row[j].p,sh.p);}stripped+=strip(row);operations++;if(sec()>85){out<<"TIMEOUT row "<<i<<" rank "<<rank<<" operations "<<operations<<" strips "<<stripped<<'\n';return 2;}}
if(i%25==0){long ds=0;for(int d:degs)if(d>0)ds+=d;cerr<<"ROW "<<i<<" RANK "<<rank<<" DEGS "<<ds<<" OPS "<<operations<<" SEC "<<sec()<<'\n';}}

// Saturate only at the four explicitly inverted frame factors.
auto insert=[&](Row row){stripped+=strip(row);while(true){auto[d,k]=lead(row);if(k<0)break;assert(!basis[k].empty());if(d<degs[k]){swap(row,basis[k]);swap(d,degs[k]);}ulong lc=nmod_poly_get_coeff_ui(row[k].p,d),b=nmod_poly_get_coeff_ui(basis[k][k].p,degs[k]);ulong c=lc*n_invmod(b,prime)%prime;for(int j=0;j<n;j++){nmod_poly_scalar_mul_nmod(tmp.p,basis[k][j].p,c);nmod_poly_shift_left(sh.p,tmp.p,d-degs[k]);nmod_poly_sub(row[j].p,row[j].p,sh.p);}stripped+=strip(row);operations++;assert(sec()<85);}};
assert(rank==n);nmod_mat_t E,K;nmod_mat_init(E,n,n,prime);nmod_mat_init(K,n,n,prime);
for(auto root:roots)for(int loop=0;loop<100;loop++){
 for(int i=0;i<n;i++)for(int j=0;j<n;j++)nmod_mat_entry(E,j,i)=nmod_poly_evaluate_nmod(basis[i][j].p,root);
 int nul=nmod_mat_nullspace(K,E);long ds=0;for(int d:degs)ds+=d;cerr<<"SAT "<<root<<" LOOP "<<loop<<" NULL "<<nul<<" DEG "<<ds<<" SEC "<<sec()<<'\n';if(!nul)break;
 vector<Row>add;for(int k=0;k<nul;k++){Row row(n);for(int i=0;i<n;i++){ulong c=nmod_mat_entry(K,i,k);if(!c)continue;for(int j=0;j<n;j++){nmod_poly_scalar_mul_nmod(tmp.p,basis[i][j].p,c);nmod_poly_add(row[j].p,row[j].p,tmp.p);}}
 for(auto&v:row){P q;assert(!nmod_poly_div_root(q.p,v.p,root));v=std::move(q);}add.push_back(std::move(row));}
 for(auto&row:add)insert(std::move(row));assert(loop<99);
}
nmod_mat_clear(E);nmod_mat_clear(K);long ds=0;for(int d:degs)ds+=d;out<<"RANK "<<rank<<" DEGREE "<<ds<<" OPERATIONS "<<operations<<" STRIPS "<<stripped<<" SECONDS "<<sec()<<'\n';out.flush();
if(ds<1000){nmod_poly_mat_t B;nmod_poly_mat_init(B,n,n,prime);for(int i=0;i<n;i++)for(int j=0;j<n;j++)nmod_poly_set(nmod_poly_mat_entry(B,i,j),basis[i][j].p);P det;nmod_poly_mat_det_interpolate(det.p,B);assert(nmod_poly_degree(det.p)==ds);out<<"DET "<<ds;for(int k=0;k<=ds;k++)out<<' '<<nmod_poly_get_coeff_ui(det.p,k);out<<'\n';nmod_poly_factor_t ff;nmod_poly_factor_init(ff);nmod_poly_factor(ff,det.p);out<<"FACTORS "<<ff->num<<'\n';for(int i=0;i<ff->num;i++){out<<ff->exp[i]<<' '<<nmod_poly_degree(ff->p+i);for(int k=0;k<=nmod_poly_degree(ff->p+i);k++)out<<' '<<nmod_poly_get_coeff_ui(ff->p+i,k);out<<'\n';}nmod_poly_factor_clear(ff);nmod_poly_mat_clear(B);}
out<<"FINAL_SECONDS "<<sec()<<'\n';}
