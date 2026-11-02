#include "s61_binary_univariate_field_20261102.hpp"
#include <flint/fmpq_vec.h>
#include <flint/fmpq_mat.h>
#include <array>
#include <chrono>
#include <numeric>
#include <algorithm>
int main(int argc,char**argv){
 assert(argc==7);ifstream mi(argv[1]),di(argv[2]);string prefix=argv[3];int bound=stoi(argv[4]);ofstream out(argv[5]);int nr,n;mi>>nr>>n;assert(nr==n);array<vector<URow>,3>A={vector<URow>(n),vector<URow>(n),vector<URow>(n)};fmpq_t q;fmpq_init(q);
 for(int i=0;i<n;i++){int cnt;mi>>cnt;for(int k=0;k<cnt;k++){int j,terms;mi>>j>>terms;for(int h=0;h<terms;h++){string s;int e,d;mi>>s>>e>>d;if(d>2)continue;assert(fmpq_set_str(q,s.c_str(),10)==0);fmpq_canonicalise(q);fmpq_poly_set_coeff_fmpq(A[d][i][j].n,e,q);}}}
 URat D;uread(di,D.n);int use=bound+1,check=use+2;fmpq*xs=_fmpq_vec_init(check),*ys=_fmpq_vec_init(check*n),*target=_fmpq_vec_init(use);fmpq_t dv,value;for(auto t:{dv,value})fmpq_init(t);fmpz_t den;fmpz_init(den);auto start=chrono::steady_clock::now();
 auto interpolate=[&](int a){ifstream si(prefix+".k"+to_string(a));int nn,N,anchor;si>>nn>>N>>anchor;assert(nn==n&&anchor==-1&&N>=check);for(int i=0;i<check;i++){long x;string ds;si>>x>>ds;fmpq_set_si(xs+i,x,1);assert(fmpz_set_str(den,ds.c_str(),10)==0);fmpq_poly_evaluate_fmpq(dv,D.n,xs+i);assert(!fmpq_is_zero(dv));for(int j=0;j<n;j++){string s;si>>s;assert(fmpz_set_str(fmpq_numref(q),s.c_str(),10)==0);fmpz_set(fmpq_denref(q),den);fmpq_canonicalise(q);fmpq_mul(ys+i*n+j,q,dv);}}
  vector<URat>V(n);long terms=0;int degree=-1,bits=0;for(int j=0;j<n;j++){for(int i=0;i<use;i++)fmpq_set(target+i,ys+i*n+j);assert(fmpq_poly_interpolate_fmpq_vec(V[j].n,xs,target,use));for(int i=use;i<check;i++){fmpq_poly_evaluate_fmpq(value,V[j].n,xs+i);assert(fmpq_equal(value,ys+i*n+j));}degree=max(degree,int(fmpq_poly_degree(V[j].n)));for(int k=0;k<fmpq_poly_length(V[j].n);k++){fmpq_poly_get_coeff_fmpq(q,V[j].n,k);if(!fmpq_is_zero(q)){terms++;bits=max(bits,int(fmpz_bits(fmpq_numref(q))));}}}
  out<<"VECTOR "<<a<<" TERMS "<<terms<<" DEGREE "<<degree<<" MAX_NUMERATOR_BITS "<<bits<<'\n';out.flush();cerr<<"INTERPOLATED "<<a<<" SEC "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';return V;};
 int checked=0;
 auto identity=[&](const vector<pair<int,const vector<URat>*>>&terms){URat tmp,sum;for(int i=0;i<n;i++){fmpq_poly_zero(sum.n);for(auto [d,V]:terms)for(auto&[j,a]:A[d][i])if(!(*V)[j].zero()){fmpq_poly_mul(tmp.n,a.n,(*V)[j].n);fmpq_poly_add(sum.n,sum.n,tmp.n);}if(!sum.zero()){cerr<<"IDENTITY_FAIL "<<checked<<" ROW "<<i<<'\n';abort();}}checked++;cerr<<"EXACT_CHAIN_IDENTITY "<<checked<<" ROWS "<<n<<" SEC "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';};
 // Only the initial four longer-chain vectors are retained after each pair.
 vector<vector<URat>>initial;
 for(int j=0;j<3;j++){auto Q0=interpolate(j),Q1=interpolate(j+3);identity({{0,&Q0}});identity({{0,&Q1},{1,&Q0}});initial.push_back(std::move(Q0));}
 auto R0=interpolate(6),R1=interpolate(7),R2=interpolate(8);identity({{0,&R0}});identity({{0,&R1},{1,&R0}});identity({{0,&R2},{1,&R1},{2,&R0}});initial.push_back(std::move(R0));
 // These are the normalized free coordinates f1=(3,11,12,13) inside
 // the old 14-dimensional kernel. Their determinant is nonzero at e=2.
 // The full 14-dimensional independence is checked separately against the
 // previously verified ten normalized simple kernel columns.
 fmpq_mat_t B;fmpq_mat_init(B,n,4);fmpq_set_si(q,2,1);for(int i=0;i<n;i++)for(int j=0;j<4;j++)fmpq_poly_evaluate_fmpq(fmpq_mat_entry(B,i,j),initial[j][i].n,q);fmpq_mat_t C;fmpq_mat_init(C,n,4);slong rank=fmpq_mat_rref(C,B);assert(rank==4);
 ifstream meta(argv[6]);int pp,mr,mc,rr,pt;meta>>pp>>mr>>mc>>rr>>pt;assert(mr==n&&mc==n&&rr==255);auto readlist=[&](){int len;meta>>len;vector<int>a(len);for(int&v:a)meta>>v;return a;};auto pr=readlist(),pc=readlist(),fc=readlist();vector<int>f1={3,11,12,13};fmpq_mat_t block;fmpq_mat_init(block,4,4);for(int i=0;i<4;i++)for(int j=0;j<4;j++)fmpq_set(fmpq_mat_entry(block,i,j),fmpq_mat_entry(B,fc[f1[i]],j));fmpq_mat_det(value,block);assert(!fmpq_is_zero(value));char*detstr=fmpq_get_str(nullptr,10,value);out<<"FOUR_INITIAL_FREE_BLOCK_DETERMINANT_E2 "<<detstr<<'\n';flint_free(detstr);fmpq_mat_clear(block);out<<"FOUR_INITIAL_VECTOR_RANK "<<rank<<'\n';out<<"EXACT_POLYNOMIAL_CHAIN_IDENTITIES "<<checked*n<<'\n';out<<"TOTAL_SECONDS "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';
 // Save only evaluations needed for the independent fourteen-vector rank.
 out<<"INITIAL_E2 "<<n<<" 4\n";for(int i=0;i<n;i++){for(int j=0;j<4;j++){char*s=fmpq_get_str(nullptr,10,fmpq_mat_entry(B,i,j));out<<s<<' ';flint_free(s);}out<<'\n';}
 fmpq_mat_clear(B);fmpq_mat_clear(C);_fmpq_vec_clear(xs,check);_fmpq_vec_clear(ys,check*n);_fmpq_vec_clear(target,use);fmpz_clear(den);for(auto t:{q,dv,value})fmpq_clear(t);flint_cleanup();}
