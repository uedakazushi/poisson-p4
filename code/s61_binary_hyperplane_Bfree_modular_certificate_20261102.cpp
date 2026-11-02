#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/nmod_mat.h>
#include <flint/ulong_extras.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <cassert>
#include <chrono>
#include <algorithm>
using namespace std;
ulong prime;
struct Poly{
 nmod_poly_t p;Poly(){nmod_poly_init(p,prime);}Poly(const Poly&a):Poly(){nmod_poly_set(p,a.p);}
 Poly(Poly&&a):Poly(){nmod_poly_swap(p,a.p);}Poly&operator=(Poly a){nmod_poly_swap(p,a.p);return *this;}~Poly(){nmod_poly_clear(p);}
};
using Row=vector<Poly>;
Row readrow(istream&in,int n){Row r(n);int count;in>>count;assert(count>=0&&count<=n);vector<bool>used(n);for(int z=0;z<count;z++){int j,len;in>>j>>len;assert(j>=0&&j<n&&!used[j]&&len>=0);used[j]=true;for(int k=0;k<len;k++){ulong c;in>>c;assert(c<prime);nmod_poly_set_coeff_ui(r[j].p,k,c);}}assert(in);return r;}
pair<int,int>leading(const Row&r){int degree=-1,pos=-1;for(int j=0;j<(int)r.size();j++){int d=nmod_poly_degree(r[j].p);if(d>=0&&(d>degree||(d==degree&&j>pos))){degree=d;pos=j;}}return {degree,pos};}
int main(int argc,char**argv){
 assert(argc==4);ifstream in(argv[1]),selection(argv[2]);ofstream out(argv[3]);assert(in&&selection&&out);
 int nr,n,nt;in>>prime>>nr>>n>>nt;assert(prime==32003&&nr==321&&n==238&&nt==7);ulong b0,b1;in>>b0>>b1;
 vector<Row>A;for(int i=0;i<nr;i++)A.push_back(readrow(in,n));for(int i=0;i<nt;i++)readrow(in,n);string extra;assert(!(in>>extra));
 int ns;selection>>ns;assert(ns==n);vector<int>rows(n);vector<bool>used(nr);for(int&i:rows){selection>>i;assert(i>=0&&i<nr&&!used[i]);used[i]=true;}assert(selection);assert(!(selection>>extra));
 auto start=chrono::steady_clock::now();auto sec=[&](){return chrono::duration<double>(chrono::steady_clock::now()-start).count();};
 nmod_poly_mat_t M;nmod_poly_mat_init(M,n,n,prime);for(int i=0;i<n;i++)for(int j=0;j<n;j++)nmod_poly_set(nmod_poly_mat_entry(M,i,j),A[rows[i]][j].p);
 Poly det;nmod_poly_mat_det_interpolate(det.p,M);int degree=nmod_poly_degree(det.p);assert(degree>=0);
 out<<"DETERMINANT_DEGREE "<<degree<<" LEADING "<<nmod_poly_get_coeff_ui(det.p,degree)<<" SECONDS "<<sec()<<'\n';out<<"DETERMINANT_COEFFICIENTS "<<degree+1;for(int i=0;i<=degree;i++)out<<' '<<nmod_poly_get_coeff_ui(det.p,i);out<<'\n';out.flush();assert(out.good());
 // Independent evaluations of the returned polynomial at three points.
 nmod_mat_t E;nmod_mat_init(E,n,n,prime);for(ulong point:{0UL,3UL,17UL}){
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)nmod_mat_entry(E,i,j)=nmod_poly_evaluate_nmod(nmod_poly_mat_entry(M,i,j),point);
  ulong value=nmod_mat_det(E);assert(value==nmod_poly_evaluate_nmod(det.p,point));out<<"EVALUATION "<<point<<' '<<value<<'\n';
 }
 nmod_mat_clear(E);nmod_poly_mat_clear(M);double detseconds=sec();
 // Polynomial row operations only; no nonconstant row scaling or division.
 vector<Row>basis(n);vector<int>degrees(n,-1);long steps=0;int rank=0;Poly tmp,shift;
 for(int i=0;i<nr;i++){
  Row r=A[i];while(true){auto[deg,pos]=leading(r);if(pos<0)break;
   if(basis[pos].empty()){basis[pos]=std::move(r);degrees[pos]=deg;rank++;break;}
   if(deg<degrees[pos]){swap(r,basis[pos]);swap(deg,degrees[pos]);}
   ulong denom=nmod_poly_get_coeff_ui(basis[pos][pos].p,degrees[pos]);assert(denom);
   ulong c=nmod_poly_get_coeff_ui(r[pos].p,deg)*n_invmod(denom,prime)%prime;
   for(int j=0;j<n;j++){
    nmod_poly_scalar_mul_nmod(tmp.p,basis[pos][j].p,c);nmod_poly_shift_left(shift.p,tmp.p,deg-degrees[pos]);nmod_poly_sub(r[j].p,r[j].p,shift.p);
   }
   steps++;assert(sec()-detseconds<30);
  }
 }
 assert(rank==n);long degree_sum=0;for(int d:degrees){assert(d==0);degree_sum+=d;}
 nmod_mat_init(E,n,n,prime);out<<"MODULE_DIAGONAL "<<n;
 for(int i=0;i<n;i++){
  assert(leading(basis[i])==make_pair(0,i));
  for(int j=0;j<n;j++){
   assert(nmod_poly_degree(basis[i][j].p)<=0);ulong c=nmod_poly_get_coeff_ui(basis[i][j].p,0);nmod_mat_entry(E,i,j)=c;if(j>i)assert(!c);
  }
  ulong d=nmod_mat_entry(E,i,i);assert(d);out<<' '<<d;
 }
 out<<'\n';ulong basisdet=nmod_mat_det(E);assert(basisdet);nmod_mat_clear(E);
 out<<"MODULE_RANK "<<rank<<" DEGREE "<<degree_sum<<" BASIS_DET "<<basisdet<<" OPERATIONS "<<steps<<" SECONDS "<<sec()-detseconds<<'\n';out<<"TOTAL_SECONDS "<<sec()<<'\n';out.flush();assert(out.good());
 cerr<<"PASS degree="<<degree<<" LC="<<nmod_poly_get_coeff_ui(det.p,degree)<<" module="<<rank<<"/"<<degree_sum<<" sec="<<sec()<<'\n';
}
