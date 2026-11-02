#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod.h>
#include <array>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include <cassert>
using namespace std;
constexpr int eqs=134,n=23,n2=276,n3=2300,m3=3082,p=32003;
struct Z {nmod_mat_t m; Z(int r,int c){nmod_mat_init(m,r,c,p);} ~Z(){nmod_mat_clear(m);} ulong& at(int i,int j){return nmod_mat_entry(m,i,j);}};
static vector<int> selectrows(Z &a,int rank){
    nmod_mat_t E;nmod_mat_init_set(E,a.m);vector<slong> perm(a.m->r);iota(perm.begin(),perm.end(),0);
    auto r=nmod_mat_lu(perm.data(),E,0);assert(r==rank);nmod_mat_clear(E);return {perm.begin(),perm.begin()+rank};
}
int main(int argc,char**argv){
 assert(argc==3);ifstream f(argv[1]);ofstream o(argv[2]);int np;f>>np;
 int pair[n][n],tri[n][n][n],ct=0;
 for(int i=0;i<n;i++)for(int j=i;j<n;j++)pair[i][j]=pair[j][i]=ct++;assert(ct==n2);ct=0;
 for(int i=0;i<n;i++)for(int j=i;j<n;j++)for(int k=j;k<n;k++){array<int,3> a{i,j,k};do{tri[a[0]][a[1]][a[2]]=ct;}while(next_permutation(a.begin(),a.end()));ct++;}assert(ct==n3);
 Z H(eqs,n2),H3(m3,n3);
 for(int r=0;r<eqs;r++)for(int i=0;i<n;i++)for(int j=i;j<n;j++){
  ulong v;f>>v;assert(f&&v<p);H.at(r,pair[i][j])=v;for(int k=0;k<n;k++)H3.at(r*n+k,tri[i][j][k])=v;
 }
 Z HT(n2,eqs);nmod_mat_transpose(HT.m,H.m);auto p2=selectrows(HT,eqs),p3=selectrows(H3,n3);
 vector<int>cf;vector<bool>used2(n2);for(int j:p2)used2[j]=true;
 for(int j=0;j<n2;j++)if(!used2[j])cf.push_back(j);
 // Fixed quadratic and cubic pivots. Basis is constant, 23 linears, 142 quadrics.
 const int m=166;
 Z A2(eqs,eqs),inv2(eqs,eqs),mat3(n3,n3),rhs3(n3,m),red3(n3,m);
 for(int i=0;i<eqs;i++)for(int j=0;j<eqs;j++)A2.at(i,j)=H.at(i,p2[j]);
 assert(nmod_mat_inv(inv2.m,A2.m));
 for(int i=0;i<n3;i++)for(int j=0;j<n3;j++)mat3.at(i,j)=H3.at(p3[i],j);
 vector<int>exponentL(n);for(int &v:exponentL)f>>v;
 o<<"pivots "<<nmod_mat_det(A2.m)<<" "<<nmod_mat_det(mat3.m)<<"\n";
 for(int sample=0;sample<np;sample++){
  int label;f>>label;Z L(eqs,n+1);for(int i=0;i<eqs;i++)for(int j=0;j<=n;j++)f>>L.at(i,j);
  Z rf(eqs,m),reduced(eqs,m),red2(n2,m);
  for(int i=0;i<eqs;i++){rf.at(i,0)=L.at(i,n);for(int k=0;k<n;k++)rf.at(i,1+k)=L.at(i,k);for(int k=0;k<int(cf.size());k++)rf.at(i,1+n+k)=H.at(i,cf[k]);}
  nmod_mat_mul(reduced.m,inv2.m,rf.m);
  for(int i=0;i<eqs;i++)for(int j=0;j<m;j++)red2.at(p2[i],j)=(p-reduced.at(i,j))%p;
  for(int k=0;k<int(cf.size());k++)red2.at(cf[k],1+n+k)=1;
  // Reduce the low parts of the selected X_j F_r rows to the retained symbols.
  for(int i=0;i<n3;i++){int r=p3[i]/n,j=p3[i]%n;
   for(int k=0;k<m;k++){ulong v= (k==j+1?L.at(r,n):0);
    for(int t=0;t<n;t++)v+=L.at(r,t)*red2.at(pair[j][t],k);
    rhs3.at(i,k)=(p-v%p)%p;
   }
  }
  assert(nmod_mat_solve(red3.m,mat3.m,rhs3.m));
  // Each row represents a basis monomial times the same fixed linear form.
  Z T(m,m),K(m,m),last(1,m),row(1,m);
  for(int i=0;i<n;i++){
   T.at(0,i+1)=exponentL[i];
   for(int j=0;j<n;j++)for(int k=0;k<m;k++)T.at(j+1,k)+=exponentL[i]*red2.at(pair[i][j],k);
   int u=0;for(int j=0;j<n;j++)for(int t=j;t<n;t++)if(!used2[pair[j][t]]){
    assert(cf[u]==pair[j][t]);for(int k=0;k<m;k++)T.at(n+1+u,k)+=exponentL[i]*red3.at(tri[i][j][t],k);u++;}
  }
  for(int i=0;i<m;i++)for(int j=0;j<m;j++)T.at(i,j)%=p;
  row.at(0,0)=1;
  for(int i=0;i<m;i++){for(int j=0;j<m;j++)K.at(i,j)=row.at(0,j);
   nmod_mat_mul(last.m,row.m,T.m);nmod_mat_set(row.m,last.m);}
  o<<sample<<' '<<label<<" cyclic "<<nmod_mat_rank(K.m)<<" det "<<nmod_mat_det(K.m)<<"\n";o.flush();
  if(nmod_mat_rank(K.m)==m){
   // Solve simultaneously for the 23 coordinate polynomials and s^166.
   Z KT(m,m),targets(m,n+1),co(m,n+1);nmod_mat_transpose(KT.m,K.m);
   for(int i=0;i<n;i++)targets.at(i+1,i)=1;
   for(int i=0;i<m;i++)targets.at(i,n)=row.at(0,i);
   assert(nmod_mat_solve(co.m,KT.m,targets.m));
   nmod_poly_t chi,cayley,G,val,tmp,h;
   for(auto* poly:{chi,cayley,G,val,tmp,h})nmod_poly_init(poly,p);
   nmod_poly_set_coeff_ui(chi,m,1);
   for(int i=0;i<m;i++)nmod_poly_set_coeff_ui(chi,i,(p-co.at(i,n))%p);
   nmod_mat_charpoly(cayley,T.m);assert(nmod_poly_equal(chi,cayley));
   nmod_poly_struct vars[n],quadrics[n2];
   for(int j=0;j<n;j++){
    nmod_poly_init(vars+j,p);
    for(int i=0;i<m;i++)nmod_poly_set_coeff_ui(vars+j,i,co.at(i,j));
    nmod_poly_scalar_addmul_nmod(h,vars+j,exponentL[j]);}
   assert(nmod_poly_degree(h)==1&&nmod_poly_get_coeff_ui(h,0)==0&&nmod_poly_get_coeff_ui(h,1)==1);
   for(int i=0;i<n;i++)for(int j=i;j<n;j++){int k=pair[i][j];nmod_poly_init(quadrics+k,p);nmod_poly_mulmod(quadrics+k,vars+i,vars+j,chi);}
   // The gcd of the actual original equations, including repeated roots.
   nmod_poly_set(G,chi);
   std::vector<int>steps;
   for(int i=0;i<eqs;i++){nmod_poly_zero(val);nmod_poly_set_coeff_ui(val,0,L.at(i,n));
    for(int j=0;j<n;j++)nmod_poly_scalar_addmul_nmod(val,vars+j,L.at(i,j));
    for(int j=0;j<n2;j++)nmod_poly_scalar_addmul_nmod(val,quadrics+j,H.at(i,j));
    auto d=nmod_poly_degree(G);nmod_poly_gcd(tmp,G,val);nmod_poly_set(G,tmp);
    if(nmod_poly_degree(G)<d)steps.push_back(i);
   }
   o<<"univariate "<<nmod_poly_degree(G)<<" gcd ";
   for(int i=0;i<=nmod_poly_degree(G);i++)o<<nmod_poly_get_coeff_ui(G,i)<<' ';
   o<<" steps ";for(int i:steps)o<<i<<' ';o<<"\n";
   o<<"chi ";for(int i=0;i<=m;i++)o<<nmod_poly_get_coeff_ui(chi,i)<<' ';o<<"\n";
   for(auto* poly:{chi,cayley,G,val,tmp,h})nmod_poly_clear(poly);
   for(int j=0;j<n;j++)nmod_poly_clear(vars+j);for(int j=0;j<n2;j++)nmod_poly_clear(quadrics+j);
  }

 }
}
