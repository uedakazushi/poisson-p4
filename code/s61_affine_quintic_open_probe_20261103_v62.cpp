// Small univariate certificate from saved1031 cyclic data; no largeGB/search.
#include <flint/nmod_poly.h>
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
constexpr ulong P=32003;
struct Poly { nmod_poly_t p; Poly(){nmod_poly_init(p,P);} ~Poly(){nmod_poly_clear(p);} };
void output(ostream&out,const char*name,const nmod_poly_t p){out<<name;for(int i=0;i<=nmod_poly_degree(p);i++)out<<' '<<nmod_poly_get_coeff_ui(p,i);out<<'\n';}
int main(int argc,char**argv){assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);int n;in>>n;assert(n==1031);Poly chi,z[29],products[435],r,g,tmp,q1,q2,rem,u,v,check,check2,tau;
for(int i=0;i<=n;i++){ulong c;in>>c;assert(in&&c<P);nmod_poly_set_coeff_ui(chi.p,i,c);}for(auto &p:z)for(int i=0;i<n;i++){ulong c;in>>c;assert(in&&c<P);nmod_poly_set_coeff_ui(p.p,i,c);}
vector<vector<ulong>>F(134,vector<ulong>(435)),L(134,vector<ulong>(29));vector<ulong>C(134);for(auto&row:F)for(auto&c:row)in>>c;for(auto&row:L)for(auto&c:row)in>>c;for(auto&c:C)in>>c;string extra;assert(!(in>>extra));
int k=0;for(int i=0;i<29;i++)for(int j=i;j<29;j++,k++){nmod_poly_mul(tmp.p,z[i].p,z[j].p);nmod_poly_rem(products[k].p,tmp.p,chi.p);}assert(k==435);
int selected=-1;vector<int>degrees;for(int i=0;i<134;i++){nmod_poly_zero(tmp.p);nmod_poly_set_coeff_ui(tmp.p,0,C[i]);for(int j=0;j<435;j++)if(F[i][j]){nmod_poly_scalar_mul_nmod(check.p,products[j].p,F[i][j]);nmod_poly_add(tmp.p,tmp.p,check.p);}for(int j=0;j<29;j++)if(L[i][j]){nmod_poly_scalar_mul_nmod(check.p,z[j].p,L[i][j]);nmod_poly_add(tmp.p,tmp.p,check.p);}nmod_poly_gcd(g.p,chi.p,tmp.p);degrees.push_back(nmod_poly_degree(g.p));if(selected<0&&nmod_poly_degree(g.p)==5&&nmod_poly_degree(tmp.p)==1030){selected=i;nmod_poly_set(r.p,tmp.p);}}
assert(selected>=0);nmod_poly_gcd(g.p,chi.p,r.p);assert(nmod_poly_degree(g.p)==5);nmod_poly_divrem(q1.p,rem.p,chi.p,g.p);assert(nmod_poly_is_zero(rem.p));nmod_poly_divrem(q2.p,rem.p,r.p,g.p);assert(nmod_poly_is_zero(rem.p));ulong resultant=nmod_poly_resultant(q1.p,q2.p);assert(resultant);
nmod_poly_xgcd(tmp.p,u.p,v.p,q1.p,q2.p);assert(nmod_poly_is_one(tmp.p));assert(nmod_poly_degree(u.p)<1025&&nmod_poly_degree(v.p)<1026);nmod_poly_mul(check.p,u.p,chi.p);nmod_poly_mul(check2.p,v.p,r.p);nmod_poly_add(check.p,check.p,check2.p);assert(nmod_poly_equal(check.p,g.p));
nmod_poly_rem(tau.p,z[0].p,g.p);assert(nmod_poly_is_one(tau.p));nmod_poly_derivative(tmp.p,g.p);ulong discriminant=nmod_poly_resultant(g.p,tmp.p);assert(discriminant);out<<"SELECTED_RESIDUAL "<<selected<<" F_DEGREE 1031 R_DEGREE 1030 G_DEGREE 5\nHIGH_BLOCK_SIZE 2051 COFACTOR_RESULTANT "<<resultant<<" G_DISCRIMINANT "<<discriminant<<" T_NORM 1\n";out<<"INDIVIDUAL_GCD_DEGREES";for(auto d:degrees)out<<' '<<d;out<<'\n';output(out,"RESIDUAL",r.p);output(out,"GCD",g.p);output(out,"BEZOUT_CHI",u.p);output(out,"BEZOUT_RESIDUAL",v.p);for(int i=0;i<29;i++){nmod_poly_rem(tmp.p,z[i].p,g.p);out<<"COORDINATE_MOD_G "<<i;for(int j=0;j<5;j++)out<<' '<<nmod_poly_get_coeff_ui(tmp.p,j);out<<'\n';}cout<<"MONIC_QUINTIC_HIGH_BLOCK_AND_BEZOUT_PASS "<<selected<<' '<<resultant<<'\n';
}
