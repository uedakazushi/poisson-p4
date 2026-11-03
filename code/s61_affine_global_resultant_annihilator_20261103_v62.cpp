#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
constexpr ulong P=32003;
int main(int argc,char**argv){
 assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);int g;in>>g;assert(g==1031);
 auto start=chrono::steady_clock::now();nmod_poly_t chi,pt,cp,q,linear,tmp,quotient,remainder,bezchi,bezpt,unit;
 for(auto x:{chi,pt,cp,q,linear,tmp,quotient,remainder,bezchi,bezpt,unit})nmod_poly_init(x,P);
 for(int i=0;i<=g;i++){ulong c;in>>c;assert(in&&c<P);nmod_poly_set_coeff_ui(chi,i,c);}assert(nmod_poly_degree(chi)==g&&nmod_poly_get_coeff_ui(chi,g)==1);
 for(int i=0;i<g;i++){ulong c;in>>c;assert(in&&c<P);nmod_poly_set_coeff_ui(pt,i,c);}string extra;assert(!(in>>extra));
 nmod_mat_t M;nmod_mat_init(M,g,g,P);vector<ulong>column(g);for(int i=0;i<g;i++)column[i]=nmod_poly_get_coeff_ui(pt,i);
 for(int j=0;j<g;j++){
  for(int i=0;i<g;i++)nmod_mat_entry(M,i,j)=column[i];ulong top=column[g-1];
  for(int i=g-1;i>=0;i--){ulong shifted=(i?column[i-1]:0);column[i]=(shifted+P-top*nmod_poly_get_coeff_ui(chi,i)%P)%P;}
 }
 nmod_mat_charpoly(cp,M);assert(nmod_poly_degree(cp)==g);
 // q(lambda)=Res_s(chi(s),lambda*pt(s)-1)
 //          =(-1)^g lambda^g charpoly_M(1/lambda).
 for(int i=0;i<=g;i++){ulong c=nmod_poly_get_coeff_ui(cp,g-i);if(g%2)c=(P-c)%P;nmod_poly_set_coeff_ui(q,i,c);}
 assert(nmod_poly_get_coeff_ui(q,0)==P-1);
 for(int lambda=0;lambda<=g;lambda++){
  nmod_poly_scalar_mul_nmod(tmp,pt,lambda);nmod_poly_set_coeff_ui(tmp,0,(nmod_poly_get_coeff_ui(tmp,0)+P-1)%P);
  ulong result=nmod_poly_resultant(chi,tmp);assert(result==nmod_poly_evaluate_nmod(q,lambda));
 }
 out<<"G "<<g<<" Q_DEGREE "<<nmod_poly_degree(q)<<" DISTINCT_RESULTANT_CHECKS "<<g+1<<"\nQ";for(int i=0;i<=g;i++)out<<' '<<nmod_poly_get_coeff_ui(q,i);out<<'\n';
 nmod_poly_set_coeff_ui(linear,0,P-1);nmod_poly_set_coeff_ui(linear,1,1);nmod_poly_set(quotient,q);int multiplicity=0;
 while(nmod_poly_degree(quotient)>0){nmod_poly_divrem(tmp,remainder,quotient,linear);if(!nmod_poly_is_zero(remainder))break;nmod_poly_set(quotient,tmp);multiplicity++;}
 out<<"MULTIPLICITY_AT_ONE "<<multiplicity<<'\n';
 for(int lambda:{1,2,3,4,5,7,11})out<<"EVALUATION "<<lambda<<' '<<nmod_poly_evaluate_nmod(q,lambda)<<'\n';
 nmod_poly_scalar_mul_nmod(tmp,pt,2);nmod_poly_set_coeff_ui(tmp,0,(nmod_poly_get_coeff_ui(tmp,0)+P-1)%P);
 nmod_poly_xgcd(unit,bezchi,bezpt,chi,tmp);assert(nmod_poly_degree(unit)==0&&nmod_poly_get_coeff_ui(unit,0)==1);
 out<<"BEZOUT_CHI";for(int i=0;i<=nmod_poly_degree(bezchi);i++)out<<' '<<nmod_poly_get_coeff_ui(bezchi,i);out<<'\n';
 out<<"BEZOUT_2PT_MINUS_ONE";for(int i=0;i<=nmod_poly_degree(bezpt);i++)out<<' '<<nmod_poly_get_coeff_ui(bezpt,i);out<<'\n';
 cerr<<"FULL_RESULTANT_POLYNOMIAL_AND_ALL1032_CHECKS_PASS "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"s\n";
 out.flush();assert(out.good());
}
