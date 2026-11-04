// Exactunivariateexception products, derivative-xgcd, squarefree and root-supportcertificate.
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;chrono::steady_clock::time_point start;
void tick(){assert(chrono::duration<double>(chrono::steady_clock::now()-start).count()<50);}
struct P{fmpq_poly_t p;P(){fmpq_poly_init(p);}P(const P&x):P(){fmpq_poly_set(p,x.p);}P(P&&x):P(){fmpq_poly_swap(p,x.p);}P&operator=(P x){fmpq_poly_swap(p,x.p);return*this;}~P(){fmpq_poly_clear(p);}};
P one(){P p;fmpq_poly_one(p.p);return p;}P add(const P&a,const P&b){P p;fmpq_poly_add(p.p,a.p,b.p);return p;}P mul(const P&a,const P&b){P p;fmpq_poly_mul(p.p,a.p,b.p);assert(fmpq_poly_length(p.p)<150000);tick();return p;}bool eq(const P&a,const P&b){return fmpq_poly_equal(a.p,b.p);}P divexact(const P&a,const P&b){P q,r;fmpq_poly_divrem(q.p,r.p,a.p,b.p);assert(fmpq_poly_is_zero(r.p)&&eq(mul(q,b),a));return q;}
void readp(istream&f,P&p){int n;f>>n;assert(n>0);fmpq_t c;fmpq_init(c);for(int i=0;i<n;i++){int d;string v;f>>d>>v;assert(d>=0&&!fmpq_set_str(c,v.c_str(),10));fmpq_canonicalise(c);fmpq_poly_set_coeff_fmpq(p.p,d,c);}fmpq_clear(c);}
void jp(ostream&o,const P&p){o<<'[';fmpq_t c;fmpq_init(c);bool first=true;for(long i=fmpq_poly_degree(p.p);i>=0;i--){fmpq_poly_get_coeff_fmpq(c,p.p,i);if(fmpq_is_zero(c))continue;if(!first)o<<',';first=false;char*s=fmpq_get_str(nullptr,10,c);o<<'['<<i<<",\""<<s<<"\"]";flint_free(s);}o<<']';fmpq_clear(c);assert(o.tellp()<31*1024*1024);}
int main(int argc,char**argv){assert(argc==3);start=chrono::steady_clock::now();ifstream f(argv[1]);int groups;f>>groups;assert(groups==2);ofstream o(argv[2]);o<<"{\"groups\":[";
 for(int g=0;g<2;g++){int n;f>>n;assert(n>0);vector<P>input(n),monic(n);P prod=one();vector<string>lcs;for(int i=0;i<n;i++){readp(f,input[i]);fmpq_t c;fmpq_init(c);fmpq_poly_get_coeff_fmpq(c,input[i].p,fmpq_poly_degree(input[i].p));char*s=fmpq_get_str(nullptr,10,c);lcs.push_back(s);flint_free(s);fmpq_clear(c);fmpq_poly_make_monic(monic[i].p,input[i].p);prod=mul(prod,monic[i]);}
  P derivative,G,s,t,S,dquot,Sprime,sg,ss,st;fmpq_poly_derivative(derivative.p,prod.p);fmpq_poly_xgcd(G.p,s.p,t.p,prod.p,derivative.p);assert(eq(add(mul(s,prod),mul(t,derivative)),G));S=divexact(prod,G);dquot=divexact(derivative,G);fmpq_poly_derivative(Sprime.p,S.p);fmpq_poly_xgcd(sg.p,ss.p,st.p,S.p,Sprime.p);assert(fmpq_poly_is_one(sg.p));assert(eq(add(mul(ss,S),mul(st,Sprime)),one()));
  vector<P>removed,after,Squot;P remaining=prod;while(fmpq_poly_degree(remaining.p)>0){P cg;fmpq_poly_gcd(cg.p,remaining.p,S.p);assert(fmpq_poly_degree(cg.p)>0);removed.push_back(cg);remaining=divexact(remaining,cg);after.push_back(remaining);Squot.push_back(divexact(S,cg));assert(removed.size()<=400);}
  assert(fmpq_poly_is_one(remaining.p));if(g)o<<',';o<<"{\"normalized_input_factors\":[";for(int i=0;i<n;i++){if(i)o<<',';jp(o,monic[i]);}o<<"],\"original_input_leading_constants\":[";for(int i=0;i<n;i++){if(i)o<<',';o<<'\"'<<lcs[i]<<'\"';}o<<"],\"monic_product_P\":";jp(o,prod);o<<",\"derivative_Pprime\":";jp(o,derivative);o<<",\"derivative_gcd_G\":";jp(o,G);o<<",\"Pprime_quotient_by_G\":";jp(o,dquot);o<<",\"P_Pprime_Bezout\":[";jp(o,s);o<<',';jp(o,t);o<<"],\"squarefree_S_P_over_G\":";jp(o,S);o<<",\"S_derivative\":";jp(o,Sprime);o<<",\"S_Sprime_Bezout1\":[";jp(o,ss);o<<',';jp(o,st);o<<"],\"root_support_removal_steps\":[";for(int i=0;i<(int)removed.size();i++){if(i)o<<',';o<<"{\"removed_factor\":";jp(o,removed[i]);o<<",\"remaining_after\":";jp(o,after[i]);o<<",\"S_quotient\":";jp(o,Squot[i]);o<<'}';}o<<"],\"squarefree_degree\":"<<fmpq_poly_degree(S.p)<<'}';
 }
 string extra;assert(!(f>>extra));o<<"],\"all_products_derivative_xgcd_root_support_identities_checked\":true,\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";flint_cleanup();}
