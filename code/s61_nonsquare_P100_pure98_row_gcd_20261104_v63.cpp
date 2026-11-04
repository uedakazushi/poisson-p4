#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;auto began=chrono::steady_clock::now();long stored=0,peak=0;
struct P{fmpq_poly_t v;P(){fmpq_poly_init(v);}P(const P&a):P(){fmpq_poly_set(v,a.v);}P(P&&a):P(){fmpq_poly_swap(v,a.v);}P&operator=(P a){fmpq_poly_swap(v,a.v);return*this;}~P(){fmpq_poly_clear(v);}};
void guard(const P&p){peak=max(peak,fmpq_poly_length(p.v));assert(peak<150000);assert(chrono::duration<double>(chrono::steady_clock::now()-began).count()<50);}
P one(){P p;fmpq_poly_one(p.v);return p;}P add(const P&a,const P&b){P p;fmpq_poly_add(p.v,a.v,b.v);guard(p);return p;}P mul(const P&a,const P&b){P p;fmpq_poly_mul(p.v,a.v,b.v);guard(p);return p;}
P quotient(const P&a,const P&b){P q,r;fmpq_poly_divrem(q.v,r.v,a.v,b.v);assert(fmpq_poly_is_zero(r.v));assert(fmpq_poly_equal(mul(q,b).v,a.v));guard(q);return q;}
P readp(istream&f){long n;f>>n;P p;fmpq_t q;fmpq_init(q);for(long i=0;i<n;i++){long e;string s;f>>e>>s;assert(!fmpq_set_str(q,s.c_str(),10));fmpq_canonicalise(q);fmpq_poly_set_coeff_fmpq(p.v,e,q);}fmpq_clear(q);guard(p);return p;}
void jp(ostream&o,const P&p){stored+=fmpq_poly_length(p.v);assert(stored<150000);o<<'[';bool first=true;fmpq_t q;fmpq_init(q);for(long j=0;j<fmpq_poly_length(p.v);j++){fmpq_poly_get_coeff_fmpq(q,p.v,j);if(fmpq_is_zero(q))continue;if(!first)o<<',';first=false;char*s=fmpq_get_str(nullptr,10,q);o<<"[["<<j<<"],\""<<s<<"\"]";flint_free(s);}fmpq_clear(q);o<<']';}
int main(int argc,char**argv){assert(argc==3);began=chrono::steady_clock::now();{
 ifstream f(argv[1]);vector<P>N(98),D(98);for(int i=0;i<98;i++){N[i]=readp(f);D[i]=readp(f);assert(!fmpq_poly_is_zero(N[i].v)&&!fmpq_poly_is_zero(D[i].v));}vector<string>names(5);vector<P>U(5);for(int i=0;i<5;i++){f>>names[i];U[i]=readp(f);assert(fmpq_poly_degree(U[i].v)==1);}string extra;assert(!(f>>extra));
 P G;for(auto&p:N){P g;fmpq_poly_gcd(g.v,G.v,p.v);G=g;}assert(!fmpq_poly_is_zero(G.v));vector<P>Q(98);for(int i=0;i<98;i++)Q[i]=quotient(N[i],G);
 vector<P>A(98);P H;vector<int>tested;for(int i=0;i<98&&!fmpq_poly_equal(H.v,G.v);i++){P h,s,t;fmpq_poly_xgcd(h.v,s.v,t.v,H.v,N[i].v);assert(fmpq_poly_equal(add(mul(s,H),mul(t,N[i])).v,h.v));for(int j=0;j<98;j++)if(!fmpq_poly_is_zero(A[j].v))A[j]=mul(s,A[j]);A[i]=add(A[i],t);H=h;tested.push_back(i);}assert(fmpq_poly_equal(H.v,G.v));P sum;for(int i=0;i<98;i++)sum=add(sum,mul(A[i],N[i]));assert(fmpq_poly_equal(sum.v,G.v));
 struct S{string name;P g,q,u;};vector<S>steps;P remaining=G;for(int k=0;k<5;k++){for(int j=0;;j++){assert(j<100);P g;fmpq_poly_gcd(g.v,remaining.v,U[k].v);if(fmpq_poly_degree(g.v)==0)break;P q=quotient(remaining,g),u=quotient(U[k],g);steps.push_back({names[k],g,q,u});remaining=q;}}
 bool unit=fmpq_poly_degree(remaining.v)==0;P reconstruct=remaining;for(auto&s:steps)reconstruct=mul(reconstruct,s.g);assert(fmpq_poly_equal(reconstruct.v,G.v));
 // All original denominators are independently checked against the same units.
 for(auto&dd:D){P p=dd;for(auto&u:U)for(int j=0;;j++){assert(j<100);P g;fmpq_poly_gcd(g.v,p.v,u.v);if(fmpq_poly_degree(g.v)==0)break;p=quotient(p,g);}assert(fmpq_poly_degree(p.v)==0);}
 ofstream o(argv[2]);o<<"{\"gcd_all98_numerators\":";jp(o,G);o<<",\"all98_quotients\":[";for(int i=0;i<98;i++){if(i)o<<',';jp(o,Q[i]);}o<<"],\"literal_Bezout_contributors\":[";bool first=true;for(int i=0;i<98;i++)if(!fmpq_poly_is_zero(A[i].v)){if(!first)o<<',';first=false;o<<"{\"position\":"<<i<<",\"coefficient\":";jp(o,A[i]);o<<'}';}o<<"],\"literal_sum_Ai_Ni_equals_gcd\":true,\"known_unit_strip_steps\":[";for(int i=0;i<int(steps.size());i++){if(i)o<<',';auto&s=steps[i];o<<"{\"unit\":\""<<s.name<<"\",\"gcd\":";jp(o,s.g);o<<",\"source_quotient\":";jp(o,s.q);o<<",\"unit_quotient\":";jp(o,s.u);o<<'}';}o<<"],\"remaining_unknown\":";jp(o,remaining);o<<",\"row_unimodular_on_whole_P100\":"<<(unit?"true":"false")<<",\"all98_denominators_curve_units_checked\":true,\"global_regular_row_dual\":[";first=true;if(unit)for(int i=0;i<98;i++)if(!fmpq_poly_is_zero(A[i].v)){if(!first)o<<',';first=false;o<<"{\"position\":"<<i<<",\"numerator\":";jp(o,mul(A[i],D[i]));o<<",\"denominator\":";jp(o,G);o<<'}';}o<<"],\"extended_gcd_input_positions_tested\":[";for(int i=0;i<int(tested.size());i++){if(i)o<<',';o<<tested[i];}o<<"],\"stored_terms\":"<<stored<<",\"peak_polynomial_length\":"<<peak<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";
 }flint_cleanup();}
