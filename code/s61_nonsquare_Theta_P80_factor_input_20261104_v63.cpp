// P80closedfactor source17x24 INPUT reductions and certifiedchartunits ONLY.
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
void readp(istream&f,P&p){int n;f>>n;assert(n>=0);fmpq_t c;fmpq_init(c);for(int i=0;i<n;i++){int d;string v;f>>d>>v;assert(d>=0&&!fmpq_set_str(c,v.c_str(),10));fmpq_canonicalise(c);fmpq_poly_set_coeff_fmpq(p.p,d,c);}fmpq_clear(c);}
void jp(ostream&o,const P&p){o<<'[';fmpq_t c;fmpq_init(c);bool first=true;for(long i=fmpq_poly_degree(p.p);i>=0;i--){fmpq_poly_get_coeff_fmpq(c,p.p,i);if(fmpq_is_zero(c))continue;if(!first)o<<',';first=false;char*s=fmpq_get_str(nullptr,10,c);o<<'['<<i<<",\""<<s<<"\"]";flint_free(s);}o<<']';fmpq_clear(c);assert(o.tellp()<31*1024*1024);}
int main(int argc,char**argv){assert(argc==3);start=chrono::steady_clock::now();ifstream f(argv[1]);int nf;f>>nf;assert(nf>0&&nf<=46);vector<P>F(nf);for(auto&p:F)readp(f,p);P U;readp(f,U);int nr,nc;f>>nr>>nc;assert(nr==17&&nc==24);vector<P>M(nr*nc);for(auto&p:M)readp(f,p);string extra;assert(!(f>>extra));ofstream o(argv[2]);o<<"{\"components\":[";
 for(int k=0;k<nf;k++){P g,s,t,ui;fmpq_poly_xgcd(g.p,s.p,t.p,F[k].p,U.p);assert(fmpq_poly_is_one(g.p));assert(eq(add(mul(s,F[k]),mul(t,U)),one()));fmpq_poly_rem(ui.p,t.p,F[k].p);P chk;fmpq_poly_mul(chk.p,U.p,ui.p);fmpq_poly_rem(chk.p,chk.p,F[k].p);assert(fmpq_poly_is_one(chk.p));vector<P>qs(nr*nc),rs(nr*nc);long terms=0,bits=0,nonzero=0;fmpq_t c;fmpq_init(c);for(int i=0;i<nr*nc;i++){fmpq_poly_divrem(qs[i].p,rs[i].p,M[i].p,F[k].p);assert(eq(add(mul(qs[i],F[k]),rs[i]),M[i]));assert(fmpq_poly_degree(rs[i].p)<fmpq_poly_degree(F[k].p));nonzero+=!fmpq_poly_is_zero(rs[i].p);for(int j=0;j<=fmpq_poly_degree(rs[i].p);j++){fmpq_poly_get_coeff_fmpq(c,rs[i].p,j);if(fmpq_is_zero(c))continue;terms++;bits=max(bits,(long)(fmpz_bits(fmpq_numref(c))+fmpz_bits(fmpq_denref(c))));}}fmpq_clear(c);assert(terms<=150000);if(k)o<<',';o<<"{\"factor_index\":"<<k<<",\"field_degree\":"<<fmpq_poly_degree(F[k].p)<<",\"factor\":";jp(o,F[k]);o<<",\"chart_unit_inverse\":";jp(o,ui);o<<",\"chart_unit_Bezout1\":[";jp(o,s);o<<',';jp(o,t);o<<"],\"matrix_quotients17x24\":[";for(int i=0;i<nr;i++){if(i)o<<',';o<<'[';for(int j=0;j<nc;j++){if(j)o<<',';jp(o,qs[i*nc+j]);}o<<']';}o<<"],\"matrix_remainders17x24\":[";for(int i=0;i<nr;i++){if(i)o<<',';o<<'[';for(int j=0;j<nc;j++){if(j)o<<',';jp(o,rs[i*nc+j]);}o<<']';}o<<"],\"reduced_rational_coefficient_count\":"<<terms<<",\"nonzero_matrix_entries\":"<<nonzero<<",\"max_numerator_plus_denominator_bits\":"<<bits<<'}';}
 o<<"],\"rank_or_rightinverse_jobs\":0,\"all_source_reductions_and_chart_unit_Bezouts_checked\":true,\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";}
