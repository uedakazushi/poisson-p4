// Five FIXED determinants, then polynomial Bezout over the existingFp only.
#define main frozen_q5_row_main_NOT_CALLED
#include "s61_nonsquare_P100_q5_pure98_row_20261104_v63.cpp"
#undef main
#include <flint/nmod_poly.h>
#include <flint/ulong_extras.h>
#include <sstream>
const ulong PRIME=32003;bool all_pden=true;long pden_checks=0;
P readpoly(istream&f){long n;f>>n;P p;fmpq_t q;fmpq_init(q);for(long i=0;i<n;i++){long e;string c;f>>e>>c;assert(!fmpq_set_str(q,c.c_str(),10));fmpq_canonicalise(q);fmpq_poly_set_coeff_fmpq(p.v,e,q);}fmpq_clear(q);guard(p);return p;}
bool same(const P&a,const P&b){return fmpq_poly_equal(a.v,b.v);}bool zz(const P&p){return fmpq_poly_is_zero(p.v);}
P subp(const P&a,const P&b){return plusp(a,negp(b));}
void pden(const P&p){fmpq_t q;fmpq_init(q);for(long i=0;i<fmpq_poly_length(p.v);i++){fmpq_poly_get_coeff_fmpq(q,p.v,i);all_pden&=fmpz_fdiv_ui(fmpq_denref(q),PRIME)!=0;pden_checks++;}fmpq_clear(q);}
struct Step{int k,swaprow;P pivot,previous;int count;};
struct Node{int id,pivot,entry,left,top,previous;long degree,terms;};
struct Norm{P full,Q,U;string c="1";vector<int>e;bool lcunit=true;};
struct Case{int basis;vector<vector<int>>ex;vector<Step>steps;vector<Node>nodes;int finalnode=-1,sign=1,zero_column=-1;Norm result;};
vector<P>F;
struct NP{nmod_poly_t v;NP(){nmod_poly_init(v,PRIME);}NP(const NP&a):NP(){nmod_poly_set(v,a.v);}NP(NP&&a):NP(){nmod_poly_swap(v,a.v);}NP&operator=(NP a){nmod_poly_swap(v,a.v);return*this;}~NP(){nmod_poly_clear(v);}};
NP modpoly(const P&p){NP r;fmpq_t q;fmpq_init(q);for(long i=0;i<fmpq_poly_length(p.v);i++){fmpq_poly_get_coeff_fmpq(q,p.v,i);ulong d=fmpz_fdiv_ui(fmpq_denref(q),PRIME);assert(d);ulong c=fmpz_fdiv_ui(fmpq_numref(q),PRIME)*n_invmod(d,PRIME)%PRIME;nmod_poly_set_coeff_ui(r.v,i,c);}fmpq_clear(q);return r;}
void jnp(ostream&o,const NP&p){o<<'[';bool first=true;for(long i=0;i<nmod_poly_length(p.v);i++){ulong c=nmod_poly_get_coeff_ui(p.v,i);if(c){if(!first)o<<',';first=false;o<<'['<<i<<','<<c<<']';}}o<<']';}
void ints(ostream&o,const vector<int>&v){o<<'[';for(int i=0;i<int(v.size());i++){if(i)o<<',';o<<v[i];}o<<']';}
struct CapBuf:streambuf{string data;streamsize xsputn(const char*s,streamsize n)override{if(data.size()+n>27*1024*1024)throw runtime_error("27MiB_JSON_guard");data.append(s,n);return n;}int overflow(int c)override{if(c!=EOF){char q=c;xsputn(&q,1);}return c;}};
int main(int argc,char**argv){assert(argc==3);began=chrono::steady_clock::now();string folder=argv[2];try{
 ifstream f(argv[1]);ulong p;f>>input_reserve>>p;assert(p==PRIME);P H=readpoly(f);pden(H);assert(fmpq_poly_degree(H.v)==15);fmpq_t lc;fmpq_init(lc);fmpq_poly_get_coeff_fmpq(lc,H.v,15);assert(fmpq_is_one(lc));fmpq_clear(lc);
 P z=shiftp(one(),1);F={z,plusp(z,pc("20")),plusp(z,pc("60")),plusp(scalar(z,"203"),pc("2700")),plusp(z,pc("12"))};Case c;f>>c.basis;assert(c.basis==9);c.ex=vector<vector<int>>(9,vector<int>(5));for(auto&row:c.ex)for(auto&v:row)f>>v;
 vector<vector<P>>A(9,vector<P>(9));vector<vector<int>>ids(9,vector<int>(9));int nextid=81;long clearcount=0,divcount=0;stage="one81_LCD_source_clearings";
 for(int i=0;i<9;i++){P den=one();for(int k=0;k<5;k++)for(int j=0;j<c.ex[i][k];j++)den=mulp(den,F[k]);
  for(int j=0;j<9;j++){P n=readpoly(f),d=readpoly(f);pden(n);pden(d);A[i][j]=mulp(n,quotient(den,d));assert(same(mulp(A[i][j],d),mulp(n,den)));pden(A[i][j]);ids[i][j]=9*i+j;clearcount++;}}
 string extra;assert(!(f>>extra));P prev=one();int prevnode=-1;bool zero=false;stage="ONE_fixed_Bareiss9";
 for(int k=0;k<8;k++){int r=k;while(r<9&&zz(A[r][k]))r++;if(r==9){zero=true;c.zero_column=k;break;}if(r!=k){swap(A[r],A[k]);swap(ids[r],ids[k]);c.sign=-c.sign;}P pivot=A[k][k];Step step{k,r,pivot,prev,0};int pivotnode=ids[k][k];
  for(int i=k+1;i<9;i++)for(int j=k+1;j<9;j++){P num=subp(mulp(pivot,A[i][j]),mulp(A[i][k],A[k][j]));P q=quotient(num,prev);assert(same(mulp(q,prev),num));pden(q);c.nodes.push_back({nextid,pivotnode,ids[i][j],ids[i][k],ids[k][j],prevnode,fmpq_poly_degree(q.v),fmpq_poly_length(q.v)});ids[i][j]=nextid++;A[i][j]=q;step.count++;divcount++;}
  for(int i=k+1;i<9;i++){A[i][k]=P();ids[i][k]=-2;}prev=pivot;prevnode=pivotnode;c.steps.push_back(move(step));}
 P N=zero?P():(c.sign==1?A[8][8]:negp(A[8][8]));c.finalnode=zero?-2:ids[8][8];pden(N);stage="literal_monicH_remainder";P Q,R;fmpq_poly_divrem(Q.v,R.v,N.v,H.v);guard(Q);guard(R);pden(Q);pden(R);assert(same(N,plusp(mulp(H,Q),R))&&fmpq_poly_degree(R.v)<15);
 NP bh,br,G,s,t;bool modrun=all_pden,unit=false;if(modrun){stage="ENTIRE_finite_special_fiber_unit_gate";bh=modpoly(H);br=modpoly(R);assert(nmod_poly_degree(bh.v)==15&&nmod_poly_get_coeff_ui(bh.v,15)==1);nmod_poly_xgcd(G.v,s.v,t.v,bh.v,br.v);NP a,b,ch;nmod_poly_mul(a.v,s.v,bh.v);nmod_poly_mul(b.v,t.v,br.v);nmod_poly_add(ch.v,a.v,b.v);assert(nmod_poly_equal(ch.v,G.v));unit=nmod_poly_is_one(G.v);}
 CapBuf buf;ostream o(&buf);o.exceptions(ios::badbit|ios::failbit);o<<"{\"all_passed\":true,\"one_fixed_case_basis9\":true,\"all81_source_LCD_identities_verified\":"<<(clearcount==81?"true":"false")<<",\"exact_determinant_N\":";jp(o,N);o<<",\"monic_H15\":";jp(o,H);o<<",\"division_quotient\":";jp(o,Q);o<<",\"finite_remainder_R\":";jp(o,R);o<<",\"literal_N_equals_HQ_plus_R_verified\":true,\"all_QQ_denominators_prime_units\":"<<(all_pden?"true":"false")<<",\"QQ_denominator_checks\":"<<pden_checks<<",\"NO_leading_coefficient_condition_on_N_or_R_used\":true,\"row_swap_sign\":"<<c.sign<<",\"zero_column_frontier\":"<<c.zero_column<<",\"final_Bareiss_node\":"<<c.finalnode<<",\"Bareiss_steps\":[";
 for(int j=0;j<int(c.steps.size());j++){if(j)o<<',';auto&v=c.steps[j];o<<"{\"k\":"<<v.k<<",\"swap_with_row\":"<<v.swaprow<<",\"pivot\":";jp(o,v.pivot);o<<",\"previous_pivot\":";jp(o,v.previous);o<<",\"literal_exact_divisions\":"<<v.count<<'}';}
 o<<"],\"exact_quotient_DAG\":[";for(int j=0;j<int(c.nodes.size());j++){if(j)o<<',';auto&v=c.nodes[j];o<<'['<<v.id<<','<<v.pivot<<','<<v.entry<<','<<v.left<<','<<v.top<<','<<v.previous<<','<<v.degree<<','<<v.terms<<']';}
 o<<"],\"prime\":32003,\"finite_modular_gate_run\":"<<(modrun?"true":"false")<<",\"H_reduction\":";jnp(o,bh);o<<",\"R_reduction\":";jnp(o,br);o<<",\"finite_modular_gcd\":";jnp(o,G);o<<",\"modular_Bezout_H_cofactor\":";jnp(o,s);o<<",\"modular_Bezout_R_cofactor\":";jnp(o,t);o<<",\"literal_modular_Bezout_verified\":"<<(modrun?"true":"false")<<",\"entire_special_fiber_quotient_dimension\":"<<(modrun?nmod_poly_degree(G.v):-1)<<",\"finite_special_fiber_unit\":"<<(unit?"true":"false")<<",\"whole_P100_quartic_lower150\":"<<(unit?"true":"false")<<",\"Bareiss_exact_divisions_total\":"<<divcount<<",\"peak_live_polynomial_terms_with_literal_reserve\":"<<peak_live+input_reserve<<",\"literal_input_reserve\":"<<input_reserve<<",\"serialized_terms\":"<<stored<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";
 ofstream of(folder+"/native.json");of<<buf.data;
 }catch(const exception&e){ofstream o(folder+"/frontier.json");o<<"{\"all_passed\":false,\"stage\":\""<<stage<<"\",\"reason\":\""<<e.what()<<"\"}\n";return 2;}flint_cleanup();}
