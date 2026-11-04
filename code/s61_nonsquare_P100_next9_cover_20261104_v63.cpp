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
Norm normalize(const P&d){Norm r;r.full=d;r.Q=d;r.U=one();r.e=vector<int>(5);pden(d);if(zz(d))return r;
 for(int k=0;k<5;k++)for(int n=0;n<1000;n++){P q,s;fmpq_poly_divrem(q.v,s.v,r.Q.v,F[k].v);guard(q);guard(s);if(!zz(s))break;assert(same(mulp(q,F[k]),r.Q));pden(q);r.Q=q;r.U=mulp(r.U,F[k]);r.e[k]++;assert(n<999);}
 fmpz_t common,content,z,tmp;fmpz_init(common);fmpz_one(common);fmpz_init(content);fmpz_init(z);fmpz_init(tmp);fmpq_t c,q;fmpq_init(c);fmpq_init(q);
 for(long i=0;i<fmpq_poly_length(r.Q.v);i++){fmpq_poly_get_coeff_fmpq(q,r.Q.v,i);fmpz_lcm(common,common,fmpq_denref(q));}
 for(long i=0;i<fmpq_poly_length(r.Q.v);i++){fmpq_poly_get_coeff_fmpq(q,r.Q.v,i);fmpz_divexact(tmp,common,fmpq_denref(q));fmpz_mul(z,tmp,fmpq_numref(q));fmpz_gcd(content,content,z);}
 fmpq_poly_get_coeff_fmpq(q,r.Q.v,fmpq_poly_degree(r.Q.v));if(fmpq_sgn(q)<0)fmpz_neg(content,content);fmpq_set_fmpz_frac(c,content,common);fmpq_canonicalise(c);char*s=fmpq_get_str(nullptr,10,c);r.c=s;flint_free(s);P prim;fmpq_poly_scalar_div_fmpq(prim.v,r.Q.v,c);guard(prim);r.Q=prim;pden(r.Q);
 fmpz_zero(content);for(long i=0;i<fmpq_poly_length(r.Q.v);i++){fmpq_poly_get_coeff_fmpq(q,r.Q.v,i);assert(fmpz_is_one(fmpq_denref(q)));fmpz_gcd(content,content,fmpq_numref(q));}assert(fmpz_is_one(content));
 fmpq_poly_get_coeff_fmpq(q,r.Q.v,fmpq_poly_degree(r.Q.v));assert(fmpq_sgn(q)>0);r.lcunit=fmpz_fdiv_ui(fmpq_numref(q),PRIME)!=0;all_pden&=fmpz_fdiv_ui(fmpq_denref(c),PRIME)!=0;assert(same(scalar(mulp(r.U,r.Q),r.c),d));
 fmpq_clear(c);fmpq_clear(q);fmpz_clear(common);fmpz_clear(content);fmpz_clear(z);fmpz_clear(tmp);return r;
}
struct NP{nmod_poly_t v;NP(){nmod_poly_init(v,PRIME);}NP(const NP&a):NP(){nmod_poly_set(v,a.v);}NP(NP&&a):NP(){nmod_poly_swap(v,a.v);}NP&operator=(NP a){nmod_poly_swap(v,a.v);return*this;}~NP(){nmod_poly_clear(v);}};
NP modpoly(const P&p){NP r;fmpq_t q;fmpq_init(q);for(long i=0;i<fmpq_poly_length(p.v);i++){fmpq_poly_get_coeff_fmpq(q,p.v,i);ulong d=fmpz_fdiv_ui(fmpq_denref(q),PRIME);assert(d);ulong c=fmpz_fdiv_ui(fmpq_numref(q),PRIME)*n_invmod(d,PRIME)%PRIME;nmod_poly_set_coeff_ui(r.v,i,c);}fmpq_clear(q);return r;}
void jnp(ostream&o,const NP&p){o<<'[';bool first=true;for(long i=0;i<nmod_poly_length(p.v);i++){ulong c=nmod_poly_get_coeff_ui(p.v,i);if(c){if(!first)o<<',';first=false;o<<'['<<i<<','<<c<<']';}}o<<']';}
void ints(ostream&o,const vector<int>&v){o<<'[';for(int i=0;i<int(v.size());i++){if(i)o<<',';o<<v[i];}o<<']';}
struct CapBuf:streambuf{string data;streamsize xsputn(const char*s,streamsize n)override{if(data.size()+n>27*1024*1024)throw runtime_error("27MiB_JSON_guard");data.append(s,n);return n;}int overflow(int c)override{if(c!=EOF){char q=c;xsputn(&q,1);}return c;}};
int main(int argc,char**argv){assert(argc==3);began=chrono::steady_clock::now();string folder=argv[2];try{
 ifstream f(argv[1]);ulong p;int count;f>>input_reserve>>p>>count;assert(p==PRIME&&count==5);P z=shiftp(one(),1);F={z,plusp(z,pc("20")),plusp(z,pc("60")),plusp(scalar(z,"203"),pc("2700")),plusp(z,pc("12"))};vector<Case>cases;long divcount=0,clearcount=0;
 for(int ca=0;ca<5;ca++){
  Case c;f>>c.basis;c.ex=vector<vector<int>>(9,vector<int>(5));for(auto&row:c.ex)for(auto&x:row)f>>x;vector<vector<P>>A(9,vector<P>(9));vector<vector<int>>ids(9,vector<int>(9));int nextid=81;
  stage="source_matrix_"+to_string(ca);for(int i=0;i<9;i++){P den=one();for(int k=0;k<5;k++)for(int j=0;j<c.ex[i][k];j++)den=mulp(den,F[k]);
   for(int j=0;j<9;j++){P n=readpoly(f),d=readpoly(f);pden(n);pden(d);A[i][j]=mulp(n,quotient(den,d));assert(same(mulp(A[i][j],d),mulp(n,den)));pden(A[i][j]);ids[i][j]=9*i+j;clearcount++;}}
  P prev=one();int prevnode=-1;bool iszero=false;stage="fixed_Bareiss9_"+to_string(ca);
  for(int k=0;k<8;k++){
   int r=k;while(r<9&&zz(A[r][k]))r++;if(r==9){iszero=true;c.zero_column=k;break;}
   if(r!=k){swap(A[r],A[k]);swap(ids[r],ids[k]);c.sign=-c.sign;}P pivot=A[k][k];Step step{k,r,pivot,prev,0};int pivotnode=ids[k][k];
   for(int i=k+1;i<9;i++)for(int j=k+1;j<9;j++){
    P num=subp(mulp(pivot,A[i][j]),mulp(A[i][k],A[k][j]));P q=quotient(num,prev);assert(same(mulp(q,prev),num));pden(q);
    c.nodes.push_back({nextid,pivotnode,ids[i][j],ids[i][k],ids[k][j],prevnode,fmpq_poly_degree(q.v),fmpq_poly_length(q.v)});ids[i][j]=nextid++;A[i][j]=q;step.count++;divcount++;
   }
   for(int i=k+1;i<9;i++){A[i][k]=P();ids[i][k]=-2;}prev=pivot;prevnode=pivotnode;c.steps.push_back(move(step));
  }
  P det=iszero?P():(c.sign==1?A[8][8]:negp(A[8][8]));c.finalnode=iszero?-2:ids[8][8];c.result=normalize(det);cases.push_back(move(c));
 }
 string extra;assert(!(f>>extra));bool lc=true,nonzero=false;for(auto&c:cases){lc&=c.result.lcunit;nonzero|=!zz(c.result.Q);}bool modrun=all_pden&&lc&&nonzero;vector<NP>bars(5),weights(5);NP G;vector<int>degrees;bool unit=false;
 if(modrun){stage="coefficient_polynomial_mod32003_Bezout";for(int i=0;i<5;i++){bars[i]=modpoly(cases[i].result.Q);assert(nmod_poly_degree(bars[i].v)==fmpq_poly_degree(cases[i].result.Q.v));}
  for(int i=0;i<5;i++){if(nmod_poly_is_zero(G.v)&&nmod_poly_is_zero(bars[i].v)){degrees.push_back(-1);continue;}NP g,s,t,l,r,check;nmod_poly_xgcd(g.v,s.v,t.v,G.v,bars[i].v);nmod_poly_mul(l.v,s.v,G.v);nmod_poly_mul(r.v,t.v,bars[i].v);nmod_poly_add(check.v,l.v,r.v);assert(nmod_poly_equal(check.v,g.v));
   for(int j=0;j<i;j++){NP q;nmod_poly_mul(q.v,weights[j].v,s.v);weights[j]=move(q);}weights[i]=t;G=g;degrees.push_back(nmod_poly_degree(G.v));}
  NP check;for(int i=0;i<5;i++){NP q;nmod_poly_mul(q.v,weights[i].v,bars[i].v);nmod_poly_add(check.v,check.v,q.v);}assert(nmod_poly_equal(check.v,G.v));unit=nmod_poly_is_one(G.v);
 }
 CapBuf buf;ostream o(&buf);o.exceptions(ios::badbit|ios::failbit);o<<"{\"all_passed\":true,\"prime\":32003,\"all405_source_LCD_identities_verified\":"<<(clearcount==405?"true":"false")<<",\"all_QQ_denominators_prime_units\":"<<(all_pden?"true":"false")<<",\"QQ_coefficient_den_checks\":"<<pden_checks<<",\"all_nonzero_primitive_LC_prime_units\":"<<(lc?"true":"false")<<",\"determinants5\":[";
 for(int ci=0;ci<5;ci++){if(ci)o<<',';auto&c=cases[ci];o<<"{\"cover_basis_index\":"<<c.basis<<",\"exact_determinant\":";jp(o,c.result.full);o<<",\"known_unit_multiplicities5\":";ints(o,c.result.e);o<<",\"known_unit_polynomial\":";jp(o,c.result.U);o<<",\"QQ_unit_content\":\""<<c.result.c<<"\",\"primitive_ZZ_polynomial\":";jp(o,c.result.Q);o<<",\"primitive_LC_prime_unit\":"<<(c.result.lcunit?"true":"false")<<",\"row_swap_sign\":"<<c.sign<<",\"zero_column_frontier\":"<<c.zero_column<<",\"final_Bareiss_node\":"<<c.finalnode<<",\"Bareiss_steps\":[";
  for(int j=0;j<int(c.steps.size());j++){if(j)o<<',';auto&s=c.steps[j];o<<"{\"k\":"<<s.k<<",\"swap_with_row\":"<<s.swaprow<<",\"pivot\":";jp(o,s.pivot);o<<",\"previous_pivot\":";jp(o,s.previous);o<<",\"literal_exact_divisions\":"<<s.count<<'}';}
  o<<"],\"exact_quotient_DAG\":[";for(int j=0;j<int(c.nodes.size());j++){if(j)o<<',';auto&n=c.nodes[j];o<<'['<<n.id<<','<<n.pivot<<','<<n.entry<<','<<n.left<<','<<n.top<<','<<n.previous<<','<<n.degree<<','<<n.terms<<']';}o<<"]}";
 }
 o<<"],\"modular_gcd_run\":"<<(modrun?"true":"false")<<",\"all_reduced_degrees_preserved\":"<<(modrun?"true":"false")<<",\"modular_gcd_is1\":"<<(unit?"true":"false")<<",\"modular_gcd\":";jnp(o,G);o<<",\"modular_gcd_degree_steps\":";ints(o,degrees);o<<",\"modular_Bezout5\":[";
 for(int i=0;i<5;i++){if(i)o<<',';jnp(o,weights[i]);}o<<"],\"all5_modular_primitive_polynomials\":[";for(int i=0;i<5;i++){if(i)o<<',';jnp(o,bars[i]);}o<<"],\"literal_modular_Bezout_verified\":"<<(modrun?"true":"false")<<",\"whole_curve_quartic_lower150\":"<<(unit&&all_pden&&lc?"true":"false")<<",\"Bareiss_exact_divisions_total\":"<<divcount<<",\"peak_live_polynomial_terms_with_literal_reserve\":"<<peak_live+input_reserve<<",\"literal_input_reserve\":"<<input_reserve<<",\"serialized_terms\":"<<stored<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";
 ofstream of(folder+"/native.json");of<<buf.data;
 }catch(const exception&e){ofstream o(folder+"/frontier.json");o<<"{\"all_passed\":false,\"stage\":\""<<stage<<"\",\"reason\":\""<<e.what()<<"\"}\n";return 2;}flint_cleanup();}
