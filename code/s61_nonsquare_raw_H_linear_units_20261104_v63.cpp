// Individual known-FULLU-unit Gaussian on15x6 raw H-linear remainder only.
#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;fmpq_mpoly_ctx_t ctx;
struct P{fmpq_mpoly_t p;P(){fmpq_mpoly_init(p,ctx);}P(const P&v):P(){fmpq_mpoly_set(p,v.p,ctx);}P(P&&v):P(){fmpq_mpoly_swap(p,v.p,ctx);}P&operator=(P v){fmpq_mpoly_swap(p,v.p,ctx);return*this;}~P(){fmpq_mpoly_clear(p,ctx);}};
struct U{string c="1";array<long,6> e{};};struct Elim{int row;P c;U divisor;};struct Step{int row,col;U pivot,scale_before;vector<Elim> elim;};
vector<P>known;chrono::steady_clock::time_point start;long maxterms=0;string pending="";
void guard(const P&p){long n=fmpq_mpoly_length(p.p,ctx);maxterms=max(maxterms,n);if(n>100000)throw runtime_error("100000_output_terms");if(chrono::duration<double>(chrono::steady_clock::now()-start).count()>50)throw runtime_error("50s_native_boundary");}
P mul(const P&a,const P&b){P z;long n=fmpq_mpoly_length(a.p,ctx),m=fmpq_mpoly_length(b.p,ctx);if(!n||!m)return z;long d=fmpq_mpoly_total_degree_si(a.p,ctx)+fmpq_mpoly_total_degree_si(b.p,ctx);long double bound=min((long double)n*m,(long double)(d+1)*(d+2)*(d+3)/6);if(bound>150000)throw runtime_error("150000_product_support_bound");fmpq_mpoly_mul(z.p,a.p,b.p,ctx);if(fmpq_mpoly_length(z.p,ctx)>150000)throw runtime_error("150000_intermediate_terms");guard(z);return z;}
U combine(const U&a,const U&b,bool divide=false){U z;fmpq_t x,y;fmpq_init(x);fmpq_init(y);assert(!fmpq_set_str(x,a.c.c_str(),10)&&!fmpq_set_str(y,b.c.c_str(),10));if(divide)fmpq_div(x,x,y);else fmpq_mul(x,x,y);char*t=fmpq_get_str(nullptr,10,x);z.c=t;flint_free(t);for(int k=0;k<6;k++)z.e[k]=a.e[k]+(divide?-b.e[k]:b.e[k]);fmpq_clear(x);fmpq_clear(y);return z;}
P upoly(const U&u){P z,t;fmpq_t c;fmpq_init(c);assert(!fmpq_set_str(c,u.c.c_str(),10));fmpq_mpoly_set_fmpq(z.p,c,ctx);for(int k=0;k<6;k++){assert(u.e[k]>=0);for(int j=0;j<u.e[k];j++)z=mul(z,known[k]);}fmpq_clear(c);return z;}
bool unit(P p,U&u){if(fmpq_mpoly_is_zero(p.p,ctx))return false;P before=p,q;for(int k=0;k<6;k++)while(fmpq_mpoly_divides(q.p,p.p,known[k].p,ctx)){p=q;u.e[k]++;}if(!fmpq_mpoly_is_fmpq(p.p,ctx))return false;fmpq_t c;fmpq_init(c);fmpq_mpoly_get_fmpq(c,p.p,ctx);char*t=fmpq_get_str(nullptr,10,c);u.c=t;flint_free(t);fmpq_clear(c);P expected=upoly(u);assert(fmpq_mpoly_equal(before.p,expected.p,ctx));return true;}
U clean(vector<P>&row){
 U u;int first=-1;for(int j=0;j<(int)row.size();j++)if(!fmpq_mpoly_is_zero(row[j].p,ctx)){first=j;break;}if(first<0)return u;
 vector<P>before=row;fmpq_t c;fmpq_init(c);fmpq_mpoly_get_term_coeff_fmpq(c,row[first].p,0,ctx);char*t=fmpq_get_str(nullptr,10,c);u.c=t;flint_free(t);for(P&p:row)fmpq_mpoly_scalar_div_fmpq(p.p,p.p,c,ctx);fmpq_clear(c);
 vector<P>q(row.size());for(int k=0;k<6;k++)while(true){bool yes=true;for(int j=0;j<(int)row.size();j++)if(!fmpq_mpoly_divides(q[j].p,row[j].p,known[k].p,ctx)){yes=false;break;}if(!yes)break;row=q;u.e[k]++;}
 P h=upoly(u);for(int j=0;j<(int)row.size();j++){P check=mul(h,row[j]);assert(fmpq_mpoly_equal(check.p,before[j].p,ctx));guard(row[j]);}return u;
}
void readp(istream&f,P&p){long n;f>>n;assert(f&&n>=0);fmpq_t c;fmpq_init(c);ulong ex[3];for(long i=0;i<n;i++){string s;f>>s>>ex[0]>>ex[1]>>ex[2];assert(f&&!fmpq_set_str(c,s.c_str(),10));fmpq_canonicalise(c);fmpq_mpoly_push_term_fmpq_ui(p.p,c,ex,ctx);}fmpq_mpoly_sort_terms(p.p,ctx);fmpq_mpoly_combine_like_terms(p.p,ctx);assert(fmpq_mpoly_length(p.p,ctx)==n);fmpq_clear(c);}
void jsonu(ostream&f,const U&u){f<<"{\"constant\":\""<<u.c<<"\",\"a_d_F_B_C_J_powers\":[";for(int k=0;k<6;k++){if(k)f<<',';f<<u.e[k];}f<<"]}";}
void jsonp(ostream&f,const P&p){f<<'[';fmpq_t c;fmpq_init(c);ulong ex[3];for(long i=0;i<fmpq_mpoly_length(p.p,ctx);i++){if(i)f<<',';fmpq_mpoly_get_term_coeff_fmpq(c,p.p,i,ctx);fmpq_mpoly_get_term_exp_ui(ex,p.p,i,ctx);char*t=fmpq_get_str(nullptr,10,c);f<<"[["<<ex[0]<<','<<ex[1]<<','<<ex[2]<<"],\""<<t<<"\"]";flint_free(t);if(f.tellp()>28*1024*1024)throw runtime_error("28MiB_checkpoint_guard");}f<<']';fmpq_clear(c);}
int main(int argc,char**argv){
 assert(argc==3);start=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);string path=argv[2],stop="";long completed=0;
 {
  ifstream f(argv[1]);int n,m,k;f>>n>>m>>k;assert(n==15&&m==6&&k==6);for(int j=0;j<6;j++){known.emplace_back();readp(f,known.back());}vector<vector<P>>A(n,vector<P>(m));for(auto&row:A)for(P&p:row)readp(f,p);string extra;assert(!(f>>extra));
  vector<int>rows(n),cols(m);for(int i=0;i<n;i++)rows[i]=i;for(int j=0;j<m;j++)cols[j]=j;vector<U>initial(n),scales(n);vector<Step>steps;U determinant;
  auto save=[&](){
   ofstream o(path+".tmp");o<<"{\"completed_pivots\":"<<steps.size()<<",\"initial_row_divisors\":[";for(int i=0;i<n;i++){if(i)o<<',';jsonu(o,initial[i]);}o<<"],\"steps\":[";
   for(int s=0;s<(int)steps.size();s++){if(s)o<<',';const Step&t=steps[s];o<<"{\"pivot_row\":"<<t.row<<",\"pivot_column\":"<<t.col<<",\"pivot_unit\":";jsonu(o,t.pivot);o<<",\"row_scale_before\":";jsonu(o,t.scale_before);o<<",\"eliminations\":[";for(int v=0;v<(int)t.elim.size();v++){if(v)o<<',';const auto&e=t.elim[v];o<<"{\"row\":"<<e.row<<",\"coefficient_before\":";jsonp(o,e.c);o<<",\"row_divisor_after\":";jsonu(o,e.divisor);o<<'}';}o<<"]}";}
   o<<"],\"original_selected_minor_unit\":";jsonu(o,determinant);o<<",\"remaining_rows\":[";for(int i=0;i<(int)rows.size();i++){if(i)o<<',';o<<rows[i];}o<<"],\"remaining_columns\":[";for(int j=0;j<(int)cols.size();j++){if(j)o<<',';o<<cols[j];}o<<"],\"matrix15x6\":[";for(int i=0;i<n;i++){if(i)o<<',';o<<'[';for(int j=0;j<m;j++){if(j)o<<',';jsonp(o,A[i][j]);}o<<']';}o<<"]}\n";o.close();rename((path+".tmp").c_str(),path.c_str());
  };
  try{
   for(int i=0;i<n;i++){initial[i]=clean(A[i]);scales[i]=combine(U{},initial[i],true);}save();
   while(!rows.empty()&&!cols.empty()){
    int ri=-1,cj=-1;long cost=1000000000;U pivot;
    for(int i:rows)for(int j:cols){long z=fmpq_mpoly_length(A[i][j].p,ctx);if(!z||z>=cost)continue;U u;if(unit(A[i][j],u)){ri=i;cj=j;cost=z;pivot=u;}}
    if(ri<0){stop="no_individual_known_FULLU_unit_entry";break;}
    pending="pivot_row"+to_string(ri)+"_col"+to_string(cj);Step st{ri,cj,pivot,scales[ri],{}};auto B=A;auto nextscale=scales;P p=A[ri][cj];
    for(int i:rows)if(i!=ri&&!fmpq_mpoly_is_zero(A[i][cj].p,ctx)){
     Elim e;e.row=i;e.c=A[i][cj];for(int j=0;j<m;j++){P x=mul(p,A[i][j]),y=mul(e.c,A[ri][j]);fmpq_mpoly_sub(B[i][j].p,x.p,y.p,ctx);guard(B[i][j]);}assert(fmpq_mpoly_is_zero(B[i][cj].p,ctx));e.divisor=clean(B[i]);nextscale[i]=combine(combine(scales[i],pivot),e.divisor,true);st.elim.push_back(move(e));
    }
    long terms=0;for(const auto&r:B)for(const P&q:r)terms+=fmpq_mpoly_length(q.p,ctx);if(terms>150000)throw runtime_error("150000_total_matrix_terms");
    determinant=combine(determinant,combine(pivot,scales[ri],true));A=move(B);scales=move(nextscale);rows.erase(find(rows.begin(),rows.end(),ri));cols.erase(find(cols.begin(),cols.end(),cj));steps.push_back(move(st));save();completed=steps.size();cerr<<"COMPLETED_UNIT "<<completed<<" ROW "<<ri<<" COL "<<cj<<" TERMS "<<terms<<'\n';
   }
  }catch(const exception&e){stop=e.what();remove((path+".tmp").c_str());}
  ofstream a(path+".audit.json");a<<"{\"completed_pivots\":"<<completed<<",\"full_remaining_H6_span\":"<<(completed==6?"true":"false")<<",\"stop_reason\":\""<<stop<<"\",\"unfinished_operation\":\""<<pending<<"\",\"max_polynomial_terms\":"<<maxterms<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
 }
 known.clear();fmpq_mpoly_ctx_clear(ctx);flint_cleanup();return 0;
}
