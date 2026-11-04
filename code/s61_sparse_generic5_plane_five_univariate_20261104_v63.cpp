// FormalQQ[T] rank witness, no sourcepoint/finitefield evaluation.
#include <algorithm>
#include <cassert>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly_q.h>
using namespace std;
struct P{fmpq_poly_t v;P(){fmpq_poly_init(v);}P(const P&a):P(){fmpq_poly_set(v,a.v);}P(P&&a):P(){fmpq_poly_swap(v,a.v);}P&operator=(P a){fmpq_poly_swap(v,a.v);return *this;}~P(){fmpq_poly_clear(v);}};
struct F{fmpz_poly_q_t v;F(){fmpz_poly_q_init(v);}F(const F&a):F(){fmpz_poly_q_set(v,a.v);}F(F&&a):F(){fmpz_poly_q_swap(v,a.v);}F&operator=(F a){fmpz_poly_q_swap(v,a.v);return *this;}~F(){fmpz_poly_q_clear(v);}};
long maxdegree=0,operations=0,divisions=0;
void check(const P&p){if(fmpq_poly_length(p.v)>150000)throw runtime_error("150k polynomial cap");maxdegree=max(maxdegree,fmpq_poly_degree(p.v));}
bool z(const P&p){return fmpq_poly_is_zero(p.v);}bool z(const F&p){return fmpz_poly_q_is_zero(p.v);}
P mul(const P&a,const P&b){P c;fmpq_poly_mul(c.v,a.v,b.v);check(c);return c;}P sub(const P&a,const P&b){P c;fmpq_poly_sub(c.v,a.v,b.v);check(c);return c;}
P divide(const P&a,const P&b){P c,r;assert(!z(b));fmpq_poly_divrem(c.v,r.v,a.v,b.v);assert(z(r));assert(fmpq_poly_equal(mul(c,b).v,a.v));divisions++;check(c);return c;}
void fc(const F&a){if(fmpz_poly_length(fmpz_poly_q_numref(a.v))>150000||fmpz_poly_length(fmpz_poly_q_denref(a.v))>150000)throw runtime_error("150k field cap");}
F cv(const P&p){F f;fmpq_poly_get_numerator(fmpz_poly_q_numref(f.v),p.v);fmpz_poly_set_fmpz(fmpz_poly_q_denref(f.v),fmpq_poly_denref(p.v));fmpz_poly_q_canonicalise(f.v);fc(f);return f;}
F fm(const F&a,const F&b){F c;fmpz_poly_q_mul(c.v,a.v,b.v);operations++;fc(c);return c;}F fs(const F&a,const F&b){F c;fmpz_poly_q_sub(c.v,a.v,b.v);operations++;fc(c);return c;}F fd(const F&a,const F&b){assert(!z(b));F c;fmpz_poly_q_div(c.v,a.v,b.v);operations++;fc(c);assert(fmpz_poly_q_equal(fm(c,b).v,a.v));return c;}
string qs(const fmpq_t q){char*t=fmpq_get_str(nullptr,10,q);string s(t);flint_free(t);return s;}string zs(const fmpz_t q){char*t=fmpz_get_str(nullptr,10,q);string s(t);flint_free(t);return s;}
void pp(ostream&f,const P&p){f<<'[';bool first=true;fmpq_t q;fmpq_init(q);for(long i=0;i<fmpq_poly_length(p.v);i++){fmpq_poly_get_coeff_fmpq(q,p.v,i);if(fmpq_is_zero(q))continue;if(!first)f<<',';first=false;f<<'['<<i<<",\""<<qs(q)<<"\"]";}f<<']';fmpq_clear(q);}
void ip(ostream&f,const fmpz_poly_t p){f<<'[';bool first=true;fmpz_t q;fmpz_init(q);for(long i=0;i<fmpz_poly_length(p);i++){fmpz_poly_get_coeff_fmpz(q,p,i);if(fmpz_is_zero(q))continue;if(!first)f<<',';first=false;f<<'['<<i<<",\""<<zs(q)<<"\"]";}f<<']';fmpz_clear(q);}
void pf(ostream&f,const F&p){f<<"{\"numerator\":";ip(f,fmpz_poly_q_numref(p.v));f<<",\"denominator\":";ip(f,fmpz_poly_q_denref(p.v));f<<'}';}
void ints(ostream&f,const vector<int>&a){f<<'[';for(size_t i=0;i<a.size();i++){if(i)f<<',';f<<a[i];}f<<']';}
void cap(const string&d){size_t n=0;for(auto&p:filesystem::directory_iterator(d))if(p.is_regular_file())n+=p.file_size();if(n>32*1024*1024)throw runtime_error("32MiB artifact cap");}
int main(int argc,char**argv){assert(argc==3);int status=0;auto start=chrono::steady_clock::now();{
try{ifstream in(argv[1]);string dir=argv[2];int n,m;in>>n>>m;assert(n==134&&m==5);vector<vector<P>>A(n,vector<P>(m));long total=0;fmpq_t q;fmpq_init(q);for(auto&row:A)for(auto&p:row){int count;in>>count;total+=count;for(int k=0;k<count;k++){int e;string c;in>>e>>c;assert(fmpq_set_str(q,c.c_str(),10)==0);fmpq_canonicalise(q);fmpq_poly_set_coeff_fmpq(p.v,e,q);}check(p);}assert(total<=150000);string extra;assert(in&&!(in>>extra));fmpq_clear(q);
 vector<vector<F>>basis;vector<int>rows,pc;vector<F>pivots;vector<int>order(n);for(int i=0;i<n;i++)order[i]=i;sort(order.begin(),order.end(),[&](int i,int j){long x=0,y=0;for(auto&p:A[i])x+=fmpq_poly_length(p.v);for(auto&p:A[j])y+=fmpq_poly_length(p.v);return x==y?i<j:x<y;});
 for(int r:order){vector<F>v(5),mult;for(int j=0;j<5;j++)v[j]=cv(A[r][j]);for(int k=0;k<int(basis.size());k++){F c=v[pc[k]];mult.push_back(c);if(!z(c))for(int j=0;j<5;j++)v[j]=fs(v[j],fm(c,basis[k][j]));assert(z(v[pc[k]]));}int p=-1;for(int j=0;j<5;j++)if(!z(v[j])){p=j;break;}if(p<0)continue;F pivot=v[p];for(int j=0;j<5;j++)v[j]=fd(v[j],pivot);int k=basis.size();rows.push_back(r);pc.push_back(p);pivots.push_back(pivot);basis.push_back(v);
  {ofstream f(dir+"/pivot-"+to_string(k)+".json");f<<"{\"original_row\":"<<r<<",\"pivot_column\":"<<p<<",\"pivot\":";pf(f,pivot);f<<",\"previous_basis_multipliers\":[";for(size_t j=0;j<mult.size();j++){if(j)f<<',';pf(f,mult[j]);}f<<"],\"normalized_row5\":[";for(int j=0;j<5;j++){if(j)f<<',';pf(f,v[j]);}f<<"]}\n";f.close();cap(dir);}cerr<<"FORMAL_UNIVARIATE_PIVOT "<<rows.size()<<" ROW "<<r<<endl;if(rows.size()==5)break;}
 assert(rows.size()==5);vector<vector<P>>B;for(int r:rows)B.push_back(A[r]);P prev;fmpq_poly_one(prev.v);int sign=1;for(int k=0;k<4;k++){int i=k;while(i<5&&z(B[i][k]))i++;assert(i<5);if(i!=k){swap(B[i],B[k]);sign=-sign;}P p=B[k][k];for(int i=k+1;i<5;i++){for(int j=k+1;j<5;j++)B[i][j]=divide(sub(mul(p,B[i][j]),mul(B[i][k],B[k][j])),prev);B[i][k]=P();}prev=p;}
 P determinant=B[4][4];if(sign<0)fmpq_poly_neg(determinant.v,determinant.v);assert(!z(determinant));F product;fmpz_poly_q_one(product.v);int cs=1;for(int i=0;i<5;i++)for(int j=i+1;j<5;j++)if(pc[i]>pc[j])cs=-cs;for(auto&p:pivots)product=fm(product,p);if(cs<0)fmpz_poly_q_neg(product.v,product.v);assert(fmpz_poly_q_equal(product.v,cv(determinant).v));
 {ofstream f(dir+"/native-certificate.json");f<<"{\"all_passed\":true,\"rank\":5,\"selected_original134_rows\":";ints(f,rows);f<<",\"pivot_column_order\":";ints(f,pc);f<<",\"selected_original_matrix5x5\":[";for(int i=0;i<5;i++){if(i)f<<',';f<<'[';for(int j=0;j<5;j++){if(j)f<<',';pp(f,A[rows[i]][j]);}f<<']';}f<<"],\"determinant_polynomial\":";pp(f,determinant);f<<",\"determinant_degree\":"<<fmpq_poly_degree(determinant.v)<<",\"all_Bareiss_divisions_remultiplied\":"<<divisions<<",\"field_pivot_product_equals_determinant\":true,\"input_total_terms\":"<<total<<",\"max_polynomial_degree\":"<<maxdegree<<",\"field_operations\":"<<operations<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";f.close();cap(dir);}
}catch(const exception&e){cerr<<"BOUNDED_UNIVARIATE_STOP "<<e.what()<<endl;status=2;}
}flint_cleanup();return status;}
