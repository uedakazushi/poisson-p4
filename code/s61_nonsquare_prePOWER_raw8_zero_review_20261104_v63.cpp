// Independent fixed190 source replay +90 QQ polynomial equalities, DEGREVLEX.
#include <flint/fmpq_mpoly.h>
#include <chrono>
#include <cassert>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <set>
using namespace std;
fmpq_mpoly_ctx_t C;auto began=chrono::steady_clock::now();long peak=0;
struct L;set<const L*>live;
struct L{fmpq_mpoly_t p;long a=0;L(){fmpq_mpoly_init(p,C);live.insert(this);}L(const L&v):L(){fmpq_mpoly_set(p,v.p,C);a=v.a;}~L(){live.erase(this);fmpq_mpoly_clear(p,C);}};
void guard(){long n=15000;for(auto x:live)n+=fmpq_mpoly_length(x->p,C);peak=max(peak,n);assert(n<=150000);assert(chrono::duration<double>(chrono::steady_clock::now()-began).count()<50);}
void polyread(istream&f,L&x){long n;f>>n;fmpq_t q;fmpq_init(q);for(long i=0;i<n;i++){ulong e[3];string s;f>>e[0]>>e[1]>>e[2]>>s;assert(!fmpq_set_str(q,s.c_str(),10));fmpq_canonicalise(q);fmpq_mpoly_set_coeff_fmpq_ui(x.p,q,e,C);}fmpq_clear(q);assert(fmpq_mpoly_is_canonical(x.p,C));guard();}
void shiftcopy(L&out,const L&in,long exp){assert(exp>=0);if(!exp){fmpq_mpoly_set(out.p,in.p,C);return;}L m;ulong e[3]={(ulong)exp,0,0};fmpq_t one;fmpq_init(one);fmpq_one(one);fmpq_mpoly_set_coeff_fmpq_ui(m.p,one,e,C);fmpq_clear(one);fmpq_mpoly_mul(out.p,in.p,m.p,C);guard();}
void multiply(L&out,const L&x,const L&y){out.a=x.a+y.a;fmpq_mpoly_mul(out.p,x.p,y.p,C);guard();}
void combine(L&out,const L&x,const L&y,bool subtract){out.a=min(x.a,y.a);L xx,yy;shiftcopy(xx,x,x.a-out.a);shiftcopy(yy,y,y.a-out.a);if(subtract)fmpq_mpoly_sub(out.p,xx.p,yy.p,C);else fmpq_mpoly_add(out.p,xx.p,yy.p,C);guard();}
bool equal(const L&x,const L&y){L xx,yy;long a=min(x.a,y.a);shiftcopy(xx,x,x.a-a);shiftcopy(yy,y,y.a-a);return fmpq_mpoly_equal(xx.p,yy.p,C);}
int main(int argc,char**argv){assert(argc==3);fmpq_mpoly_ctx_init(C,3,ORD_DEGREVLEX);{
 ifstream in(argv[1]);long n;in>>n;assert(n==190);vector<unique_ptr<L>>v;vector<long>oldids;for(long i=0;i<n;i++){long old;char op;in>>old>>op;oldids.push_back(old);auto x=make_unique<L>();if(op=='p'){long shift;string unit;in>>shift>>unit;polyread(in,*x);x->a=shift;fmpq_t q;fmpq_init(q);assert(!fmpq_set_str(q,unit.c_str(),10));fmpq_canonicalise(q);fmpq_mpoly_scalar_div_fmpq(x->p,x->p,q,C);fmpq_clear(q);}else{long a,b;in>>a>>b;assert(a<i);if(op=='k'){fmpq_mpoly_set(x->p,v[a]->p,C);x->a=v[a]->a+b;}else{assert(b<i);if(op=='m')multiply(*x,*v[a],*v[b]);else{assert(op=='a'||op=='s');combine(*x,*v[a],*v[b],op=='s');}}}v.push_back(move(x));guard();}
 vector<long>roots(24);for(auto&i:roots){in>>i;assert(i>=0&&i<n);}vector<unique_ptr<L>>saved;for(long i=0;i<24;i++){auto p=make_unique<L>();polyread(in,*p);assert(equal(*v[roots[i]],*p));saved.push_back(move(p));}
 string extra;assert(!(in>>extra));long minors=0;vector<vector<long>>labels;for(long i=0;i<4;i++)for(long j=i+1;j<4;j++)for(long a=0;a<6;a++)for(long b=a+1;b<6;b++){L x,y,d;multiply(x,*saved[6*i+a],*saved[6*j+b]);multiply(y,*saved[6*i+b],*saved[6*j+a]);combine(d,x,y,true);assert(fmpq_mpoly_is_zero(d.p,C));labels.push_back({i,j,a,b});minors++;}
 assert(minors==90);ofstream out(argv[2]);out<<"{\"all_passed\":true,\"fixed_source_operations\":190,\"literal_cached24_values_match\":true,\"all90_QQ_two_minors_zero\":true,\"different_monomial_order\":\"DEGREVLEX\",\"peak_visible_terms_including15000reserve\":"<<peak<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<",\"zero_minor_labels\":[";for(long k=0;k<90;k++){if(k)out<<',';out<<'[';for(int j=0;j<4;j++){if(j)out<<',';out<<labels[k][j];}out<<']';}out<<"]}\n";
 }fmpq_mpoly_ctx_clear(C);flint_cleanup();return 0;}
