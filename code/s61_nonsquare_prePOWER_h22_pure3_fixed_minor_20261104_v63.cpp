// Exact one fixed3x3 minor, Laurent a only. No alternate minor or factor search.
#include <flint/flint.h>
#include <flint/fmpq_mpoly.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <climits>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>
using namespace std;fmpq_mpoly_ctx_t ctx;struct L;set<const L*>alive;long peak=0,reserve_terms=1120;chrono::steady_clock::time_point began;
struct L{fmpq_mpoly_t p;long a=0;L(){fmpq_mpoly_init(p,ctx);alive.insert(this);}L(const L&b):L(){fmpq_mpoly_set(p,b.p,ctx);a=b.a;}L(L&&b):L(){fmpq_mpoly_swap(p,b.p,ctx);a=b.a;b.a=0;}L&operator=(L b){fmpq_mpoly_swap(p,b.p,ctx);swap(a,b.a);return*this;}~L(){alive.erase(this);fmpq_mpoly_clear(p,ctx);}};
void guard(){long n=reserve_terms;for(auto*t:alive)n+=fmpq_mpoly_length(t->p,ctx);peak=max(peak,n);assert(n<=150000);assert(chrono::duration<double>(chrono::steady_clock::now()-began).count()<50);}
bool zero(const L&b){return fmpq_mpoly_is_zero(b.p,ctx);}void norm(L&b){if(zero(b)){b.a=0;return;}ulong ex[3];long m=LONG_MAX;for(long i=0;i<fmpq_mpoly_length(b.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(ex,b.p,i,ctx);m=min(m,long(ex[0]));}if(m){L r;fmpq_t q;fmpq_init(q);for(long i=0;i<fmpq_mpoly_length(b.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(ex,b.p,i,ctx);ex[0]-=m;fmpq_mpoly_get_term_coeff_fmpq(q,b.p,i,ctx);fmpq_mpoly_set_coeff_fmpq_ui(r.p,q,ex,ctx);}b.a+=m;fmpq_mpoly_swap(b.p,r.p,ctx);fmpq_clear(q);}guard();}
L shifted(const L&b,long k){assert(k>=0);L r;ulong ex[3];fmpq_t q;fmpq_init(q);for(long i=0;i<fmpq_mpoly_length(b.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(ex,b.p,i,ctx);ex[0]+=k;fmpq_mpoly_get_term_coeff_fmpq(q,b.p,i,ctx);fmpq_mpoly_set_coeff_fmpq_ui(r.p,q,ex,ctx);}fmpq_clear(q);guard();return r;}
L one(){L r;fmpq_mpoly_one(r.p,ctx);return r;}
L add(const L&b,const L&c,int sign=1){if(zero(b)){L r=c;if(sign<0)fmpq_mpoly_neg(r.p,r.p,ctx);return r;}if(zero(c))return b;long s=min(b.a,c.a);L bb=shifted(b,b.a-s),cc=shifted(c,c.a-s),r;if(sign>0)fmpq_mpoly_add(r.p,bb.p,cc.p,ctx);else fmpq_mpoly_sub(r.p,bb.p,cc.p,ctx);r.a=s;norm(r);guard();return r;}
L mul(const L&b,const L&c){L r;fmpq_mpoly_mul(r.p,b.p,c.p,ctx);r.a=b.a+c.a;norm(r);guard();return r;}
bool eq(const L&b,const L&c){return b.a==c.a&&fmpq_mpoly_equal(b.p,c.p,ctx);}
L readpoly(istream&f){L r;int n;f>>n;assert(n>=0&&n<=150000);fmpq_t q;fmpq_init(q);for(int i=0;i<n;i++){ulong e[3];string c;f>>e[0]>>e[1]>>e[2]>>c;assert(!fmpq_set_str(q,c.c_str(),10));fmpq_canonicalise(q);fmpq_mpoly_set_coeff_fmpq_ui(r.p,q,e,ctx);}fmpq_clear(q);norm(r);return r;}
void jp(ostream&o,const L&b){o<<'[';ulong ex[3];fmpq_t q;fmpq_init(q);for(long i=0;i<fmpq_mpoly_length(b.p,ctx);i++){if(i)o<<',';fmpq_mpoly_get_term_exp_ui(ex,b.p,i,ctx);fmpq_mpoly_get_term_coeff_fmpq(q,b.p,i,ctx);char*s=fmpq_get_str(nullptr,10,q);o<<"[["<<long(ex[0])+b.a<<','<<ex[1]<<','<<ex[2]<<"],\""<<s<<"\"]";flint_free(s);}fmpq_clear(q);o<<']';}
int main(int argc,char**argv){assert(argc==3);began=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);{ifstream f(argv[1]);long n;f>>reserve_terms>>n;assert(n==108);vector<L>v;v.reserve(n);vector<int>oldids;for(int i=0;i<n;i++){int id;char op;f>>id>>op;oldids.push_back(id);L x;if(op=='z'){}else if(op=='1')x=one();else if(op=='p'){long aden;string c;f>>aden>>c;assert(aden>=0);x=readpoly(f);fmpq_t q;fmpq_init(q);assert(!fmpq_set_str(q,c.c_str(),10));fmpq_canonicalise(q);fmpq_mpoly_scalar_div_fmpq(x.p,x.p,q,ctx);fmpq_clear(q);x.a-=aden;norm(x);}else{int a,b;f>>a>>b;assert(a<i&&b<i);if(op=='m')x=mul(v[a],v[b]);else x=add(v[a],v[b],op=='s'?-1:1);}v.push_back(move(x));guard();}
 int pos[3][4];for(auto&row:pos)for(auto&i:row)f>>i;vector<L>Q;for(int i=0;i<4;i++)Q.push_back(readpoly(f));string extra;assert(!(f>>extra));vector<vector<L>>M(3);for(int i=0;i<3;i++)for(int j=0;j<3;j++){M[i].push_back(add(mul(Q[3],v[pos[i][j]]),mul(v[pos[i][3]],Q[j]),-1));assert(M[i][j].a>=0);}
 vector<L>cof;cof.push_back(add(mul(M[1][1],M[2][2]),mul(M[1][2],M[2][1]),-1));cof.push_back(add(mul(M[1][0],M[2][2]),mul(M[1][2],M[2][0]),-1));cof.push_back(add(mul(M[1][0],M[2][1]),mul(M[1][1],M[2][0]),-1));L det=add(add(mul(M[0][0],cof[0]),mul(M[0][1],cof[1]),-1),mul(M[0][2],cof[2]));L alt;vector<int>perm={0,1,2};do{int parity=0;for(int i=0;i<3;i++)for(int j=i+1;j<3;j++)parity+=perm[i]>perm[j];L t=one();for(int i=0;i<3;i++)t=mul(t,M[i][perm[i]]);alt=add(alt,t,parity%2?-1:1);}while(next_permutation(perm.begin(),perm.end()));assert(eq(det,alt));guard();assert(fmpq_mpoly_length(det.p,ctx)<=10080);
 ofstream o(argv[2]);o<<"{\"all_passed\":true,\"fixed_source_actual135\":[18,39,68],\"fixed_columns_original32\":[11,13,25],\"h_column_original32\":22,\"node_values\":[";for(int i=0;i<n;i++){if(i)o<<',';o<<"{\"old_DAG_root\":"<<oldids[i]<<",\"Laurent_polynomial\":";jp(o,v[i]);o<<'}';}o<<"],\"nine_actual_entries\":[";for(int i=0;i<3;i++){if(i)o<<',';o<<'[';for(int j=0;j<3;j++){if(j)o<<',';jp(o,M[i][j]);}o<<']';}o<<"],\"three_signed_cofactor_terms\":[";for(int i=0;i<3;i++){if(i)o<<',';jp(o,cof[i]);}o<<"],\"exact_fixed3_determinant\":";jp(o,det);o<<",\"cofactor_and_six_permutation_identity_verified\":true,\"determinant_is_zero\":"<<(zero(det)?"true":"false")<<",\"determinant_terms\":"<<fmpq_mpoly_length(det.p,ctx)<<",\"determinant_common_a_shift\":"<<det.a<<",\"peak_live_terms_with_input_reserve\":"<<peak<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";}
 fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
