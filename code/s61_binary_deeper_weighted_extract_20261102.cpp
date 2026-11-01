#include "s61_binary_deeper_quartic_field_20261101.hpp"
#include <array>
#include <set>
using Mon=array<int,9>;
using PolyRow=map<Mon,Rat>;
int nw[9]={1,1,1,1,1,2,1,2,1};
int weight(const Mon&m){int a=0;for(int j=0;j<9;j++)a+=nw[j]*m[j];return a;}
int cw[45]={1,2,3,1,2,1,3,2,1,2,3,3,2,4,4,3,4,3,2,0,1,2,1,0,1,2,2,1,3,3,2,3,2,1,1,2,2,1,3,3,2,3,2,1,3};
int rw[135]={2,3,3,2,4,4,3,4,3,2,5,5,4,5,4,3,5,4,3,2,3,4,4,3,5,5,4,5,4,3,6,6,5,6,5,4,6,5,4,3,3,4,4,3,5,5,4,5,4,3,6,6,5,6,5,4,6,5,4,3,2,3,2,4,3,2,5,4,3,2,1,1,2,1,2,1,5,4,4,3,3,2,3,2,2,3,3,2,4,4,3,4,3,2,5,5,4,5,4,3,5,4,3,2,6,6,5,6,5,4,6,5,4,3,2,3,3,2,4,4,3,4,3,2,5,5,4,5,4,3,5,4,3,2,3};
void read11(istream&in,PolyRow&out,int expected,bool keep){long n;in>>n;fmpq_t c;fmpq_init(c);ulong e[2];Mon m;for(long k=0;k<n;k++){string a;in>>a>>e[0]>>e[1];for(auto&v:m)in>>v;assert(in&&weight(m)==expected);if(keep){assert(fmpq_set_str(c,a.c_str(),10)==0);fmpq_canonicalise(c);fmpq_mpoly_push_term_fmpq_ui(out[m].n,c,e,ctx);}}if(keep)for(auto&[m,v]:out){fmpq_mpoly_sort_terms(v.n,ctx);fmpq_mpoly_combine_like_terms(v.n,ctx);}fmpq_clear(c);}
void monomials(vector<Mon>&out,int degree,int pos=0,Mon m={}){if(pos==9){if(!degree)out.push_back(m);return;}for(int j=0;j*nw[pos]<=degree;j++){m[pos]=j;monomials(out,degree-j*nw[pos],pos+1,m);}}
int main(int argc,char**argv){assert(argc==6);ifstream weights(argv[5]);for(auto&v:cw)weights>>v;for(auto&v:rw)weights>>v;for(auto&v:nw)weights>>v;assert(weights);fmpq_mpoly_ctx_init(ctx,2,ORD_LEX);{ifstream in(argv[1]),co(argv[2]);ofstream out(argv[3]),meta(argv[4]);assert(in&&co&&out&&meta);int F,np,n,D;in>>F>>np>>n>>D;assert(F==135&&np==2&&n==9&&D==4);Rat den;readpoly(co,den.n);int nk;co>>nk;vector<Rat>known(nk);for(auto&v:known)readpoly(co,v.n);vector<PolyRow>A(F);long checked=0;
for(int d=2;d<=D;d++){for(int i=0;i<45;i++){PolyRow discarded;read11(in,discarded,cw[i],false);}for(int i=0;i<F;i++){PolyRow row;read11(in,row,rw[i],rw[i]<=4);if(rw[i]>4)continue;Rat scale;fmpq_mpoly_pow_ui(scale.n,den.n,4-d,ctx);for(auto&[m,v]:row){Rat z=mul(v,scale);fmpq_mpoly_add(A[i][m].n,A[i][m].n,z.n,ctx);checked++;}}}cerr<<"ALL RAW WEIGHTS CHECKED; LOW COEFFICIENTS "<<checked<<'\n';
int rows=0;for(auto&row:A)if(!row.empty())rows++;out<<rows<<'\n';meta<<"row weight coefficient_count common_content_remaining\n";
for(int i=0;i<F;i++){auto&row=A[i];if(row.empty())continue;Rat g,allowed,q,tmp;for(auto&[m,v]:row){assert(!v.zero());if(g.zero())g=v;else{assert(fmpq_mpoly_gcd(tmp.n,g.n,v.n,ctx));g=tmp;}}fmpq_mpoly_one(allowed.n,ctx);for(auto&factor:known){while(fmpq_mpoly_divides(q.n,g.n,factor.n,ctx)){g=q;fmpq_mpoly_mul(tmp.n,allowed.n,factor.n,ctx);allowed=tmp;}}fmpq_t lc;fmpq_init(lc);fmpq_mpoly_get_term_coeff_fmpq(lc,row.begin()->second.n,0,ctx);fmpq_mpoly_scalar_mul_fmpq(allowed.n,allowed.n,lc,ctx);fmpq_clear(lc);out<<i<<' '<<rw[i]<<' '<<row.size()<<'\n';for(auto&[m,v]:row){assert(fmpq_mpoly_divides(q.n,v.n,allowed.n,ctx));for(int a:m)out<<a<<' ';out<<'\n';writepoly(out,q.n);}meta<<i<<' '<<rw[i]<<' '<<row.size()<<' ';writepoly(meta,g.n);}
vector<Mon>cols;monomials(cols,4);assert(cols.size()==269);cerr<<"WEIGHTED ROWS "<<rows<<" COLUMNS "<<cols.size()<<'\n';}fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
