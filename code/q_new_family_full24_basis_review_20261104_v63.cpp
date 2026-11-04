// Independent formal coefficient engine, (t,s)/LEX, no author includes.
#include <algorithm>
#include <cassert>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <flint/nmod_mpoly.h>
using namespace std;
nmod_mpoly_ctx_t ctx;
constexpr ulong prime=32003;
struct P {
 nmod_mpoly_t v;P(){nmod_mpoly_init(v,ctx);}
 P(const P&a):P(){nmod_mpoly_set(v,a.v,ctx);}P(P&&a):P(){nmod_mpoly_swap(v,a.v,ctx);}
 P&operator=(P a){nmod_mpoly_swap(v,a.v,ctx);return *this;}~P(){nmod_mpoly_clear(v,ctx);}
};
P cnst(ulong z){P a;nmod_mpoly_set_ui(a.v,z,ctx);return a;}
P add(const P&a,const P&b){P c;nmod_mpoly_add(c.v,a.v,b.v,ctx);return c;}
P scale(const P&a,long z){P b;z%=long(prime);if(z<0)z+=prime;nmod_mpoly_scalar_mul_ui(b.v,a.v,z,ctx);return b;}
P mul(const P&a,const P&b){P c;nmod_mpoly_mul(c.v,a.v,b.v,ctx);assert(nmod_mpoly_length(c.v,ctx)<=325);return c;}
P powr(P a,int n){P z=cnst(1);while(n){if(n&1)z=mul(z,a);n>>=1;if(n)a=mul(a,a);}return z;}
using M=vector<P>;
M mm(const M&a,const M&b,int n){M c(n*n);for(int i=0;i<n;i++)for(int j=0;j<n;j++)for(int k=0;k<n;k++)c[n*i+j]=add(c[n*i+j],mul(a[n*i+k],b[n*k+j]));return c;}
P trace(const M&a,int n){P z;for(int i=0;i<n;i++)z=add(z,a[n*i+i]);return z;}
P pairtr(const M&a,const M&b,int n){P z;for(int i=0;i<n;i++)for(int j=0;j<n;j++)z=add(z,mul(a[n*i+j],b[n*j+i]));return z;}
ulong inv(ulong x){ulong z=1;for(ulong n=prime-2;n;n>>=1,x=x*x%prime)if(n&1)z=z*x%prime;return z;}
vector<string> words(int n){vector<string> r;for(int k=0;k<(1<<n);k++){string w;for(int j=n-1;j>=0;j--)w+=char('0'+((k>>j)&1));bool low=true;for(int j=1;j<n;j++)if(w.substr(j)+w.substr(0,j)<w)low=false;if(low)r.push_back(w);}return r;}
map<string,M> prefixes(const M&a,const M&b,int n){map<string,M> r;M eye(n*n);for(int i=0;i<n;i++)eye[n*i+i]=cnst(1);r[""]=eye;for(int d=1;d<=4;d++)for(int k=0;k<(1<<d);k++){string w;for(int j=d-1;j>=0;j--)w+=char('0'+((k>>j)&1));r[w]=mm(r.at(w.substr(0,d-1)),w.back()=='0'?a:b,n);}return r;}
vector<P> traces(const map<string,M>&r,int length,int n){vector<P> ans;for(auto&w:words(length))ans.push_back(pairtr(r.at(w.substr(0,length/2)),r.at(w.substr(length/2)),n));return ans;}
M readmat(ifstream&f,int n){M a(n*n);ulong exps[6][2]={{0,0},{0,1},{1,0},{0,2},{1,1},{2,0}};for(auto&p:a)for(auto&e:exps){ulong z;f>>z;assert(f&&z<prime);if(z)nmod_mpoly_set_coeff_ui_ui(p.v,z,e,ctx);}return a;}
void output(ostream&o,const vector<P>&a,int d){o<<'[';for(int k=0;k<(int)a.size();k++){if(k)o<<',';o<<'[';bool comma=false;for(ulong i=0;i<=ulong(d);i++)for(ulong j=0;j<=ulong(d)-i;j++){if(comma)o<<',';comma=true;ulong e[2]={j,i};o<<nmod_mpoly_get_coeff_ui_ui(a[k].v,e,ctx);}o<<']';}o<<']';}
int main(int argc,char**argv){assert(argc==3);nmod_mpoly_ctx_init(ctx,2,ORD_LEX,prime);{
 ifstream f(argv[1]);M C=readmat(f,4),D=readmat(f,4),T=readmat(f,16),H=readmat(f,16);string extra;assert(!(f>>extra));
 auto cd=prefixes(C,D,4),th=prefixes(T,H,16);auto cd6=traces(cd,6,4),cd8=traces(cd,8,4),th6=traces(th,6,16),th8=traces(th,8,16);
 P qC=trace(cd.at("00"),4),qD=trace(cd.at("11"),4),qCD=trace(cd.at("01"),4);
 P F=add(add(qC,scale(qCD,-4)),scale(qD,4));M cross(16);for(int i=0;i<16;i++)cross[i]=add(cd.at("01")[i],cd.at("10")[i]);
 vector<M> mats={cd.at("00"),cross,cd.at("11")};vector<P> tr={trace(mats[0],4),trace(mats[1],4),trace(mats[2],4)},pencil;
 for(int d=0;d<5;d++){P aa,bb;for(int i=0;i<3;i++){int j=d-i;if(j>=0&&j<3){aa=add(aa,mul(tr[i],tr[j]));bb=add(bb,pairtr(mats[i],mats[j],4));}}pencil.push_back(scale(add(aa,scale(bb,-2)),inv(8)));}
 P N=add(add(pairtr(cd.at("011"),cd.at("011"),4),scale(mul(qD,pairtr(cd.at("00"),cd.at("11"),4)),-long(inv(2)))),mul(pencil[4],qC));
 vector<P>E={mul(F,F),pencil[1],pencil[2],pencil[3],pencil[4]},d12;
 for(auto&e:E)d12.push_back(mul(F,e));d12.push_back(N);for(auto&p:cd6)d12.push_back(p);for(auto&p:th6)d12.push_back(p);
 vector<int>b12={0,1,2,3,4,5,9,10,21};vector<P>basis12,primitive12;for(int i:b12)basis12.push_back(d12[i]);for(int i=5;i<9;i++)primitive12.push_back(basis12[i]);
 vector<P>d16;for(auto&p:basis12)d16.push_back(mul(F,p));for(int i=0;i<5;i++)for(int j=i;j<5;j++)d16.push_back(mul(E[i],E[j]));for(auto&p:cd8)d16.push_back(p);for(auto&p:th8)d16.push_back(p);
 P primitive16=th8[1];vector<P>d24;
 for(int count=0;count<=3;count++){vector<int>indices;auto rec=[&](auto&&self,int lo,int left)->void{if(!left){P p=powr(F,6-2*count);for(int j:indices)p=mul(p,E[j+1]);d24.push_back(p);return;}for(int j=lo;j<4;j++){indices.push_back(j);self(self,j,left-1);indices.pop_back();}};rec(rec,0,count);}
 vector<P>s12={powr(F,3)};for(int j=1;j<5;j++)s12.push_back(mul(F,E[j]));for(auto&p:primitive12)for(auto&q:s12)d24.push_back(mul(p,q));for(int i=0;i<4;i++)for(int j=i;j<4;j++)d24.push_back(mul(primitive12[i],primitive12[j]));for(auto&p:E)d24.push_back(mul(primitive16,p));assert(d24.size()==70);
 ofstream out(argv[2]);out<<"{\"all_passed\":true,\"degree12_columns\":";output(out,d12,12);out<<",\"degree16_columns\":";output(out,d16,16);out<<",\"degree24_columns\":";output(out,d24,24);out<<"}\n";
 }nmod_mpoly_ctx_clear(ctx);flint_cleanup();}
