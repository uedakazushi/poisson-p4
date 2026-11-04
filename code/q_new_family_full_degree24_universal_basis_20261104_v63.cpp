// One formal ambient two-plane, coefficient ranks only; no NEW25 source points.
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <gmp.h>
#include <flint/nmod_mpoly.h>
using namespace std;
constexpr ulong PRIME=32003;
nmod_mpoly_ctx_t ctx;
long multiplications=0,maxterms=0;
struct P {
 nmod_mpoly_t v;P(){nmod_mpoly_init(v,ctx);}P(const P&p):P(){nmod_mpoly_set(v,p.v,ctx);}
 P(P&&p):P(){nmod_mpoly_swap(v,p.v,ctx);}P&operator=(P p){nmod_mpoly_swap(v,p.v,ctx);return *this;}
 ~P(){nmod_mpoly_clear(v,ctx);}
};
void guard(const P&p){long n=nmod_mpoly_length(p.v,ctx);maxterms=max(maxterms,n);if(n>325)throw runtime_error("formal degree24 support exceeded");}
P add(const P&a,const P&b){P c;nmod_mpoly_add(c.v,a.v,b.v,ctx);guard(c);return c;}
P scale(long k,const P&a){P c;k%=long(PRIME);if(k<0)k+=PRIME;nmod_mpoly_scalar_mul_ui(c.v,a.v,k,ctx);guard(c);return c;}
P mul(const P&a,const P&b){P c;nmod_mpoly_mul(c.v,a.v,b.v,ctx);++multiplications;guard(c);return c;}
P cnst(ulong k){P a;nmod_mpoly_set_ui(a.v,k,ctx);return a;}
P power(P a,int n){P b=cnst(1);while(n){if(n&1)b=mul(b,a);n>>=1;if(n)a=mul(a,a);}return b;}
using M=vector<P>;
M mm(const M&a,const M&b){M c(16);for(int i=0;i<4;i++)for(int j=0;j<4;j++)for(int k=0;k<4;k++)c[4*i+j]=add(c[4*i+j],mul(a[4*i+k],b[4*k+j]));return c;}
P tr(const M&a){P p;for(int i=0;i<4;i++)p=add(p,a[5*i]);return p;}
P trpair(const M&a,const M&b){P p;for(int i=0;i<4;i++)for(int j=0;j<4;j++)p=add(p,mul(a[4*i+j],b[4*j+i]));return p;}
vector<string> necklaces(int n){vector<string> ws;for(int k=0;k<(1<<n);k++){string w;for(int j=n-1;j>=0;j--)w+=char('0'+((k>>j)&1));string low=w;for(int j=1;j<n;j++)low=min(low,w.substr(j)+w.substr(0,j));if(w==low)ws.push_back(w);}return ws;}
vector<ulong> coeffs(const P&p,int degree){vector<ulong> v;ulong e[2];for(int i=0;i<=degree;i++)for(int j=0;j<=degree-i;j++){e[0]=i;e[1]=j;v.push_back(nmod_mpoly_get_coeff_ui_ui(p.v,e,ctx));}return v;}
ulong modpow(ulong a,ulong n){ulong b=1;while(n){if(n&1)b=b*a%PRIME;a=a*a%PRIME;n>>=1;}return b;}
struct Rank {
 map<int,vector<ulong>> basis;vector<int> selected,pivots;
 bool push(const vector<ulong>&input,int index){auto r=input;for(const auto&pv:basis){ulong q=r[pv.first];if(q)for(int j=pv.first;j<(int)r.size();j++)r[j]=(r[j]+PRIME-q*pv.second[j]%PRIME)%PRIME;}
  int p=0;while(p<(int)r.size()&&!r[p])p++;if(p==(int)r.size())return false;
  ulong inv=modpow(r[p],PRIME-2);for(auto&x:r)x=x*inv%PRIME;basis[p]=move(r);selected.push_back(index);pivots.push_back(p);return true;}
 int rank()const{return basis.size();}
};
void ints(ostream&o,const vector<int>&v){o<<'[';for(int i=0;i<(int)v.size();i++){if(i)o<<',';o<<v[i];}o<<']';}
void poly(ostream&o,const P&p,int degree){auto v=coeffs(p,degree);o<<'[';for(int i=0;i<(int)v.size();i++){if(i)o<<',';o<<v[i];}o<<']';}
void polys(ostream&o,const vector<P>&v,int degree){o<<'[';for(int i=0;i<(int)v.size();i++){if(i)o<<',';poly(o,v[i],degree);}o<<']';}
int main(int argc,char**argv){assert(argc==3);nmod_mpoly_ctx_init(ctx,2,ORD_DEGLEX,PRIME);int result=0;auto start=chrono::steady_clock::now();{
 try{
  ifstream in(argv[1]);int n,width;ulong p;in>>n>>width>>p;assert(n==16&&width==136&&p==PRIME);
  vector<P>b(16);ulong e[2]={0,0};
  for(int j=0;j<3;j++)for(int i=0;i<16;i++){ulong z;in>>z;assert(z<PRIME);e[0]=j==1;e[1]=j==2;nmod_mpoly_set_coeff_ui_ui(b[i].v,z,e,ctx);}
  vector<pair<int,int>>mon(width);for(auto&ij:mon){in>>ij.first>>ij.second;assert(0<=ij.first&&ij.first<=ij.second&&ij.second<16);}
  vector<vector<long>>cc(32,vector<long>(width));for(auto&r:cc)for(auto&x:r)in>>x;assert(in);string extra;assert(!(in>>extra));
  vector<P>bp;for(auto ij:mon)bp.push_back(mul(b[ij.first],b[ij.second]));M C(16),D(16);ulong inv4=modpow(4,PRIME-2);
  for(int r=0;r<32;r++){P sum;for(int k=0;k<width;k++)if(cc[r][k])sum=add(sum,scale(cc[r][k],bp[k]));(r<16?C[r]:D[r-16])=scale(inv4,sum);}
  map<string,M> words;M eye(16);for(int i=0;i<4;i++)eye[5*i]=cnst(1);words[""]=eye;
  for(int len=1;len<=4;len++)for(int k=0;k<(1<<len);k++){string w;for(int j=len-1;j>=0;j--)w+=char('0'+((k>>j)&1));words[w]=mm(words[w.substr(0,len-1)],w.back()=='0'?C:D);}
  auto traceword=[&](const string&w){int k=w.size()/2;return trpair(words.at(w.substr(0,k)),words.at(w.substr(k)));};
  P qC=tr(words["00"]),qD=tr(words["11"]),qCD=tr(words["01"]);
  assert(nmod_mpoly_is_zero(add(qC,scale(-2,qD)).v,ctx));assert(nmod_mpoly_is_zero(add(qCD,qD).v,ctx));
  P F=add(add(qC,scale(-4,qCD)),scale(4,qD));vector<P>pencil(5);vector<int>perm={0,1,2,3};
  do{int sign=1;for(int i=0;i<4;i++)for(int j=i+1;j<4;j++)if(perm[i]>perm[j])sign=-sign;vector<P>cur={cnst(sign==1?1:PRIME-1)};
   for(int i=0;i<4;i++){vector<P>next(cur.size()+1);for(int k=0;k<(int)cur.size();k++){next[k]=add(next[k],mul(cur[k],C[4*i+perm[i]]));next[k+1]=add(next[k+1],mul(cur[k],D[4*i+perm[i]]));}cur=move(next);}for(int k=0;k<5;k++)pencil[k]=add(pencil[k],cur[k]);
  }while(next_permutation(perm.begin(),perm.end()));
  assert(nmod_mpoly_is_zero(add(add(mul(F,F),scale(-400,pencil[0])),scale(-1600,pencil[4])).v,ctx));
  vector<P>E={mul(F,F),pencil[1],pencil[2],pencil[3],pencil[4]};
  P N=add(add(traceword("011011"),scale(-long(modpow(2,PRIME-2)),mul(qD,traceword("0011")))),mul(pencil[4],qC));
  auto w6=necklaces(6),w8=necklaces(8);assert(w6.size()==14&&w8.size()==36);vector<P>t6,t8;
  for(auto w:w6)t6.push_back(traceword(w));for(auto w:w8)t8.push_back(traceword(w));
  vector<P>deg12;for(const auto&p:E)deg12.push_back(mul(F,p));deg12.push_back(N);for(const auto&p:t6)deg12.push_back(p);
  Rank r12;int rank_lower12=0;for(int i=0;i<(int)deg12.size();i++){r12.push(coeffs(deg12[i],12),i);if(i==4)rank_lower12=r12.rank();}
  bool gate12=rank_lower12==5&&r12.selected.size()>=6&&r12.selected[5]==5&&r12.rank()==9;
  vector<P>basis12,prim12,product16,deg16,all24;Rank r16product,r16,r24;int chosen8=-1;bool gate16product=false,gate16=false,gate24=false;
  if(gate12){for(int i:r12.selected)basis12.push_back(deg12[i]);for(int i=5;i<9;i++)prim12.push_back(basis12[i]);
   for(auto&p:basis12)product16.push_back(mul(F,p));for(int i=0;i<5;i++)for(int j=i;j<5;j++)product16.push_back(mul(E[i],E[j]));assert(product16.size()==24);
   deg16=product16;for(auto&p:t8)deg16.push_back(p);for(int i=0;i<24;i++)r16product.push(coeffs(product16[i],16),i);gate16product=r16product.rank()==19;
   if(gate16product){for(int i=0;i<(int)deg16.size();i++)r16.push(coeffs(deg16[i],16),i);gate16=r16.rank()==20;
    if(gate16){assert(r16.selected.size()==20);chosen8=r16.selected.back()-24;assert(chosen8>=0&&chosen8<36);
     for(int r=0;r<=3;r++){vector<int>indices;auto rec=[&](auto&&self,int next,int left)->void{if(!left){P q=power(F,6-2*r);for(int k:indices)q=mul(q,pencil[k+1]);all24.push_back(q);return;}for(int k=next;k<4;k++){indices.push_back(k);self(self,k,left-1);indices.pop_back();}};rec(rec,0,r);}assert(all24.size()==35);
     vector<P>spectral12={power(F,3)};for(int i=1;i<5;i++)spectral12.push_back(mul(F,pencil[i]));
     for(auto&q:prim12)for(auto&s:spectral12)all24.push_back(mul(q,s));for(int i=0;i<4;i++)for(int j=i;j<4;j++)all24.push_back(mul(prim12[i],prim12[j]));for(auto&q:E)all24.push_back(mul(t8[chosen8],q));assert(all24.size()==70);
     for(int i=0;i<70;i++)r24.push(coeffs(all24[i],24),i);gate24=r24.rank()==67;
    }
   }
  }
  string status=!gate12?"degree12_rank_gate_failed":!gate16product?"degree16_product_rank_gate_failed":!gate16?"degree16_rank_gate_failed":!gate24?"degree24_rank_gate_failed":"complete_Inv24_basis_certified";
  ofstream out(argv[2]);out<<"{\"all_passed\":true,\"status\":\""<<status<<"\",\"prime\":32003,\"formal_variable_order\":[\"s\",\"t\"],\"coefficient_order\":\"i=0..d,j=0..d-i: s^i*t^j\",\"F4\":";poly(out,F,4);out<<",\"E5\":";polys(out,E,8);out<<",\"N12\":";poly(out,N,12);
  out<<",\"trace6_columns\":";polys(out,t6,12);out<<",\"trace8_columns\":";polys(out,t8,16);
  out<<",\"degree12_columns\":";polys(out,deg12,12);out<<",\"degree12_lower_rank\":"<<rank_lower12<<",\"degree12_rank\":"<<r12.rank()<<",\"degree12_selected_columns\":";ints(out,r12.selected);out<<",\"degree12_selected_coefficient_rows\":";ints(out,r12.pivots);
  out<<",\"degree16_product_rank\":"<<r16product.rank()<<",\"degree16_columns\":";polys(out,deg16,16);out<<",\"degree16_rank\":"<<r16.rank()<<",\"degree16_selected_columns\":";ints(out,r16.selected);out<<",\"degree16_selected_coefficient_rows\":";ints(out,r16.pivots);out<<",\"primitive16_trace8_index\":"<<chosen8;
  out<<",\"degree24_columns\":";polys(out,all24,24);out<<",\"degree24_rank\":"<<r24.rank()<<",\"degree24_selected_columns\":";ints(out,r24.selected);out<<",\"degree24_selected_coefficient_rows\":";ints(out,r24.pivots);
  out<<",\"max_intermediate_terms\":"<<maxterms<<",\"polynomial_multiplications\":"<<multiplications<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";out.close();cerr<<status<<" RANKS "<<r12.rank()<<' '<<r16product.rank()<<' '<<r16.rank()<<' '<<r24.rank()<<endl;
 }catch(const exception&e){cerr<<"BOUNDED_STOP "<<e.what()<<endl;result=2;}
 }nmod_mpoly_ctx_clear(ctx);flint_cleanup();return result;
}
