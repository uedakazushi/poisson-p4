#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <flint/fmpq_mpoly_factor.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <chrono>
#include <algorithm>
using namespace std;fmpq_mpoly_ctx_t ctx;
struct Poly{fmpq_mpoly_t p;Poly(){fmpq_mpoly_init(p,ctx);}Poly(const Poly&a):Poly(){fmpq_mpoly_set(p,a.p,ctx);}Poly(Poly&&a):Poly(){fmpq_mpoly_swap(p,a.p,ctx);}Poly&operator=(Poly a){fmpq_mpoly_swap(p,a.p,ctx);return *this;}~Poly(){fmpq_mpoly_clear(p,ctx);}};
void write(ostream&o,const Poly&a){fmpq_t q;fmpq_init(q);ulong ex[3];int n=fmpq_mpoly_length(a.p,ctx);o<<n;for(int i=0;i<n;i++){fmpq_mpoly_get_term_coeff_fmpq(q,a.p,i,ctx);fmpq_mpoly_get_term_exp_ui(ex,a.p,i,ctx);char*s=fmpq_get_str(nullptr,10,q);o<<' '<<s;for(ulong x:ex)o<<' '<<x;flint_free(s);}o<<'\n';fmpq_clear(q);}


vector<Poly> known;
bool unit(Poly p){if(fmpq_mpoly_is_zero(p.p,ctx))return false;Poly q;for(auto&f:known)while(fmpq_mpoly_divides(q.p,p.p,f.p,ctx))p=q;return fmpq_mpoly_is_fmpq(p.p,ctx);}
void clean(vector<Poly>&row){int first=-1;for(int j=0;j<(int)row.size();j++)if(!fmpq_mpoly_is_zero(row[j].p,ctx)){first=j;break;}if(first<0)return;fmpq_t c;fmpq_init(c);fmpq_mpoly_get_term_coeff_fmpq(c,row[first].p,0,ctx);for(auto&p:row)fmpq_mpoly_scalar_div_fmpq(p.p,p.p,c,ctx);fmpq_clear(c);vector<Poly> quot(row.size());for(auto&f:known)while(true){bool yes=true;for(int j=first;j<(int)row.size();j++)if(!fmpq_mpoly_divides(quot[j].p,row[j].p,f.p,ctx)){yes=false;break;}if(!yes)break;row=quot;}}
int main(int argc,char**argv){assert(argc==3);fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);{const char*names[]={"a","e","d"};for(auto text:{"a","500*a+90*e^2+3*e-24","50*a*d+18*d*e^2-9*d*e-2*d+60*e-40","100*a+18*e^2-9*e-2"}){Poly p;assert(!fmpq_mpoly_set_str_pretty(p.p,text,names,ctx));known.push_back(p);}ifstream in(argv[1]);ofstream out(argv[2]);int n,m,r;in>>n>>m>>r;assert(n==63&&m==91);vector<vector<Poly>>A(n,vector<Poly>(m));fmpq_t q;fmpq_init(q);ulong ex[3];for(auto&row:A)for(auto&p:row){int l;in>>l;for(int k=0;k<l;k++){string t;in>>t>>ex[0]>>ex[1]>>ex[2];assert(!fmpq_set_str(q,t.c_str(),10));fmpq_canonicalise(q);fmpq_mpoly_push_term_fmpq_ui(p.p,q,ex,ctx);}fmpq_mpoly_sort_terms(p.p,ctx);fmpq_mpoly_combine_like_terms(p.p,ctx);}for(auto&row:A)clean(row);vector<int>used,origin(n);for(int i=0;i<n;i++)origin[i]=i;auto start=chrono::steady_clock::now();int k=0;for(;k<42;k++){int bi=-1,bj=-1,cost=1000000;for(int i=k;i<42;i++)for(int j=0;j<m;j++){if(find(used.begin(),used.end(),j)!=used.end())continue;int l=fmpq_mpoly_length(A[i][j].p,ctx);if(l&&l<cost&&unit(A[i][j])){bi=i;bj=j;cost=l;}}if(bi<0)break;swap(A[k],A[bi]);swap(origin[k],origin[bi]);used.push_back(bj);Poly pivot=A[k][bj],x,y;for(int i=k+1;i<n;i++){if(fmpq_mpoly_is_zero(A[i][bj].p,ctx))continue;Poly mult=A[i][bj];for(int j=0;j<m;j++){fmpq_mpoly_mul(x.p,pivot.p,A[i][j].p,ctx);fmpq_mpoly_mul(y.p,mult.p,A[k][j].p,ctx);fmpq_mpoly_sub(A[i][j].p,x.p,y.p,ctx);}clean(A[i]);assert(fmpq_mpoly_is_zero(A[i][bj].p,ctx));}long terms=0;for(int i=k+1;i<n;i++)for(auto&p:A[i])terms+=fmpq_mpoly_length(p.p,ctx);cerr<<"PIVOT "<<k<<" COL "<<bj<<" ROW "<<origin[k]<<" TERMS "<<terms<<" SEC "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';assert(chrono::duration<double>(chrono::steady_clock::now()-start).count()<60 && terms<3000000);}out<<n-k<<' '<<m-k<<' '<<k<<'\n';for(int j:used)out<<j<<' ';out<<'\n';for(int i=k;i<n;i++)out<<origin[i]<<' ';out<<'\n';for(int i=k;i<n;i++)for(int j=0;j<m;j++)if(find(used.begin(),used.end(),j)==used.end())write(out,A[i][j]);cerr<<"PASS "<<k<<'\n';fmpq_clear(q);known.clear();}fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
