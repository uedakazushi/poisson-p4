// Fraction-free exact rank with core-only row contents and pivots.
#include <cassert>
#include <fstream>
#include <iostream>
#include <chrono>
#include <map>
#include <vector>
#include <string>
#include <flint/fmpq_mpoly.h>
using namespace std;
fmpq_mpoly_ctx_t ctx;int nv;
struct P {fmpq_mpoly_struct p[1];P(){fmpq_mpoly_init(p,ctx);}P(const P&o):P(){fmpq_mpoly_set(p,o.p,ctx);}P(P&&o):P(){fmpq_mpoly_swap(p,o.p,ctx);}P&operator=(P o){fmpq_mpoly_swap(p,o.p,ctx);return *this;}~P(){fmpq_mpoly_clear(p,ctx);}bool zero()const{return fmpq_mpoly_is_zero(p,ctx);}long size()const{return fmpq_mpoly_length(p,ctx);}};
void read(istream&in,P&a){int n;in>>n;vector<ulong>ex(nv);fmpq_t c;fmpq_init(c);for(int k=0;k<n;k++){string v;in>>v;assert(fmpq_set_str(c,v.c_str(),10)==0);fmpq_canonicalise(c);for(auto&e:ex)in>>e;fmpq_mpoly_push_term_fmpq_ui(a.p,c,ex.data(),ctx);}fmpq_mpoly_sort_terms(a.p,ctx);fmpq_mpoly_combine_like_terms(a.p,ctx);fmpq_clear(c);assert(in);}
void write(ostream&out,const P&a){out<<a.size()<<'\n';vector<ulong>ex(nv);fmpq_t c;fmpq_init(c);for(int k=0;k<a.size();k++){fmpq_mpoly_get_term_coeff_fmpq(c,a.p,k,ctx);char*t=fmpq_get_str(nullptr,10,c);out<<t;flint_free(t);fmpq_mpoly_get_term_exp_ui(ex.data(),a.p,k,ctx);for(auto e:ex)out<<' '<<e;out<<'\n';}fmpq_clear(c);}
bool core(const P&a){vector<ulong>ex(nv);for(int k=0;k<a.size();k++){fmpq_mpoly_get_term_exp_ui(ex.data(),a.p,k,ctx);for(int i=0;i<nv;i++)if(i!=0&&i!=3&&ex[i])return false;}return true;}
using Row=map<int,P>;
vector<P> removed_factors;
void normalize(Row&row){if(row.empty())return;auto it=row.begin();for(auto k=row.begin();k!=row.end();++k)if(core(k->second)){it=k;break;}fmpq_t c;fmpq_init(c);fmpq_mpoly_get_term_coeff_fmpq(c,it->second.p,0,ctx);for(auto&[j,p]:row)fmpq_mpoly_scalar_div_fmpq(p.p,p.p,c,ctx);fmpq_clear(c);}
void primitive(Row&row){
 if(row.empty())return;P g,tmp;bool started=false;
 for(auto&[j,p]:row)if(core(p)){g=p;started=true;break;}
 if(started){for(auto&[j,p]:row){assert(fmpq_mpoly_gcd(tmp.p,g.p,p.p,ctx));g=tmp;if(fmpq_mpoly_is_fmpq(g.p,ctx)){normalize(row);return;}}}
 else{
  vector<ulong>ex(nv);fmpq_t co;fmpq_init(co);
  for(auto&[j,p]:row){map<vector<ulong>,P>coeffs;for(int k=0;k<p.size();k++){fmpq_mpoly_get_term_exp_ui(ex.data(),p.p,k,ctx);fmpq_mpoly_get_term_coeff_fmpq(co,p.p,k,ctx);auto key=ex;key[0]=key[3]=0;for(int i=0;i<nv;i++)if(i!=0&&i!=3)ex[i]=0;fmpq_mpoly_push_term_fmpq_ui(coeffs[key].p,co,ex.data(),ctx);}for(auto&[key,c]:coeffs){fmpq_mpoly_sort_terms(c.p,ctx);fmpq_mpoly_combine_like_terms(c.p,ctx);if(!started){g=c;started=true;}else{assert(fmpq_mpoly_gcd(tmp.p,g.p,c.p,ctx));g=tmp;}if(fmpq_mpoly_is_fmpq(g.p,ctx)){fmpq_clear(co);normalize(row);return;}}}fmpq_clear(co);
 }
 assert(started&&core(g));removed_factors.push_back(g);for(auto&[j,p]:row){assert(fmpq_mpoly_divides(tmp.p,p.p,g.p,ctx));p=tmp;}normalize(row);
}
int main(int argc,char**argv){assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);assert(in&&out);int n;in>>nv>>n;assert(nv==8);fmpq_mpoly_ctx_init(ctx,nv,ORD_LEX);flint_set_num_threads(2);map<int,Row>rows;for(int k=0;k<n;k++){int i,m;in>>i>>m;for(int l=0;l<m;l++){int j;in>>j;read(in,rows[i][j]);}}
 auto st=chrono::steady_clock::now();auto secs=[&](){return chrono::duration<double>(chrono::steady_clock::now()-st).count();};for(auto&[i,row]:rows)primitive(row);cerr<<"ready "<<secs()<<'\n';vector<pair<pair<int,int>,Row>>basis;
 while(!rows.empty()){
  long best=LONG_MAX;int ri=-1,cj=-1;
  for(auto&[i,row]:rows)for(auto&[j,p]:row)if(core(p)&&p.size()<best){best=p.size();ri=i;cj=j;}
  if(ri<0){cerr<<"NO CORE PIVOT "<<rows.size()<<'\n';break;}
  Row row=std::move(rows[ri]);rows.erase(ri);P pivot=row[cj];cerr<<"begin "<<basis.size()+1<<' '<<ri<<' '<<cj<<' '<<secs()<<'\n';
  for(auto it=rows.begin();it!=rows.end();){Row&rr=it->second;auto pos=rr.find(cj);if(pos!=rr.end()){P c=pos->second;Row newer;for(auto&[j,p]:rr)fmpq_mpoly_mul(newer[j].p,p.p,pivot.p,ctx);for(auto&[j,p]:row){P prod;fmpq_mpoly_mul(prod.p,c.p,p.p,ctx);fmpq_mpoly_sub(newer[j].p,newer[j].p,prod.p,ctx);}for(auto q=newer.begin();q!=newer.end();)if(q->second.zero())q=newer.erase(q);else++q;primitive(newer);rr=std::move(newer);}if(rr.empty())it=rows.erase(it);else++it;}
  basis.push_back({{ri,cj},std::move(row)});cerr<<"end "<<basis.size()<<" rows "<<rows.size()<<" seconds "<<secs()<<'\n';
 }
 out<<basis.size()<<' '<<rows.size()<<' '<<removed_factors.size()<<'\n';
 for(auto&b:basis){out<<b.first.first<<' '<<b.first.second<<'\n';write(out,b.second.at(b.first.second));}
 for(auto&g:removed_factors)write(out,g);

 cerr<<"FINAL RANK "<<(rows.empty()?int(basis.size()):-1)<<" seconds "<<secs()<<'\n';assert(out);
}
