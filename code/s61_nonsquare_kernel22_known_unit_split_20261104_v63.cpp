// One bounded known-unit Gaussian job on the frozen174x22 source kernel.
#include "s61_binary_hyperplane_cubic_qq3_20261102.hpp"
#include <chrono>
#include <climits>
#include <cstdio>
#include <functional>
using Matrix=vector<vector<Rat>>;
vector<Rat> units;
auto began=chrono::steady_clock::now();
double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-began).count();}
long terms(const Matrix&B){long z=0;for(const auto&r:B)for(const auto&c:r)if(!c.zero())z+=fmpq_mpoly_length(c.n,ctx)+fmpq_mpoly_length(c.d,ctx);return z;}
bool knownunit(const Rat&c,vector<int>&exps,string&scalar){
 if(c.zero())return false;
 fmpq_mpoly_t p,t;fmpq_mpoly_init(p,ctx);fmpq_mpoly_init(t,ctx);fmpq_t n,d;fmpq_init(n);fmpq_init(d);exps.assign(6,0);bool ok=true;
 for(int side=0;side<2&&ok;side++){
  fmpq_mpoly_set(p,side?c.d:c.n,ctx);
  for(int j=0;j<6;j++)while(fmpq_mpoly_divides(t,p,units[j].n,ctx)){fmpq_mpoly_swap(p,t,ctx);exps[j]+=side?-1:1;}
  if(!fmpq_mpoly_is_fmpq(p,ctx)){ok=false;break;}
  fmpq_mpoly_get_fmpq(side?d:n,p,ctx);
 }
 if(ok){fmpq_div(n,n,d);char*x=fmpq_get_str(nullptr,10,n);scalar=x;flint_free(x);}
 fmpq_clear(n);fmpq_clear(d);fmpq_mpoly_clear(p,ctx);fmpq_mpoly_clear(t,ctx);return ok;
}
Rat one(){Rat c;fmpq_mpoly_one(c.n,ctx);return c;}
void writerat(ostream&o,const Rat&c){writepoly(o,c.n);writepoly(o,c.d);}
void write_state(const string&path,const Matrix&O,const Matrix&B,const Matrix&LC,
 const vector<int>&rr,const vector<int>&cc,const string&status,bool verified){
 string tmp=path+".tmp";ofstream f(tmp);int n=O.size(),m=O[0].size(),r=rr.size();
 f<<n<<' '<<m<<' '<<r<<'\n';for(int i:rr)f<<i<<' ';f<<'\n';for(int j:cc)f<<j<<' ';f<<'\n';
 f<<status<<' '<<verified<<' '<<elapsed()<<' '<<terms(B)<<'\n';
 for(int i=0;i<r;i++)for(int j=0;j<r;j++)writerat(f,O[rr[i]][cc[j]]);
 for(int i=0;i<r;i++)for(int j=0;j<r;j++){Rat z;if(i==j)z=one();else if(i>j)z=LC[rr[i]][j];writerat(f,z);}
 for(int i=0;i<r;i++)for(int j=0;j<r;j++){Rat z;if(i<=j)z=B[rr[i]][cc[j]];writerat(f,z);}
 for(const auto&row:B)for(const auto&c:row)writerat(f,c);
 f.close();ifstream sizecheck(tmp,ios::binary|ios::ate);assert(sizecheck.tellg()<32*1024*1024);sizecheck.close();assert(!rename(tmp.c_str(),path.c_str()));
}
int main(int argc,char**argv){
 assert(argc==4);fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);{
 const char*names[]={"a","e","d"};
 for(const char*t:{"a","d","500*a+90*e^2+3*e-24","50*a*d+18*d*e^2-9*d*e-2*d+60*e-40","100*a+18*e^2-9*e-2","5*a*d-6*e-4"}){Rat q;assert(!fmpq_mpoly_set_str_pretty(q.n,t,names,ctx));units.push_back(q);}
 ifstream in(argv[1]);ofstream log(argv[3]);int n,m;in>>n>>m;assert(n==174&&m==22);Matrix O(n,vector<Rat>(m));
 for(auto&row:O)for(auto&c:row){readpoly(in,c.n);readpoly(in,c.d);assert(!fmpq_mpoly_is_zero(c.d,ctx));}
 string extra;assert(!(in>>extra));Matrix B=O,LC(n,vector<Rat>(m));vector<int>rr,cc;vector<bool>ur(n),uc(m);string status="started";
 write_state(argv[2],O,B,LC,rr,cc,status,false);
 while((int)rr.size()<m){
  int bi=-1,bj=-1;long best=LONG_MAX;vector<int>be;string bs;
  for(int i=0;i<n;i++)if(!ur[i])for(int j=0;j<m;j++)if(!uc[j]&&!B[i][j].zero()){
   long cost=fmpq_mpoly_length(B[i][j].n,ctx)+fmpq_mpoly_length(B[i][j].d,ctx);if(cost>=best)continue;
   vector<int>e;string s;if(knownunit(B[i][j],e,s)){bi=i;bj=j;best=cost;be=e;bs=s;}
   if(elapsed()>35)break;
  }
  if(elapsed()>35){status="search_time_frontier";break;}
  if(bi<0){status="no_individual_known_unit_entry";break;}
  Rat pivot=B[bi][bj];Matrix next=B;vector<Rat>newL(n);bool abort=false;
  for(int i=0;i<n&&!abort;i++)if(!ur[i]&&i!=bi&&!B[i][bj].zero()){
   newL[i]=divi(B[i][bj],pivot);
   for(int j=0;j<m;j++)if(!uc[j]&&j!=bj&&!B[bi][j].zero()){
    next[i][j]=sub(B[i][j],mul(newL[i],B[bi][j]));
    if(elapsed()>35){abort=true;break;}
   }
   next[i][bj]=Rat();
  }
  if(abort){status="update_time_frontier_previous_pivot_retained";break;}
  long nt=terms(next);if(nt>300000){status="term_frontier_previous_pivot_retained";break;}
  int r=rr.size();for(int i=0;i<n;i++)if(!ur[i]&&i!=bi)LC[i][r]=newL[i];
  B=move(next);rr.push_back(bi);cc.push_back(bj);ur[bi]=uc[bj]=true;status="unit_pivot_checkpoint";
  log<<"PIVOT "<<r<<" ROW "<<bi<<" COL "<<bj<<" TERMS "<<nt<<" SEC "<<elapsed()<<" SCALAR "<<bs<<" EXPONENTS";for(int z:be)log<<' '<<z;log<<'\n';log.flush();
  write_state(argv[2],O,B,LC,rr,cc,status,false);
 }
 bool verified=true;int r=rr.size();
 for(int i=0;i<r&&verified;i++)for(int j=0;j<r&&verified;j++){
  Rat sum;
  for(int k=0;k<=min(i,j);k++){Rat l=(i==k)?one():LC[rr[i]][k];sum=sub(sum,sub(Rat(),mul(l,B[rr[k]][cc[j]])));}
  if(!sub(sum,O[rr[i]][cc[j]]).zero()){verified=false;log<<"A_LU_MISMATCH "<<i<<' '<<j<<'\n';}
  if(elapsed()>52){verified=false;status="LU_verification_time_frontier";break;}
 }
 if(r==22&&verified)status="FULL_GLOBAL_LEFT_INVERSE_CIRCUIT_CERTIFIED";
 log<<"FINAL_RANK "<<r<<" A_EQUALS_LU "<<verified<<" STATUS "<<status<<" SEC "<<elapsed()<<'\n';log.flush();
 write_state(argv[2],O,B,LC,rr,cc,status,verified);assert(verified);
 units.clear();
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();
}
