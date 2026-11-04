// ONE exact polynomial restriction, ambient16 derivatives before substitution.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <flint/fmpq_mpoly.h>
using namespace std;
fmpq_mpoly_ctx_t ctx;
const ulong PRIME=32003;
long products=0,maxterms=0;size_t savedbytes=0;
struct P{
 fmpq_mpoly_t v;P(){fmpq_mpoly_init(v,ctx);}P(const P&a):P(){fmpq_mpoly_set(v,a.v,ctx);}
 P(P&&a):P(){fmpq_mpoly_swap(v,a.v,ctx);}P&operator=(P a){fmpq_mpoly_swap(v,a.v,ctx);return *this;}
 ~P(){fmpq_mpoly_clear(v,ctx);}
};
void check(const P&a){long n=fmpq_mpoly_length(a.v,ctx);maxterms=max(n,maxterms);if(n>150000)throw runtime_error("intermediate150k cap");}
bool zero(const P&a){return fmpq_mpoly_is_zero(a.v,ctx);}
P add(const P&a,const P&b){P c;fmpq_mpoly_add(c.v,a.v,b.v,ctx);check(c);return c;}
P mul(const P&a,const P&b){P c;fmpq_mpoly_mul(c.v,a.v,b.v,ctx);++products;check(c);return c;}
P scaleq(const P&a,const string&z){P c;fmpq_t q;fmpq_init(q);assert(fmpq_set_str(q,z.c_str(),10)==0);fmpq_canonicalise(q);fmpq_mpoly_scalar_mul_fmpq(c.v,a.v,q,ctx);fmpq_clear(q);return c;}
P scale(const P&a,long n,long d=1){return scaleq(a,to_string(n)+"/"+to_string(d));}
bool equal(const P&a,const P&b){return fmpq_mpoly_equal(a.v,b.v,ctx);}
using M=vector<P>;
M mm(const M&a,const M&b,int n){M c(n*n);for(int i=0;i<n;i++)for(int k=0;k<n;k++)if(!zero(a[n*i+k]))for(int j=0;j<n;j++)if(!zero(b[n*k+j]))c[n*i+j]=add(c[n*i+j],mul(a[n*i+k],b[n*k+j]));return c;}
P trpair(const M&a,const M&b,int n){P p;for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(!zero(a[n*i+j])&&!zero(b[n*j+i]))p=add(p,mul(a[n*i+j],b[n*j+i]));return p;}
string qstr(const fmpq_t q){char*p=fmpq_get_str(nullptr,10,q);string z(p);flint_free(p);return z;}
string zstr(const fmpz_t z){char*p=fmpz_get_str(nullptr,10,z);string t(p);flint_free(p);return t;}
void writepoly(ostream&f,const P&p){
 f<<"{\"content\":\""<<qstr(p.v->content)<<"\",\"terms\":[";ulong e[4];fmpz_t z;fmpz_init(z);
 for(long k=0;k<fmpq_mpoly_length(p.v,ctx);k++){
  fmpz_mpoly_get_term_coeff_fmpz(z,p.v->zpoly,k,ctx->zctx);fmpq_mpoly_get_term_exp_ui(e,p.v,k,ctx);
  if(k)f<<',';f<<"[\""<<zstr(z)<<"\"";for(auto x:e)f<<','<<x;f<<']';
 }f<<"]}";fmpz_clear(z);
}
void saved(const string&path){savedbytes+=filesystem::file_size(path);if(savedbytes>32*1024*1024)throw runtime_error("total32MiB cap");}
void writepolys(const string&path,const vector<P>&ps){ofstream f(path);f<<'[';for(size_t i=0;i<ps.size();i++){if(i)f<<',';writepoly(f,ps[i]);}f<<"]\n";f.close();saved(path);}
void ints(ostream&f,const vector<int>&x){f<<'[';for(size_t i=0;i<x.size();i++){if(i)f<<',';f<<x[i];}f<<']';}
ulong power(ulong a,ulong n){ulong x=1;while(n){if(n&1)x=x*a%PRIME;a=a*a%PRIME;n>>=1;}return x;}
ulong qmod(const fmpq_t q){ulong d=fmpz_fdiv_ui(fmpq_denref(q),PRIME);assert(d);return fmpz_fdiv_ui(fmpq_numref(q),PRIME)*power(d,PRIME-2)%PRIME;}
using E=array<ulong,4>;
int main(int argc,char**argv){
 assert(argc==3);auto started=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,4,ORD_DEGLEX);int status=0;{
 try{
  string dir=argv[2];ifstream in(argv[1]);int dim,width;in>>dim>>width;assert(dim==16&&width==136);vector<P>beta(16);fmpq_t z;fmpq_init(z);ulong ex[4];
  for(auto&p:beta){int count;in>>count;for(int k=0;k<count;k++){string c;in>>c;for(auto&e:ex)in>>e;assert(fmpq_set_str(z,c.c_str(),10)==0);fmpq_canonicalise(z);fmpq_mpoly_set_coeff_fmpq_ui(p.v,z,ex,ctx);}}
  vector<pair<int,int>>mon(136);for(auto&x:mon)in>>x.first>>x.second;
  vector<vector<long>>cc(32,vector<long>(136));for(auto&r:cc)for(auto&x:r)in>>x;
  vector<vector<long>>hc(512,vector<long>(136));for(auto&r:hc)for(auto&x:r)in>>x;
  vector<string>Oi(256);for(auto&x:Oi)in>>x;
  vector<int>selected(45);for(auto&x:selected)in>>x;string extra;assert(in&&!(in>>extra));
  vector<P>quad;for(auto[i,j]:mon)quad.push_back(mul(beta[i],beta[j]));
  M C(16),D(16),Theta(256),Hess(256);vector<vector<M>>der(2,vector<M>(16,M(16)));
  for(int row=0;row<32;row++)for(int k=0;k<136;k++)if(cc[row][k]){
   auto[i,j]=mon[k];int letter=row/16,pos=row%16;P term=scale(quad[k],cc[row][k],4);(letter==0?C[pos]:D[pos])=add(letter==0?C[pos]:D[pos],term);
   // Differentiate the GLOBAL beta_i beta_j, then substitute polynomial beta.
   der[letter][i][pos]=add(der[letter][i][pos],scale(beta[j],cc[row][k],4));
   der[letter][j][pos]=add(der[letter][j][pos],scale(beta[i],cc[row][k],4));
  }
  for(int row=0;row<512;row++)for(int k=0;k<136;k++)if(hc[row][k]){P term=scale(quad[k],hc[row][k],2);(row<256?Theta[row]:Hess[row-256])=add(row<256?Theta[row]:Hess[row-256],term);}
  for(int a=0;a<2;a++)for(int p=0;p<16;p++){P e;for(int k=0;k<16;k++)e=add(e,mul(beta[k],der[a][k][p]));assert(equal(e,scale(a==0?C[p]:D[p],2)));}
  vector<string>w4={"0000","0001","0011","0101","0111"},w6={"000101","000111","011011"},words=w4;words.insert(words.end(),w6.begin(),w6.end());set<string>needed;
  for(auto&w:words){int n=w.size();for(int a=0;a<n;a++){string rest=w.substr(a+1)+w.substr(0,a);for(int k=1;k<=int(rest.size());k++)needed.insert(rest.substr(0,k));}for(int half=0;half<2;half++){string h=w.substr(half*n/2,n/2);for(int k=1;k<=int(h.size());k++)needed.insert(h.substr(0,k));}}
  vector<string>order(needed.begin(),needed.end());sort(order.begin(),order.end(),[](auto&a,auto&b){return a.size()==b.size()?a<b:a.size()<b.size();});map<string,M>table;table["0"]=C;table["1"]=D;
  for(auto&w:order)if(w.size()>1)table[w]=mm(table.at(w.substr(0,w.size()-1)),w.back()=='0'?C:D,4);
  vector<P>traces;vector<vector<P>>grads;
  for(auto&w:words){int n=w.size();traces.push_back(trpair(table.at(w.substr(0,n/2)),table.at(w.substr(n/2)),4));vector<M>partial(2,M(16));
   for(int r=0;r<n;r++){string rest=w.substr(r+1)+w.substr(0,r);for(int p=0;p<16;p++)partial[w[r]-'0'][p]=add(partial[w[r]-'0'][p],table.at(rest)[p]);}
   vector<P>g(16);for(int k=0;k<16;k++)for(int a=0;a<2;a++)for(int i=0;i<4;i++)for(int j=0;j<4;j++)if(!zero(partial[a][4*i+j])&&!zero(der[a][k][4*j+i]))g[k]=add(g[k],mul(partial[a][4*i+j],der[a][k][4*j+i]));
   P euler;for(int k=0;k<16;k++)euler=add(euler,mul(beta[k],g[k]));assert(equal(euler,scale(traces.back(),2*n)));grads.push_back(move(g));
   cerr<<"TRACE_GRAD "<<w<<" TERMS "<<fmpq_mpoly_length(traces.back().v,ctx)<<endl;
  }
  vector<P>flat;for(auto&g:grads)flat.insert(flat.end(),g.begin(),g.end());writepolys(dir+"/ambient-gradients128.json",flat);writepolys(dir+"/trace-values8.json",traces);
  auto poisson=[&](int f,int g){P p;for(int i=0;i<16;i++)for(int j=0;j<16;j++)if(Oi[16*i+j]!="0"&&!zero(grads[f][i])&&!zero(grads[g][j]))p=add(p,scaleq(mul(grads[f][i],grads[g][j]),Oi[16*i+j]));return p;};
  vector<P>Q14={poisson(0,1),poisson(0,2),poisson(1,2)},Q18;for(int i=0;i<2;i++)for(int j=5;j<8;j++)Q18.push_back(poisson(i,j));
  M mu(16);for(int p=0;p<16;p++)mu[p]=add(scale(C[p],-1),scale(D[p],2));P F4=trpair(mu,mu,4);P c2=trpair(C,C,4),d2=trpair(D,D,4);assert(equal(scale(c2,5),F4)&&equal(scale(d2,10),F4)&&equal(scale(trpair(C,D,4),-10),F4));
  M D2=mm(D,D,4),C2=mm(C,C,4),CD2=mm(C,D2,4);P Pdet=add(scale(mul(d2,d2),1,8),scale(trpair(D2,D2,4),-1,4));
  P N12=add(add(trpair(CD2,CD2,4),scale(mul(d2,trpair(C2,D2,4)),-1,2)),mul(Pdet,c2));
  // Full physical16x16 Hessian value, deliberately distinct from C/D.
  P H12;for(int column=0;column<16;column++){
   vector<P>level(16);for(int i=0;i<16;i++)level[i]=Hess[16*i+column];
   for(int degree=0;degree<5;degree++){vector<P>next(16);for(int i=0;i<16;i++)for(int j=0;j<16;j++)if(!zero(Theta[16*i+j])&&!zero(level[j]))next[i]=add(next[i],mul(Theta[16*i+j],level[j]));level=move(next);}H12=add(H12,level[column]);
  }
  vector<P>Inv12;for(int i=0;i<5;i++)Inv12.push_back(mul(F4,traces[i]));Inv12.push_back(N12);Inv12.push_back(traces[5]);Inv12.push_back(traces[6]);Inv12.push_back(H12);
  vector<P>low=Q14;low.insert(low.end(),Q18.begin(),Q18.end());low.push_back(F4);low.insert(low.end(),Inv12.begin(),Inv12.end());writepolys(dir+"/low-values19.json",low);
  map<E,vector<ulong>>coeffs;vector<long>lengths;int column=0;fmpz_t zz;fmpz_init(zz);
  for(int index:selected){P p;if(index<27)p=mul(Q14[index/9],Inv12[index%9]);else p=mul(Q18[(index-27)/5],traces[(index-27)%5]);long len=fmpq_mpoly_length(p.v,ctx);assert(len<=2874&&len<=100000);lengths.push_back(len);
   ulong scalar=qmod(p.v->content);for(long k=0;k<len;k++){
    fmpq_mpoly_get_term_exp_ui(ex,p.v,k,ctx);assert(ex[1]%2==0&&(2*ex[0]+3*ex[1])%4==0);assert(ex[1]<=8&&ex[0]+ex[1]<=26&&ex[2]+ex[3]<=(26-3*ex[1])/2);
    E e={(2*ex[0]+3*ex[1])/4,ex[1]/2,ex[2],ex[3]};auto&v=coeffs[e];if(v.empty())v.resize(45);fmpz_mpoly_get_term_coeff_fmpz(zz,p.v->zpoly,k,ctx->zctx);v[column]=scalar*fmpz_fdiv_ui(zz,PRIME)%PRIME;
   }
   ofstream f(dir+"/column"+to_string(column)+".json");f<<"{\"universal_candidate_index\":"<<index<<",\"xi0_power\":13,\"HLsq_polynomial\":";writepoly(f,p);f<<"}\n";f.close();saved(dir+"/column"+to_string(column)+".json");cerr<<"COLUMN "<<column<<" ORIGINAL "<<index<<" TERMS "<<len<<endl;column++;
  }
  vector<vector<ulong>>basis;vector<int>pivots;vector<E>rows;vector<vector<ulong>>rawrows;ulong det=1;
  for(auto&[e,raw]:coeffs){vector<ulong>v=raw;for(size_t k=0;k<basis.size();k++){ulong z=v[pivots[k]];if(z)for(int j=pivots[k];j<45;j++)v[j]=(v[j]+PRIME-z*basis[k][j]%PRIME)%PRIME;}
   int pivot=0;while(pivot<45&&!v[pivot])pivot++;if(pivot==45)continue;ulong lead=v[pivot];det=det*lead%PRIME;ulong inv=power(lead,PRIME-2);for(int j=pivot;j<45;j++)v[j]=v[j]*inv%PRIME;basis.push_back(v);pivots.push_back(pivot);rows.push_back(e);rawrows.push_back(raw);if(basis.size()==45)break;
  }
  // Sort rank rows' pivot order is automatically increasing only for row-echelon insertion;
  // the reduction above finds each FIRST free pivot, which can be out of order.
  // Recompute an ordinary dense elimination determinant/rank on saved independent rows.
  auto dense=rawrows;int rank=0;ulong determinant=1;
  for(int col=0;col<45&&rank<int(dense.size());col++){int p=rank;while(p<int(dense.size())&&!dense[p][col])p++;if(p==int(dense.size()))continue;if(p!=rank){swap(dense[p],dense[rank]);determinant=(PRIME-determinant)%PRIME;}ulong lead=dense[rank][col];determinant=determinant*lead%PRIME;ulong inv=power(lead,PRIME-2);for(int j=col;j<45;j++)dense[rank][j]=dense[rank][j]*inv%PRIME;for(int k=rank+1;k<int(dense.size());k++){ulong c=dense[k][col];if(c)for(int j=col;j<45;j++)dense[k][j]=(dense[k][j]+PRIME-c*dense[rank][j]%PRIME)%PRIME;}rank++;}
  assert(rank==int(rawrows.size()));
  ofstream out(dir+"/native-certificate.json");out<<"{\"all_passed\":true,\"prime\":32003,\"rank\":"<<rank<<",\"coefficient_minor_mod32003\":"<<(rank==45?determinant:0)<<",\"selected_universal_candidates\":";ints(out,selected);out<<",\"coefficient_monomials_r0_B0_s0_q0\":[";
  for(size_t i=0;i<rows.size();i++){if(i)out<<',';out<<'[';for(int j=0;j<4;j++){if(j)out<<',';out<<rows[i][j];}out<<']';}out<<"],\"selected_coefficient_matrix\":[";for(size_t i=0;i<rawrows.size();i++){if(i)out<<',';out<<'[';for(int j=0;j<45;j++){if(j)out<<',';out<<rawrows[i][j];}out<<']';}out<<"],\"column_term_counts\":[";for(size_t i=0;i<lengths.size();i++){if(i)out<<',';out<<lengths[i];}out<<"],\"all32_quadratic_and_all8_trace_ambient_Euler\":true,\"all45_finite_parity_and_support\":true,\"exactQQ_columns_saved\":45,\"max_intermediate_terms\":"<<maxterms<<",\"products\":"<<products<<",\"bytes_saved\":"<<savedbytes<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"}\n";out.close();saved(dir+"/native-certificate.json");fmpz_clear(zz);fmpq_clear(z);cerr<<"FINAL_RANK "<<rank<<" MINOR "<<determinant<<endl;
 }catch(const exception&e){cerr<<"BOUNDED_SOURCE26_STOP "<<e.what()<<endl;status=2;}
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();return status;
}
