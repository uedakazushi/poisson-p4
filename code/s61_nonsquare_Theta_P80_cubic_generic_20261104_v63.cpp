// ONE source-bound univariate substitution and exact fraction-free rank test.
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_mat.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>
#include <unistd.h>
using namespace std;
struct Z{fmpz_t v;Z(){fmpz_init(v);}Z(const Z&x):Z(){fmpz_set(v,x.v);}Z(Z&&x):Z(){fmpz_swap(v,x.v);}Z&operator=(Z x){fmpz_swap(v,x.v);return *this;}~Z(){fmpz_clear(v);}};
struct P{fmpz_poly_t v;P(){fmpz_poly_init(v);}P(const P&x):P(){fmpz_poly_set(v,x.v);}P(P&&x):P(){fmpz_poly_swap(v,x.v);}P&operator=(P x){fmpz_poly_swap(v,x.v);return *this;}~P(){fmpz_poly_clear(v);}};
P product(const P&a,const P&b){P p;fmpz_poly_mul(p.v,a.v,b.v);return p;}
P power(const P&a,long k){P p;fmpz_poly_pow(p.v,a.v,k);return p;}
string zs(const fmpz_t x){char*s=fmpz_get_str(nullptr,10,x);string q(s);flint_free(s);return q;}
void jp(ostream&o,const fmpz_poly_t p){o<<'[';bool first=true;for(long i=0;i<fmpz_poly_length(p);i++){auto c=fmpz_poly_get_coeff_ptr(p,i);if(fmpz_is_zero(c))continue;if(!first)o<<',';first=false;o<<'['<<i<<",\""<<zs(c)<<"\"]";}o<<']';}
long terms(const fmpz_poly_t p){long n=0;for(long i=0;i<fmpz_poly_length(p);i++)n+=!fmpz_is_zero(fmpz_poly_get_coeff_ptr(p,i));return n;}
struct Term{Z c;int i,j,k;};
int main(int argc,char**argv){assert(argc==4);alarm(55);auto begin=chrono::steady_clock::now();ifstream in(argv[1]);string prefix=argv[2],resultfile=argv[3];
 P a,e,d,b;fmpz_poly_set_coeff_si(a.v,0,-640);fmpz_poly_set_coeff_si(a.v,1,-56);fmpz_poly_set_coeff_si(e.v,0,80);fmpz_poly_set_coeff_si(e.v,1,2);fmpz_poly_set_coeff_si(d.v,0,-60);fmpz_poly_set_coeff_si(d.v,1,-3);fmpz_poly_set_coeff_si(b.v,1,300);fmpz_poly_set_coeff_si(b.v,2,7);
 vector<P>ap(14),ep(24),dp(7);for(int i=0;i<14;i++)ap[i]=power(a,i);for(int j=0;j<24;j++)ep[j]=power(e,j);for(int k=0;k<7;k++)dp[k]=power(d,k);
 vector<Z>ca(14),ce(24),cd(7);for(int i=0;i<14;i++)fmpz_ui_pow_ui(ca[i].v,45,i);for(int j=0;j<24;j++)fmpz_ui_pow_ui(ce[j].v,9,j);for(int k=0;k<7;k++)fmpz_ui_pow_ui(cd[k].v,80,k);
 map<tuple<int,int,int>,P>cache;long cachedterms=0;
 auto get=[&](int i,int j,int k)->const P&{auto key=make_tuple(i,j,k);auto it=cache.find(key);if(it==cache.end()){P p=product(product(ap[i],ep[j]),dp[k]);cachedterms+=terms(p.v);assert(cachedterms<=150000);it=cache.emplace(key,move(p)).first;}return it->second;};
 long rows,cols;in>>rows>>cols;assert(rows==17&&cols==24);fmpz_poly_mat_t A,B;fmpz_poly_mat_init(A,rows,cols);fmpz_poly_mat_init(B,rows,cols);
 vector<Z>rowlcms(rows),rowcontents(rows);long inputterms=0,maxdeg=0,maxbits=0;vector<vector<Z>> entrydens(rows,vector<Z>(cols));
 for(long r=0;r<rows;r++){
  vector<P>numer(cols);vector<Z>dens(cols);fmpz_one(rowlcms[r].v);
  for(long c=0;c<cols;c++){
   long count,kind;string scale;in>>kind>>scale>>count;assert(in&&kind>=0&&kind<=2&&count>=0);vector<Term>ts(count);int I=0,J=0,K=0;
   for(auto&t:ts){string s;in>>s>>t.i>>t.j>>t.k;assert(!fmpz_set_str(t.c.v,s.c_str(),10));assert(t.i>=0&&t.i<14&&t.j>=0&&t.j<24&&t.k>=0&&t.k<7);I=max(I,t.i);J=max(J,t.j);K=max(K,t.k);}
   fmpz_mul(dens[c].v,ca[I].v,ce[J].v);fmpz_mul(dens[c].v,dens[c].v,cd[K].v);P temp;Z co;
   for(auto&t:ts){fmpz_mul(co.v,t.c.v,ca[I-t.i].v);fmpz_mul(co.v,co.v,ce[J-t.j].v);fmpz_mul(co.v,co.v,cd[K-t.k].v);
    auto&v=get(t.i,t.j,t.k);int shift=24+t.k-2*t.i-t.j;assert(shift>=0);fmpz_poly_scalar_mul_fmpz(temp.v,v.v,co.v);fmpz_poly_shift_left(temp.v,temp.v,shift);fmpz_poly_add(numer[c].v,numer[c].v,temp.v);}
   // common scalar is zeta^24*a*B, and original q=scale*(1,a,B).
   if(kind==0||kind==2){fmpz_poly_mul(numer[c].v,numer[c].v,a.v);assert((fmpz_poly_length(numer[c].v)<1||fmpz_is_zero(fmpz_poly_get_coeff_ptr(numer[c].v,0)))&&(fmpz_poly_length(numer[c].v)<2||fmpz_is_zero(fmpz_poly_get_coeff_ptr(numer[c].v,1))));fmpz_poly_shift_right(numer[c].v,numer[c].v,2);fmpz_mul_ui(dens[c].v,dens[c].v,45);}
   if(kind==0||kind==1){fmpz_poly_mul(numer[c].v,numer[c].v,b.v);fmpz_mul_ui(dens[c].v,dens[c].v,60);}
   Z qcontent,common;assert(!fmpz_set_str(qcontent.v,scale.c_str(),10));assert(!fmpz_is_zero(qcontent.v));fmpz_mul(dens[c].v,dens[c].v,qcontent.v);
   if(fmpz_sgn(dens[c].v)<0){fmpz_neg(dens[c].v,dens[c].v);fmpz_poly_neg(numer[c].v,numer[c].v);}
   if(fmpz_poly_is_zero(numer[c].v)){fmpz_one(dens[c].v);}else{fmpz_poly_content(qcontent.v,numer[c].v);fmpz_gcd(common.v,qcontent.v,dens[c].v);fmpz_poly_scalar_divexact_fmpz(numer[c].v,numer[c].v,common.v);fmpz_divexact(dens[c].v,dens[c].v,common.v);}
   entrydens[r][c]=dens[c];fmpz_lcm(rowlcms[r].v,rowlcms[r].v,dens[c].v);
  }
  for(long c=0;c<cols;c++){Z factor,content;fmpz_divexact(factor.v,rowlcms[r].v,dens[c].v);fmpz_poly_scalar_mul_fmpz(fmpz_poly_mat_entry(A,r,c),numer[c].v,factor.v);fmpz_poly_content(content.v,fmpz_poly_mat_entry(A,r,c));fmpz_gcd(rowcontents[r].v,rowcontents[r].v,content.v);}
  assert(!fmpz_is_zero(rowcontents[r].v));
  for(long c=0;c<cols;c++){auto p=fmpz_poly_mat_entry(A,r,c);fmpz_poly_scalar_divexact_fmpz(p,p,rowcontents[r].v);inputterms+=terms(p);maxdeg=max(maxdeg,fmpz_poly_degree(p));for(long j=0;j<fmpz_poly_length(p);j++)maxbits=max(maxbits,long(fmpz_bits(fmpz_poly_get_coeff_ptr(p,j))));}
 }
 string extra;assert(!(in>>extra));assert(inputterms<=150000&&maxdeg<=34);
 {ofstream o(prefix);o<<"{\"variable\":\"zeta\",\"common_source_scalar\":\"zeta^24*a*B\",\"row_LCMs\":[";for(int i=0;i<rows;i++){if(i)o<<',';o<<'\"'<<zs(rowlcms[i].v)<<'\"';}o<<"],\"row_integer_contents\":[";for(int i=0;i<rows;i++){if(i)o<<',';o<<'\"'<<zs(rowcontents[i].v)<<'\"';}o<<"],\"integer_matrix17x24\":[";
  for(int i=0;i<rows;i++){if(i)o<<',';o<<'[';for(int j=0;j<cols;j++){if(j)o<<',';jp(o,fmpz_poly_mat_entry(A,i,j));}o<<']';}o<<"],\"ordinary_polynomial_terms\":"<<inputterms<<",\"max_degree\":"<<maxdeg<<",\"max_coefficient_bits\":"<<maxbits<<",\"cache_terms\":"<<cachedterms<<"}\n";}
 cout<<"SOURCE17x24_CHECKPOINT terms="<<inputterms<<" degree="<<maxdeg<<" bits="<<maxbits<<endl;
 vector<slong>perm(rows);for(int i=0;i<rows;i++)perm[i]=i;P den;slong rank=fmpz_poly_mat_fflu(B,den.v,perm.data(),A,0);assert(rank>=0&&rank<=rows);
 vector<long>piv;long next=0;for(int i=0;i<rank;i++){while(next<cols&&fmpz_poly_is_zero(fmpz_poly_mat_entry(B,i,next)))next++;assert(next<cols);piv.push_back(next++);}
 fmpz_poly_mat_t minor;fmpz_poly_mat_init(minor,rank,rank);for(int i=0;i<rank;i++)for(int j=0;j<rank;j++)fmpz_poly_set(fmpz_poly_mat_entry(minor,i,j),fmpz_poly_mat_entry(A,perm[i],piv[j]));
 P det,neg;fmpz_poly_mat_det_fflu(det.v,minor);assert(!fmpz_poly_is_zero(det.v));fmpz_poly_neg(neg.v,den.v);assert(fmpz_poly_equal(det.v,den.v)||fmpz_poly_equal(det.v,neg.v));
 {ofstream o(resultfile);o<<"{\"exact_QQ_zeta_rank\":"<<rank<<",\"selected_local_rows\":[";for(int i=0;i<rank;i++){if(i)o<<',';o<<perm[i];}o<<"],\"selected_local_columns\":[";for(int j=0;j<rank;j++){if(j)o<<',';o<<piv[j];}o<<"],\"nonzero_selected_minor\":";jp(o,det.v);o<<",\"fflu_denominator\":";jp(o,den.v);o<<",\"minor_degree\":"<<fmpz_poly_degree(det.v)<<",\"full_row_rank\":"<<(rank==rows?"true":"false")<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-begin).count()<<"}\n";}
 cout<<"EXACT_RANK "<<rank<<endl;fmpz_poly_mat_clear(minor);fmpz_poly_mat_clear(B);fmpz_poly_mat_clear(A);flint_cleanup();
}
