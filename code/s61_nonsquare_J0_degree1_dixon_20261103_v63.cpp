// Modular-guided single-target Dixon variant; exact QQ acceptance only.
// Degree1 polynomial-cofactor extension of frozenV9 primitive-entry span.
// Original polynomial unit identity is verified after every scalar/unit lift.
// Args: input output unit-polynomial-file has-metadata(0/1) [transpose] [colclean] [witness-file target-row].
// Witness: column-count, then one sparse polynomial per column in the same 3-variable format.
// target-row is zero-based AFTER the optional transpose. The exact resulting entry must be a unit.
// Three variables are a,e,d. Input: n m metadata_count, optional metadata_count+n
// integer labels, then n*m sparse polynomials: term_count,(coefficient,aexp,eexp,dexp)*.
// Units are one or more nonconstant polynomials inverted on the intended open set.
// Only nonzero QQ constants and these units are inverted. Added span pivots use
// invertible constant row/column operations, verified by exact target equality.
// Output preserves the sparse format with metadata (pivot columns, remaining rows).
// A partial pivot count is NOT a full-rank certificate. Call under a time limit.
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/nmod_mat.h>
#include <flint/nmod.h>
#include <map>
#include <functional>
#include <array>
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


int primitiveTargetFactorBound=8;vector<Poly> known;vector<Poly>primitiveWitness,primitiveTarget;
bool unit(Poly p){if(fmpq_mpoly_is_zero(p.p,ctx))return false;Poly q;for(auto&f:known)while(fmpq_mpoly_divides(q.p,p.p,f.p,ctx))p=q;return fmpq_mpoly_is_fmpq(p.p,ctx);}
void clean(vector<Poly>&row){int first=-1;for(int j=0;j<(int)row.size();j++)if(!fmpq_mpoly_is_zero(row[j].p,ctx)){first=j;break;}if(first<0)return;fmpq_t c;fmpq_init(c);fmpq_mpoly_get_term_coeff_fmpq(c,row[first].p,0,ctx);for(auto&p:row)fmpq_mpoly_scalar_div_fmpq(p.p,p.p,c,ctx);fmpq_clear(c);vector<Poly> quot(row.size());for(auto&f:known)while(true){bool yes=true;for(int j=first;j<(int)row.size();j++)if(!fmpq_mpoly_divides(quot[j].p,row[j].p,f.p,ctx)){yes=false;break;}if(!yes)break;row=quot;}}

// A modular rank choice is only a search heuristic.  Every accepted solution
// is reconstructed over QQ with Dixon and checked in the complete polynomial row.
bool combo(vector<vector<Poly>>&A,int k,const vector<int>&used,int &bi,int &bj){
 assert(A.size()==1&&k==0&&used.empty());int nc=A[0].size();map<array<ulong,3>,int>mons;vector<array<ulong,3>>monkeys;ulong ex[3];long maxdeg=0;
 for(auto&p:A[0]){maxdeg=max(maxdeg,fmpq_mpoly_total_degree_si(p.p,ctx));for(int t=0;t<fmpq_mpoly_length(p.p,ctx);t++){fmpq_mpoly_get_term_exp_ui(ex,p.p,t,ctx);array<ulong,3>v={ex[0],ex[1],ex[2]};if(!mons.count(v)){mons[v]=mons.size();monkeys.push_back(v);}}}
 int nr=mons.size();const ulong prime=2147483647; nmod_t mod;nmod_init(&mod,prime);
 fmpq_mat_t M;fmpq_mat_init(M,nr,nc);nmod_mat_t N,E;nmod_mat_init(N,nr,nc,prime);nmod_mat_init(E,nr,nc,prime);fmpq_t q;fmpq_init(q);
 auto reduce=[&](const fmpq_t v)->ulong{ulong d=fmpz_fdiv_ui(fmpq_denref(v),prime);assert(d);return nmod_mul(fmpz_fdiv_ui(fmpq_numref(v),prime),nmod_inv(d,mod),mod);};
 for(int j=0;j<nc;j++)for(int t=0;t<fmpq_mpoly_length(A[0][j].p,ctx);t++){fmpq_mpoly_get_term_exp_ui(ex,A[0][j].p,t,ctx);fmpq_mpoly_get_term_coeff_fmpq(q,A[0][j].p,t,ctx);int i=mons.at({ex[0],ex[1],ex[2]});fmpq_set(fmpq_mat_entry(M,i,j),q);nmod_mat_entry(N,i,j)=reduce(q);}
 nmod_mat_set(E,N);int rr=nmod_mat_rref(E);cerr<<"MODULAR_SEARCH_RANK "<<rr<<" MONS "<<nr<<" COLS "<<nc<<endl;vector<int>cols;
 for(int i=0;i<rr;i++){int j=0;while(j<nc&&!nmod_mat_entry(E,i,j))j++;assert(j<nc);cols.push_back(j);}
 nmod_mat_t MT; nmod_mat_init(MT,rr,nr,prime);for(int i=0;i<rr;i++)for(int j=0;j<nr;j++)nmod_mat_entry(MT,i,j)=nmod_mat_entry(N,j,cols[i]);assert(nmod_mat_rref(MT)==rr);vector<int>rows;
 for(int i=0;i<rr;i++){int j=0;while(j<nr&&!nmod_mat_entry(MT,i,j))j++;assert(j<nr);rows.push_back(j);}
 nmod_mat_t B,BI,b,x;nmod_mat_init(B,rr,rr,prime);nmod_mat_init(BI,rr,rr,prime);nmod_mat_init(b,rr,1,prime);nmod_mat_init(x,rr,1,prime);
 fmpq_mat_t QB,Qb,Qx;fmpq_mat_init(QB,rr,rr);fmpq_mat_init(Qb,rr,1);fmpq_mat_init(Qx,rr,1);
 for(int i=0;i<rr;i++)for(int j=0;j<rr;j++){nmod_mat_entry(B,i,j)=nmod_mat_entry(N,rows[i],cols[j]);fmpq_set(fmpq_mat_entry(QB,i,j),fmpq_mat_entry(M,rows[i],cols[j]));}assert(nmod_mat_inv(BI,B));
 long tested=0,eligible=0,modpassed=0,exactsolves=0;
 auto test=[&](const Poly&t)->bool{
  tested++;for(int h=0;h<fmpq_mpoly_length(t.p,ctx);h++){fmpq_mpoly_get_term_exp_ui(ex,t.p,h,ctx);if(!mons.count({ex[0],ex[1],ex[2]}))return false;}eligible++;
  for(int i=0;i<rr;i++){auto v=monkeys[rows[i]];fmpq_mpoly_get_coeff_fmpq_ui(fmpq_mat_entry(Qb,i,0),t.p,v.data(),ctx);nmod_mat_entry(b,i,0)=reduce(fmpq_mat_entry(Qb,i,0));}nmod_mat_mul(x,BI,b);
  for(int i=0;i<nr;i++){ulong s=0;for(int j=0;j<rr;j++)s=nmod_addmul(s,nmod_mat_entry(N,i,cols[j]),nmod_mat_entry(x,j,0),mod);auto v=monkeys[i];fmpq_mpoly_get_coeff_fmpq_ui(q,t.p,v.data(),ctx);if(s!=reduce(q))return false;}
  modpassed++;cerr<<"MODULAR_TARGET_PASS "<<modpassed<<" EXACT_DIXON_START RANK "<<rr<<endl;exactsolves++;
  assert(fmpq_mat_solve_dixon(Qx,QB,Qb));Poly sum,tmp;
  for(int i=0;i<rr;i++){fmpq_mpoly_scalar_mul_fmpq(tmp.p,A[0][cols[i]].p,fmpq_mat_entry(Qx,i,0),ctx);fmpq_mpoly_add(sum.p,sum.p,tmp.p,ctx);}
  if(!fmpq_mpoly_equal(sum.p,t.p,ctx)){cerr<<"MODULAR_TARGET_FALSE_POSITIVE_EXCLUDED\n";return false;}
  primitiveWitness=vector<Poly>(nc);for(int i=0;i<rr;i++)fmpq_mpoly_set_fmpq(primitiveWitness[cols[i]].p,fmpq_mat_entry(Qx,i,0),ctx);primitiveTarget={t};bi=0;bj=cols[0];cerr<<"EXACT_DIXON_FULL_POLYNOMIAL_IDENTITY_PASS\n";return true;
 };
 Poly one;fmpq_mpoly_one(one.p,ctx);bool found=test(one);vector<long>degrees;for(auto&f:known)degrees.push_back(fmpq_mpoly_total_degree_si(f.p,ctx));
 function<bool(int,int,long,const Poly&)>extend=[&](int pos,int left,long degree,const Poly&prod)->bool{
  if(pos==(int)known.size()-1){if(degree+left*degrees[pos]>maxdeg)return false;Poly v=prod,tmp;for(int z=0;z<left;z++){fmpq_mpoly_mul(tmp.p,v.p,known[pos].p,ctx);v=tmp;}return test(v);}
  Poly v=prod,tmp;for(int z=0;z<=left&&degree+z*degrees[pos]<=maxdeg;z++){if(extend(pos+1,left-z,degree+z*degrees[pos],v))return true;fmpq_mpoly_mul(tmp.p,v.p,known[pos].p,ctx);v=tmp;}return false;};
 for(int degree=1;!found&&degree<=primitiveTargetFactorBound;degree++)found=extend(0,degree,0,one);
 cerr<<"DIXON_SEARCH_FINISHED TESTED "<<tested<<" ELIGIBLE "<<eligible<<" MODULAR_PASSES "<<modpassed<<" EXACT_SOLVES "<<exactsolves<<" FOUND "<<found<<endl;
 fmpq_clear(q);fmpq_mat_clear(M);fmpq_mat_clear(QB);fmpq_mat_clear(Qb);fmpq_mat_clear(Qx);nmod_mat_clear(N);nmod_mat_clear(E);nmod_mat_clear(MT);nmod_mat_clear(B);nmod_mat_clear(BI);nmod_mat_clear(b);nmod_mat_clear(x);return found;
}

void rd(istream&in,Poly&p){int nt;in>>nt;assert(in&&nt>=0);fmpq_t q;fmpq_init(q);ulong ex[3];for(int t=0;t<nt;t++){string c;in>>c>>ex[0]>>ex[1]>>ex[2];assert(in&&!fmpq_set_str(q,c.c_str(),10));fmpq_canonicalise(q);fmpq_mpoly_push_term_fmpq_ui(p.p,q,ex,ctx);}fmpq_mpoly_sort_terms(p.p,ctx);fmpq_mpoly_combine_like_terms(p.p,ctx);fmpq_clear(q);}
int main(int argc,char**argv){assert(argc==7);primitiveTargetFactorBound=stoi(argv[6]);fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);{ifstream in(argv[1]),uf(argv[2]);int n,m,k;in>>n>>m>>k;int metadata=stoi(argv[3]),target=stoi(argv[4]);assert(target>=0&&target<n);if(metadata){int discard;for(int z=0;z<k+n;z++)in>>discard;}else assert(metadata==0);const char*names[]={"a","e","d"};string text;while(getline(uf,text)){if(text.empty())continue;Poly q;assert(!fmpq_mpoly_set_str_pretty(q.p,text.c_str(),names,ctx));assert(!fmpq_mpoly_is_fmpq(q.p,ctx));known.push_back(q);}assert(!known.empty());vector<Poly>original(m);for(int i=0;i<n;i++)for(int j=0;j<m;j++){Poly q;rd(in,q);if(i==target)original[j]=q;}
vector<vector<int>>val(m,vector<int>(known.size()));vector<Poly>base=original,scales(m);Poly quotient;fmpq_t scalar;fmpq_init(scalar);
for(int j=0;j<m;j++){fmpq_mpoly_one(scales[j].p,ctx);if(fmpq_mpoly_is_zero(base[j].p,ctx))continue;for(int t=0;t<(int)known.size();t++)while(fmpq_mpoly_divides(quotient.p,base[j].p,known[t].p,ctx)){base[j]=quotient;val[j][t]++;}fmpq_mpoly_get_term_coeff_fmpq(scalar,base[j].p,0,ctx);fmpq_mpoly_set_fmpq(scales[j].p,scalar,ctx);fmpq_mpoly_scalar_div_fmpq(base[j].p,base[j].p,scalar,ctx);}
vector<Poly>monomials(3);fmpq_mpoly_one(monomials[0].p,ctx);assert(!fmpq_mpoly_set_str_pretty(monomials[1].p,"a",names,ctx));assert(!fmpq_mpoly_set_str_pretty(monomials[2].p,"e",names,ctx));
vector<vector<Poly>>primitive(1,vector<Poly>(3*m));for(int r=0;r<3;r++)for(int j=0;j<m;j++)fmpq_mpoly_mul(primitive[0][r*m+j].p,base[j].p,monomials[r].p,ctx);
cerr<<"EXACT_DEGREE1_COFACTOR_MATRIX ROW "<<target<<" ORIGINAL_COLUMNS "<<m<<" EXPANDED_COLUMNS "<<3*m<<endl;
vector<int>empty;int bi,bj;if(!combo(primitive,0,empty,bi,bj)){cerr<<"NO_LOCALIZED_CONSTANT_SPAN_WITNESS\n";return 2;}assert(primitiveWitness.size()==(size_t)(3*m)&&primitiveTarget.size()==1);vector<Poly>coefficients(m);Poly addend;
for(int r=0;r<3;r++)for(int j=0;j<m;j++){fmpq_mpoly_mul(addend.p,primitiveWitness[r*m+j].p,monomials[r].p,ctx);fmpq_mpoly_add(coefficients[j].p,coefficients[j].p,addend.p,ctx);}
for(int j=0;j<m;j++){fmpq_mpoly_get_fmpq(scalar,scales[j].p,ctx);fmpq_mpoly_scalar_div_fmpq(coefficients[j].p,coefficients[j].p,scalar,ctx);}fmpq_clear(scalar);
vector<int>maxval(known.size());for(int j=0;j<m;j++)if(!fmpq_mpoly_is_zero(coefficients[j].p,ctx))for(int t=0;t<(int)known.size();t++)maxval[t]=max(maxval[t],val[j][t]);
Poly common,tmp,sum,expected;fmpq_mpoly_one(common.p,ctx);for(int t=0;t<(int)known.size();t++)for(int h=0;h<maxval[t];h++){fmpq_mpoly_mul(tmp.p,common.p,known[t].p,ctx);common=tmp;}vector<Poly>lift(m);long terms=0,degree=0;int active=0;for(int j=0;j<m;j++){if(fmpq_mpoly_is_zero(coefficients[j].p,ctx))continue;lift[j]=coefficients[j];for(int t=0;t<(int)known.size();t++)for(int h=val[j][t];h<maxval[t];h++){fmpq_mpoly_mul(tmp.p,lift[j].p,known[t].p,ctx);lift[j]=tmp;}fmpq_mpoly_mul(tmp.p,original[j].p,lift[j].p,ctx);fmpq_mpoly_add(sum.p,sum.p,tmp.p,ctx);terms+=fmpq_mpoly_length(lift[j].p,ctx);degree=max(degree,fmpq_mpoly_total_degree_si(lift[j].p,ctx));active++;}
fmpq_mpoly_mul(expected.p,common.p,primitiveTarget[0].p,ctx);assert(fmpq_mpoly_equal(sum.p,expected.p,ctx)&&unit(expected));ofstream out(argv[5]);out<<m<<'\n';for(auto&p:lift)write(out,p);cerr<<"EXACT_KNOWN_UNIT_DENOMINATOR_LIFT_PASS ROW "<<target<<" ACTIVE "<<active<<" TERMS "<<terms<<" DEGREE "<<degree<<" COMMON_DEGREE "<<fmpq_mpoly_total_degree_si(common.p,ctx)<<"\n";primitiveWitness.clear();primitiveTarget.clear();known.clear();}fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
