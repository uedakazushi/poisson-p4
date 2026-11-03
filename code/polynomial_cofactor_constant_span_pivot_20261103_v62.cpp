// Optional eighth CLI argument: rational constant-combination file (column count, then coefficients).
// It exposes a constant in the cofactor span by an invertible QQ column change.
// The complete polynomial image is independently recomputed before/after that change.
// One exact known-unit witness pivot; no greedy continuation.
// CLI: input output known-units metadata(0/1) witness target-row.
// Checks the full polynomial image of the witness and a known-unit target.
// Replacement uses a known-unit witness coefficient; otherwise appends a redundant column.
// Exact polynomial unit pivots, including constant column/row span search.
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


vector<Poly> known;int forcedWitnessColumn=-1;
bool unit(Poly p){if(fmpq_mpoly_is_zero(p.p,ctx))return false;Poly q;for(auto&f:known)while(fmpq_mpoly_divides(q.p,p.p,f.p,ctx))p=q;return fmpq_mpoly_is_fmpq(p.p,ctx);}
void clean(vector<Poly>&row){int first=-1;for(int j=0;j<(int)row.size();j++)if(!fmpq_mpoly_is_zero(row[j].p,ctx)){first=j;break;}if(first<0)return;fmpq_t c;fmpq_init(c);fmpq_mpoly_get_term_coeff_fmpq(c,row[first].p,0,ctx);for(auto&p:row)fmpq_mpoly_scalar_div_fmpq(p.p,p.p,c,ctx);fmpq_clear(c);vector<Poly> quot(row.size());for(auto&f:known)while(true){bool yes=true;for(int j=first;j<(int)row.size();j++)if(!fmpq_mpoly_divides(quot[j].p,row[j].p,f.p,ctx)){yes=false;break;}if(!yes)break;row=quot;}}

bool combo(vector<vector<Poly>>&A,int k,const vector<int>&used,int &bi,int &bj){
 int n=A.size(),m=A[0].size();vector<int> cols;for(int j=0;j<m;j++)if(find(used.begin(),used.end(),j)==used.end())cols.push_back(j);
 vector<Poly>targets;Poly one;fmpq_mpoly_one(one.p,ctx);targets.push_back(one);

 function<void(int,int,const Poly&)>extend=[&](int pos,int left,const Poly&prod){
  if(pos==(int)known.size()-1){Poly v=prod,tmp;for(int e=0;e<left;e++){fmpq_mpoly_mul(tmp.p,v.p,known[pos].p,ctx);v=tmp;}targets.push_back(v);return;}
  Poly v=prod,tmp;for(int e=0;e<=left;e++){extend(pos+1,left-e,v);fmpq_mpoly_mul(tmp.p,v.p,known[pos].p,ctx);v=tmp;}
 };
 for(int deg=1;deg<=3;deg++)extend(0,deg,one);

 for(int i=k;i<n;i++){
  map<array<ulong,3>,int>mons;ulong ex[3];
  auto collect=[&](const Poly&p){for(int t=0;t<fmpq_mpoly_length(p.p,ctx);t++){fmpq_mpoly_get_term_exp_ui(ex,p.p,t,ctx);array<ulong,3>key={ex[0],ex[1],ex[2]};if(!mons.count(key))mons[key]=mons.size();}};
  for(int j:cols)collect(A[i][j]);for(auto&t:targets)collect(t);
  int nc=cols.size(),nt=targets.size(),nr=mons.size();fmpq_mat_t M,E;fmpq_mat_init(M,nr,nc+nt);fmpq_mat_init(E,nr,nc+nt);fmpq_t q;fmpq_init(q);
  auto put=[&](const Poly&p,int c){for(int t=0;t<fmpq_mpoly_length(p.p,ctx);t++){fmpq_mpoly_get_term_exp_ui(ex,p.p,t,ctx);fmpq_mpoly_get_term_coeff_fmpq(q,p.p,t,ctx);array<ulong,3>key={ex[0],ex[1],ex[2]};fmpq_set(fmpq_mat_entry(M,mons.at(key),c),q);}};
  for(int j=0;j<nc;j++)put(A[i][cols[j]],j);for(int j=0;j<nt;j++)put(targets[j],nc+j);
  fmpq_mat_rref(E,M);vector<int> pivots;
  for(int r=0;r<nr;r++){int c=0;for(;c<nc;c++)if(!fmpq_is_zero(fmpq_mat_entry(E,r,c)))break;if(c==nc)break;pivots.push_back(c);}
  int rr=pivots.size();int found=-1;
  for(int j=0;j<nt;j++){bool yes=true;for(int r=rr;r<nr;r++)if(!fmpq_is_zero(fmpq_mat_entry(E,r,nc+j))){yes=false;break;}if(yes){found=j;break;}}
  if(found>=0){int pivot=-1;for(int r=0;r<rr;r++)if(!fmpq_is_zero(fmpq_mat_entry(E,r,nc+found))){pivot=r;break;}assert(pivot>=0);bj=cols[pivots[pivot]];bi=i;
   vector<Poly>newcol(n);Poly tmp;for(int r=0;r<rr;r++){auto c=fmpq_mat_entry(E,r,nc+found);if(fmpq_is_zero(c))continue;for(int l=0;l<n;l++){fmpq_mpoly_scalar_mul_fmpq(tmp.p,A[l][cols[pivots[r]]].p,c,ctx);fmpq_mpoly_add(newcol[l].p,newcol[l].p,tmp.p,ctx);}}
   assert(fmpq_mpoly_equal(newcol[i].p,targets[found].p,ctx));for(int l=0;l<n;l++)A[l][bj]=newcol[l];cerr<<"COMBO row "<<i<<" rank "<<rr<<" target "<<found<<" col "<<bj<<"\n";
   fmpq_clear(q);fmpq_mat_clear(M);fmpq_mat_clear(E);return true;
  }
  cerr<<"NO_COMBO row "<<i<<" rank "<<rr<<" mons "<<nr<<"\n";fmpq_clear(q);fmpq_mat_clear(M);fmpq_mat_clear(E);
 }
 return false;
}

bool combo_rows(vector<vector<Poly>>&A,int k,const vector<int>&used,int &bi,int &bj){
 int n=A.size(),m=A[0].size();vector<int>cols;for(int j=0;j<m;j++)if(find(used.begin(),used.end(),j)==used.end())cols.push_back(j);
 vector<vector<Poly>>T(cols.size(),vector<Poly>(n-k));for(int i=k;i<n;i++)for(int j=0;j<(int)cols.size();j++)T[j][i-k]=A[i][cols[j]];
 int ti,tj;vector<int>none;bool yes=combo(T,0,none,ti,tj);if(!yes)return false;
 for(int i=k;i<n;i++)for(int j=0;j<(int)cols.size();j++)A[i][cols[j]]=T[j][i-k];bi=k+tj;bj=cols[ti];return true;
}

string coefficientComboPath;
void apply_witness(vector<vector<Poly>>&A,const char*filename,int target){
 int n=A.size(),m=A[0].size();assert(target>=0&&target<n);ifstream in(filename);int count;in>>count;assert(count==m);vector<Poly>c(m);fmpq_t q;fmpq_init(q);ulong ex[3];
 for(auto&p:c){int terms;in>>terms;assert(terms>=0);for(int t=0;t<terms;t++){string text;in>>text>>ex[0]>>ex[1]>>ex[2];assert(!fmpq_set_str(q,text.c_str(),10));fmpq_canonicalise(q);fmpq_mpoly_push_term_fmpq_ui(p.p,q,ex,ctx);}fmpq_mpoly_sort_terms(p.p,ctx);fmpq_mpoly_combine_like_terms(p.p,ctx);}

 vector<Poly>oldImage;
 if(!coefficientComboPath.empty()){
  ifstream combo(coefficientComboPath);int nc;combo>>nc;assert(nc==m);vector<Poly>alpha(m);int replace=-1,active=0;Poly exposed,tmp;
  for(int j=0;j<m;j++){string text;combo>>text;assert(!fmpq_set_str(q,text.c_str(),10));fmpq_canonicalise(q);fmpq_mpoly_set_fmpq(alpha[j].p,q,ctx);if(!fmpq_is_zero(q)){active++;if(replace<0)replace=j;fmpq_mpoly_mul(tmp.p,alpha[j].p,c[j].p,ctx);fmpq_mpoly_add(exposed.p,exposed.p,tmp.p,ctx);}}
  assert(replace>=0&&!fmpq_mpoly_is_zero(exposed.p,ctx)&&fmpq_mpoly_is_fmpq(exposed.p,ctx));
  oldImage=vector<Poly>(n);for(int i=0;i<n;i++)for(int j=0;j<m;j++){fmpq_mpoly_mul(tmp.p,A[i][j].p,c[j].p,ctx);fmpq_mpoly_add(oldImage[i].p,oldImage[i].p,tmp.p,ctx);}
  fmpq_t ai,aj,ratio;fmpq_init(ai);fmpq_init(aj);fmpq_init(ratio);fmpq_mpoly_get_fmpq(ai,alpha[replace].p,ctx);
  for(int i=0;i<n;i++){Poly old=A[i][replace],scaled;fmpq_mpoly_scalar_div_fmpq(A[i][replace].p,old.p,ai,ctx);for(int j=0;j<m;j++)if(j!=replace&&!fmpq_mpoly_is_zero(alpha[j].p,ctx)){fmpq_mpoly_get_fmpq(aj,alpha[j].p,ctx);fmpq_div(ratio,aj,ai);fmpq_mpoly_scalar_mul_fmpq(scaled.p,old.p,ratio,ctx);fmpq_mpoly_sub(A[i][j].p,A[i][j].p,scaled.p,ctx);}}
  c[replace]=exposed;fmpq_clear(ai);fmpq_clear(aj);fmpq_clear(ratio);cerr<<"EXACT_COFACTOR_CONSTANT_COMBINATION "<<active<<" EXPOSED_COLUMN "<<replace<<"\n";
 }
 vector<Poly>v(n);Poly tmp;for(int j=0;j<m;j++)if(!fmpq_mpoly_is_zero(c[j].p,ctx))for(int i=0;i<n;i++){fmpq_mpoly_mul(tmp.p,A[i][j].p,c[j].p,ctx);fmpq_mpoly_add(v[i].p,v[i].p,tmp.p,ctx);}
 if(!oldImage.empty())for(int i=0;i<n;i++)assert(fmpq_mpoly_equal(oldImage[i].p,v[i].p,ctx));
 for(int i=0;i<n;i++)cerr<<"IMAGE_ROW "<<i<<" DEGREE "<<fmpq_mpoly_total_degree_si(v[i].p,ctx)<<" TERMS "<<fmpq_mpoly_length(v[i].p,ctx)<<"\n";assert(unit(v[target]));int replace=-1;for(int j=0;j<m;j++)if(unit(c[j])){replace=j;break;}
 forcedWitnessColumn=replace>=0?replace:m;if(replace>=0){for(int i=0;i<n;i++)A[i][replace]=v[i];cerr<<"EXACT_WITNESS_REPLACE "<<replace<<" TARGET "<<target<<"\n";}
 else{for(int i=0;i<n;i++)A[i].push_back(v[i]);cerr<<"EXACT_WITNESS_APPEND "<<m<<" TARGET "<<target<<"\n";}
 fmpq_clear(q);
}

void clean_columns(vector<vector<Poly>>&A,int k,const vector<int>&used){
 int n=A.size(),m=A[0].size();for(int j=0;j<m;j++){if(find(used.begin(),used.end(),j)!=used.end())continue;vector<Poly>col;for(int i=k;i<n;i++)col.push_back(A[i][j]);clean(col);for(int i=k;i<n;i++)A[i][j]=col[i-k];}
}
int main(int argc,char**argv){if(argc==8){coefficientComboPath=argv[7];argc--;}assert(argc>=5);fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);{const char*names[]={"a","e","d"};ifstream unitfile(argv[3]);string unittext;while(getline(unitfile,unittext)){if(unittext.empty())continue;Poly p;assert(!fmpq_mpoly_set_str_pretty(p.p,unittext.c_str(),names,ctx));assert(!fmpq_mpoly_is_fmpq(p.p,ctx));known.push_back(p);}assert(!known.empty());ifstream in(argv[1]);ofstream out(argv[2]);int n,m,r;in>>n>>m>>r;if(string(argv[4])=="1"){int discard;for(int l=0;l<r+n;l++)in>>discard;}else assert(string(argv[4])=="0");vector<vector<Poly>>A(n,vector<Poly>(m));fmpq_t q;fmpq_init(q);ulong ex[3];for(auto&row:A)for(auto&p:row){int l;in>>l;for(int k=0;k<l;k++){string t;in>>t>>ex[0]>>ex[1]>>ex[2];assert(!fmpq_set_str(q,t.c_str(),10));fmpq_canonicalise(q);fmpq_mpoly_push_term_fmpq_ui(p.p,q,ex,ctx);}fmpq_mpoly_sort_terms(p.p,ctx);fmpq_mpoly_combine_like_terms(p.p,ctx);}bool do_column_clean=false;int argindex=5;if(argc>5&&string(argv[5])=="transpose"){argindex++;vector<vector<Poly>>T(m,vector<Poly>(n));for(int i=0;i<n;i++)for(int j=0;j<m;j++)T[j][i]=A[i][j];A=move(T);swap(n,m);}if(argc>argindex&&string(argv[argindex])=="colclean"){do_column_clean=true;argindex++;}for(auto&row:A)clean(row);if(argc>argindex){assert(argc==argindex+2);apply_witness(A,argv[argindex],stoi(argv[argindex+1]));m=A[0].size();}vector<int>used;if(do_column_clean){clean_columns(A,0,used);for(auto&row:A)clean(row);}vector<int>origin(n);for(int i=0;i<n;i++)origin[i]=i;auto start=chrono::steady_clock::now();int k=0;assert(argc==7&&forcedWitnessColumn>=0);for(;k<1;k++){int bi=stoi(argv[6]),bj=forcedWitnessColumn;assert(bi>=0&&bi<n&&unit(A[bi][bj]));swap(A[k],A[bi]);swap(origin[k],origin[bi]);used.push_back(bj);Poly pivot=A[k][bj],x,y;for(int i=k+1;i<n;i++){if(fmpq_mpoly_is_zero(A[i][bj].p,ctx))continue;Poly mult=A[i][bj];for(int j=0;j<m;j++){fmpq_mpoly_mul(x.p,pivot.p,A[i][j].p,ctx);fmpq_mpoly_mul(y.p,mult.p,A[k][j].p,ctx);fmpq_mpoly_sub(A[i][j].p,x.p,y.p,ctx);}clean(A[i]);assert(fmpq_mpoly_is_zero(A[i][bj].p,ctx));}if(do_column_clean){clean_columns(A,k+1,used);for(int i=k+1;i<n;i++)clean(A[i]);}long terms=0;for(int i=k+1;i<n;i++)for(auto&p:A[i])terms+=fmpq_mpoly_length(p.p,ctx);cerr<<"PIVOT "<<k<<" COL "<<bj<<" ROW "<<origin[k]<<" TERMS "<<terms<<" SEC "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';assert(chrono::duration<double>(chrono::steady_clock::now()-start).count()<180 && terms<3000000);}out<<n-k<<' '<<m-k<<' '<<k<<'\n';for(int j:used)out<<j<<' ';out<<'\n';for(int i=k;i<n;i++)out<<origin[i]<<' ';out<<'\n';for(int i=k;i<n;i++)for(int j=0;j<m;j++)if(find(used.begin(),used.end(),j)==used.end())write(out,A[i][j]);cerr<<"PASS "<<k<<'\n';fmpq_clear(q);known.clear();}fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
