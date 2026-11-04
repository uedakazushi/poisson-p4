// ONE full67 point-free coefficient test; ONE prime, conditional exact QQ lift.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpz_mpoly.h>
using namespace std;
fmpz_mpoly_ctx_t ctx;constexpr ulong PRIME=32003;constexpr int NC=67,NG=10;long maxterms=0,muls=0;
struct P{fmpz_mpoly_t v;P(){fmpz_mpoly_init(v,ctx);}P(const P&p):P(){fmpz_mpoly_set(v,p.v,ctx);}P(P&&p):P(){fmpz_mpoly_swap(v,p.v,ctx);}P&operator=(P p){fmpz_mpoly_swap(v,p.v,ctx);return*this;}~P(){fmpz_mpoly_clear(v,ctx);}};
struct Z{fmpz_t v;Z(){fmpz_init(v);}Z(const Z&z):Z(){fmpz_set(v,z.v);}~Z(){fmpz_clear(v);}};
void guard(const P&p){long n=fmpz_mpoly_length(p.v,ctx);maxterms=max(maxterms,n);if(n>150000)throw runtime_error("150k intermediate term cap");}
P mul(const P&a,const P&b){P c;fmpz_mpoly_mul(c.v,a.v,b.v,ctx);muls++;guard(c);return c;}
P one(){P p;fmpz_mpoly_one(p.v,ctx);return p;}
void readpoly(istream&in,P&p){long n;in>>n;assert(n>=0&&n<=100000);Z z;ulong e[4];for(long k=0;k<n;k++){string c;in>>c;for(auto&i:e)in>>i;assert(in&&fmpz_set_str(z.v,c.c_str(),10)==0);fmpz_mpoly_push_term_fmpz_ui(p.v,z.v,e,ctx);}fmpz_mpoly_sort_terms(p.v,ctx);fmpz_mpoly_combine_like_terms(p.v,ctx);guard(p);}
string str(const fmpz_t z){char*p=fmpz_get_str(nullptr,10,z);string s(p);flint_free(p);return s;}
ulong invmod(ulong x){assert(x);ulong ans=1,p=PRIME-2;while(p){if(p&1)ans=ans*x%PRIME;x=x*x%PRIME;p>>=1;}return ans;}
ulong detmod(vector<vector<ulong>>A){ulong det=1;int n=A.size();for(int k=0;k<n;k++){int p=k;while(p<n&&!A[p][k])p++;assert(p<n);if(p!=k){swap(A[p],A[k]);det=(PRIME-det)%PRIME;}det=det*A[k][k]%PRIME;ulong inv=invmod(A[k][k]);for(int i=k+1;i<n;i++){ulong q=A[i][k]*inv%PRIME;for(int j=k;j<n;j++)A[i][j]=(A[i][j]+PRIME-q*A[k][j]%PRIME)%PRIME;}}return det;}
void frontier(const string&path,const string&stage,const vector<long>&counts,int rank=-1){ofstream f(path);f<<"{\"stage\":\""<<stage<<"\",\"completed_columns\":"<<counts.size()<<",\"column_term_counts\":[";for(int i=0;i<int(counts.size());i++){if(i)f<<',';f<<counts[i];}f<<"],\"max_intermediate_terms\":"<<maxterms<<",\"polynomial_multiplications\":"<<muls<<",\"coefficient_rank_modp_if_finished\":"<<rank<<"}\n";}
int main(int argc,char**argv){assert(argc==4);fmpz_mpoly_ctx_init(ctx,4,ORD_DEGLEX);int exitcode=0;auto started=chrono::steady_clock::now();vector<long>counts;{
 try{
  ifstream in(argv[1]);int vars,ng,nc;ulong prime;in>>vars>>ng>>nc>>prime;assert(vars==4&&ng==NG&&nc==NC&&prime==PRIME);
  vector<P>G(NG);vector<Z>den(NG);vector<int>bp(NG),dp(NG);for(int i=0;i<NG;i++){string d;in>>d>>bp[i]>>dp[i];assert(fmpz_set_str(den[i].v,d.c_str(),10)==0&&fmpz_fdiv_ui(den[i].v,PRIME));readpoly(in,G[i]);}
  vector<array<int,NG>>exps(NC);for(auto&e:exps)for(auto&x:e)in>>x;assert(in);string extra;assert(!(in>>extra));
  vector<vector<P>>powers(NG);for(int i=0;i<NG;i++){int n=0;for(auto&e:exps)n=max(n,e[i]);powers[i].push_back(one());for(int k=1;k<=n;k++)powers[i].push_back(mul(powers[i].back(),G[i]));}
  vector<P>cols;vector<Z>colDen(NC);Z q;
  frontier(argv[3],"generator_powers_complete",counts);
  for(int j=0;j<NC;j++){P p=one();int bsum=0,dsum=0;fmpz_one(colDen[j].v);
   for(int i=0;i<NG;i++)if(exps[j][i]){p=mul(p,powers[i][exps[j][i]]);bsum+=bp[i]*exps[j][i];dsum+=dp[i]*exps[j][i];fmpz_pow_ui(q.v,den[i].v,exps[j][i]);fmpz_mul(colDen[j].v,colDen[j].v,q.v);}
   assert(dsum==24&&bsum<=12);if(bsum<12){P B,powB;fmpz_mpoly_gen(B.v,3,ctx);assert(fmpz_mpoly_pow_ui(powB.v,B.v,12-bsum,ctx));p=mul(p,powB);}guard(p);long terms=fmpz_mpoly_length(p.v,ctx);if(terms>100000)throw runtime_error("100k column term cap");counts.push_back(terms);cols.push_back(move(p));frontier(argv[3],"integer_columns_in_progress",counts);cerr<<"COLUMN "<<j<<" TERMS "<<terms<<endl;
  }
  using Key=array<ulong,4>;using Row=array<ulong,NC>;map<Key,Row>coeff;Z z;Key key;
  for(int j=0;j<NC;j++)for(long k=0;k<counts[j];k++){fmpz_mpoly_get_term_exp_ui(key.data(),cols[j].v,k,ctx);fmpz_mpoly_get_term_coeff_fmpz(z.v,cols[j].v,k,ctx);auto[it,added]=coeff.try_emplace(key);if(added)it->second.fill(0);it->second[j]=fmpz_fdiv_ui(z.v,PRIME);}
  array<Row,NC>basis;array<bool,NC>used{};vector<Key>selected;int rank=0;
  for(const auto&kv:coeff){Row row=kv.second;for(int p=0;p<NC;p++)if(used[p]&&row[p]){ulong f=row[p];for(int j=p;j<NC;j++)row[j]=(row[j]+PRIME-f*basis[p][j]%PRIME)%PRIME;}int p=0;while(p<NC&&!row[p])p++;if(p==NC)continue;ulong inv=invmod(row[p]);for(auto&c:row)c=c*inv%PRIME;basis[p]=row;used[p]=true;selected.push_back(kv.first);rank++;if(rank==NC)break;}
  vector<int>pivotcols;for(int p=0;p<NC;p++)if(used[p])pivotcols.push_back(p);vector<vector<ulong>>minor(rank,vector<ulong>(rank));for(int i=0;i<rank;i++)for(int j=0;j<rank;j++)minor[i][j]=coeff.at(selected[i])[pivotcols[j]];ulong det=detmod(minor);assert(det);frontier(argv[3],"one_GF_coefficient_rank_complete",counts,rank);
  ostringstream out;out<<"{\"all_passed\":true,\"prime\":32003,\"coefficient_rank_modp\":"<<rank<<",\"coefficient_union_rows\":"<<coeff.size()<<",\"selected_minor_modp\":"<<det<<",\"selected_rows\":[";
  for(int i=0;i<rank;i++){if(i)out<<',';out<<'[';for(int k=0;k<4;k++){if(k)out<<',';out<<selected[i][k];}out<<']';}out<<"],\"selected_columns\":[";for(int i=0;i<rank;i++){if(i)out<<',';out<<pivotcols[i];}out<<"],\"column_term_counts\":[";for(int i=0;i<NC;i++){if(i)out<<',';out<<counts[i];}out<<"],\"selected_integer_matrix\":[";
  fmpq_mat_t A,R;fmpq_mat_init(A,rank,NC);fmpq_mat_init(R,rank,NC);
  for(int i=0;i<rank;i++){if(i)out<<',';out<<'[';for(int j=0;j<NC;j++){if(j)out<<',';fmpz_mpoly_get_coeff_fmpz_ui(z.v,cols[j].v,selected[i].data(),ctx);fmpq_set_fmpz(fmpq_mat_entry(A,i,j),z.v);out<<'"'<<str(z.v)<<'"';}out<<']';}out<<']';
  bool identity=false,augmented=false;
  if(rank==NC-1){assert(fmpq_mat_rref(R,A)==NC-1);vector<int>pc;for(int i=0;i<rank;i++){int j=0;while(j<NC&&fmpq_is_zero(fmpq_mat_entry(R,i,j)))j++;assert(j<NC);pc.push_back(j);}int freecol=0;while(find(pc.begin(),pc.end(),freecol)!=pc.end())freecol++;assert(freecol<NC);
   vector<Z>kernel(NC),universal(NC);Z lcm,gcd,tmp;fmpz_one(lcm.v);for(int i=0;i<rank;i++)fmpz_lcm(lcm.v,lcm.v,fmpq_denref(fmpq_mat_entry(R,i,freecol)));fmpz_set(kernel[freecol].v,lcm.v);
   for(int i=0;i<rank;i++){fmpz_divexact(tmp.v,lcm.v,fmpq_denref(fmpq_mat_entry(R,i,freecol)));fmpz_mul(kernel[pc[i]].v,tmp.v,fmpq_numref(fmpq_mat_entry(R,i,freecol)));fmpz_neg(kernel[pc[i]].v,kernel[pc[i]].v);}
   P sum,term;for(int j=0;j<NC;j++)if(!fmpz_is_zero(kernel[j].v)){fmpz_mpoly_scalar_mul_fmpz(term.v,cols[j].v,kernel[j].v,ctx);fmpz_mpoly_add(sum.v,sum.v,term.v,ctx);guard(sum);}identity=fmpz_mpoly_is_zero(sum.v,ctx);
   fmpz_zero(gcd.v);for(int j=0;j<NC;j++){fmpz_mul(universal[j].v,kernel[j].v,colDen[j].v);fmpz_gcd(gcd.v,gcd.v,universal[j].v);}assert(!fmpz_is_zero(gcd.v));for(auto&c:universal)fmpz_divexact(c.v,c.v,gcd.v);int first=0;while(first<NC&&fmpz_is_zero(universal[first].v))first++;assert(first<NC);if(fmpz_sgn(universal[first].v)<0)for(auto&c:universal)fmpz_neg(c.v,c.v);
   out<<",\"exact_QQ_selected_rank\":66,\"exact_integer_coefficient_identity_verified\":"<<(identity?"true":"false")<<",\"nonzero_remainder_terms\":"<<fmpz_mpoly_length(sum.v,ctx)<<",\"integer_column_kernel\":[";for(int j=0;j<NC;j++){if(j)out<<',';out<<'"'<<str(kernel[j].v)<<'"';}out<<"],\"universal_primitive_integer_relation\":[";for(int j=0;j<NC;j++){if(j)out<<',';out<<'"'<<str(universal[j].v)<<'"';}out<<']';
   if(!identity){Key witness;fmpz_mpoly_get_term_exp_ui(witness.data(),sum.v,0,ctx);fmpz_mpoly_get_term_coeff_fmpz(z.v,sum.v,0,ctx);out<<",\"countermonomial\":[";for(int k=0;k<4;k++){if(k)out<<',';out<<witness[k];}out<<"],\"countercoefficient\":\""<<str(z.v)<<"\",\"counterrow_integer_coefficients\":[";Z dot;fmpz_zero(dot.v);for(int j=0;j<NC;j++){if(j)out<<',';fmpz_mpoly_get_coeff_fmpz_ui(z.v,cols[j].v,witness.data(),ctx);out<<'"'<<str(z.v)<<'"';fmpz_addmul(dot.v,z.v,kernel[j].v);}assert(!fmpz_is_zero(dot.v));out<<"],\"augmented_QQ_rank67_by_kernel_pairing\":true";augmented=true;}
  }
  string status=rank==NC?"all_Inv24_QQ_restrictions_injective":rank==NC-1?(identity?"verified_NEW25_degree24_relation":augmented?"all_Inv24_QQ_injective_by_counterrow":"QQ_candidate_failed"):"GF_rank_below66_inconclusive";
  out<<",\"status\":\""<<status<<"\",\"max_intermediate_terms\":"<<maxterms<<",\"polynomial_multiplications\":"<<muls<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"}\n";
  string text=out.str();if(text.size()>32L*1024*1024)throw runtime_error("32MiB certificate cap");ofstream file(argv[2]);file<<text;file.close();frontier(argv[3],status,counts,rank);fmpq_mat_clear(A);fmpq_mat_clear(R);cerr<<"COEFFICIENT_RANK "<<rank<<" IDENTITY "<<identity<<" STATUS "<<status<<endl;
 }catch(const exception&e){cerr<<"BOUNDED_TEST_STOP "<<e.what()<<endl;frontier(argv[3],"bounded_failure",counts);exitcode=2;}
 }fmpz_mpoly_ctx_clear(ctx);flint_cleanup();return exitcode;
}
