// Finite coefficient-degree unit lift over QQ[x0,x1,x2].
// CLI: input max-degree output row-bounds. Plain input n*m polynomial entries;
// the last row is the negative target, with coefficient-degree bound zero.
// A fixed-prime echelon chooses a square exact system. Rational
// reconstruction or exact QQ Dixon solve supplies ONE candidate; success is
// reported only after every original QQ polynomial identity is checked.
// Exit4 is a refusal, not a QQ nonmembership assertion. No complete-kernel
// or all-degree assertion is made. The last output coefficient is exactly1.
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/nmod_mat.h>
#include <flint/ulong_extras.h>
#include <flint/fmpq_mpoly.h>
#include <map>
#include <array>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cassert>
#include <random>
#include <numeric>
#include <algorithm>
#include <deque>
#include <chrono>
using namespace std;fmpq_mpoly_ctx_t ctx;
struct P{fmpq_mpoly_t p;P(){fmpq_mpoly_init(p,ctx);}P(const P&x):P(){fmpq_mpoly_set(p,x.p,ctx);}P(P&&x):P(){fmpq_mpoly_swap(p,x.p,ctx);}P&operator=(P x){fmpq_mpoly_swap(p,x.p,ctx);return *this;}~P(){fmpq_mpoly_clear(p,ctx);}};
void rd(istream&in,P&p){int nt;in>>nt;fmpq_t q;fmpq_init(q);ulong ex[3];for(int t=0;t<nt;t++){string c;in>>c>>ex[0]>>ex[1]>>ex[2];assert(!fmpq_set_str(q,c.c_str(),10));fmpq_canonicalise(q);fmpq_mpoly_push_term_fmpq_ui(p.p,q,ex,ctx);}fmpq_mpoly_sort_terms(p.p,ctx);fmpq_mpoly_combine_like_terms(p.p,ctx);fmpq_clear(q);}
void wr(ostream&o,const P&p){int nt=fmpq_mpoly_length(p.p,ctx);o<<nt;fmpq_t q;fmpq_init(q);ulong ex[3];for(int t=0;t<nt;t++){fmpq_mpoly_get_term_exp_ui(ex,p.p,t,ctx);fmpq_mpoly_get_term_coeff_fmpq(q,p.p,t,ctx);char*c=fmpq_get_str(nullptr,10,q);o<<' '<<c;flint_free(c);for(auto x:ex)o<<' '<<x;}o<<'\n';fmpq_clear(q);}
int main(int argc,char**argv){assert(argc==5);fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);{ifstream in(argv[1]);int deg=stoi(argv[2]),n,m,meta;in>>n>>m>>meta;string line;getline(in,line);vector<vector<P>>A(n,vector<P>(m));for(auto&row:A)for(auto&p:row)rd(in,p);assert(in.good()||in.eof());vector<array<ulong,3>>mus;for(int total=0;total<=deg;total++)for(int a=0;a<=total;a++)for(int e=0;e<=total-a;e++)mus.push_back({(ulong)a,(ulong)e,(ulong)(total-a-e)});int nm=mus.size();vector<int>rowCount(n,nm),bounds(n,deg),offset(n+1);{ifstream bf(argv[4]);for(int i=0;i<n;i++){bf>>bounds[i];assert(bf&&bounds[i]>=0&&bounds[i]<=deg);rowCount[i]=0;for(auto mu:mus)if((int)(mu[0]+mu[1]+mu[2])<=bounds[i])rowCount[i]++;}}for(int i=0;i<n;i++)offset[i+1]=offset[i]+rowCount[i];assert(bounds[n-1]==0);int nu0=offset[n];vector<int>unknownRow(nu0),unknownMu(nu0);for(int i=0;i<n;i++)for(int u=0;u<rowCount[i];u++){unknownRow[offset[i]+u]=i;unknownMu[offset[i]+u]=u;}cerr<<"ROW_COEFFICIENT_DEGREE_BOUNDS ";for(int b:bounds)cerr<<b<<' ';cerr<<'\n';map<array<ulong,4>,int>eq;ulong ex[3];for(int j=0;j<m;j++)for(int i=0;i<n;i++)for(int t=0;t<fmpq_mpoly_length(A[i][j].p,ctx);t++){fmpq_mpoly_get_term_exp_ui(ex,A[i][j].p,t,ctx);for(int u=0;u<rowCount[i];u++){auto mu=mus[u];array<ulong,4>key={(ulong)j,ex[0]+mu[0],ex[1]+mu[1],ex[2]+mu[2]};if(!eq.count(key))eq[key]=eq.size();}}

const ulong prime=(UWORD(1)<<61)-1;assert(n_is_prime(prime));
auto mul=[&](ulong a,ulong b)->ulong{unsigned __int128 t=(unsigned __int128)a*b;ulong r=(ulong)(t&prime)+(ulong)(t>>61);return r>=prime?r-prime:r;};
auto sub=[&](ulong a,ulong b)->ulong{return a>=b?a-b:a+prime-b;};
using Row=vector<pair<int,ulong>>;vector<Row>equations(eq.size());vector<int>freq(nu0);fmpq_t q;fmpq_init(q);ulong ex2[3];
for(int i=0;i<n;i++)for(int j=0;j<m;j++)for(int t=0;t<fmpq_mpoly_length(A[i][j].p,ctx);t++){fmpq_mpoly_get_term_exp_ui(ex2,A[i][j].p,t,ctx);fmpq_mpoly_get_term_coeff_fmpq(q,A[i][j].p,t,ctx);ulong den=fmpz_fdiv_ui(fmpq_denref(q),prime);assert(den);ulong co=mul(fmpz_fdiv_ui(fmpq_numref(q),prime),n_invmod(den,prime));if(!co)continue;for(int u=0;u<rowCount[i];u++){auto mu=mus[u];int e=eq.at({(ulong)j,ex2[0]+mu[0],ex2[1]+mu[1],ex2[2]+mu[2]});int v=offset[i]+u;equations[e].push_back({v,co});freq[v]++;}}
vector<int>order(nu0),perm(nu0);iota(order.begin(),order.end(),0);sort(order.begin(),order.end(),[&](int i,int j){if(i==nu0-1||j==nu0-1)return j==nu0-1&&i!=nu0-1;return freq[i]!=freq[j]?freq[i]<freq[j]:i<j;});for(int i=0;i<nu0;i++)perm[order[i]]=i;for(auto&row:equations){for(auto&term:row)term.first=perm[term.first];sort(row.begin(),row.end());for(int i=1;i<(int)row.size();i++)assert(row[i-1].first!=row[i].first);}
vector<int>eqorder(equations.size());iota(eqorder.begin(),eqorder.end(),0);sort(eqorder.begin(),eqorder.end(),[&](int i,int j){return equations[i].size()!=equations[j].size()?equations[i].size()<equations[j].size():i<j;});
vector<Row>basis(nu0);vector<int>chosen;long basisTerms=0,operations=0;int processed=0;auto start=chrono::steady_clock::now();cerr<<"SPARSE_FINITE_MATRIX "<<eq.size()<<' '<<nu0<<" PRIME "<<prime<<'\n';
vector<vector<P>>certifiedS;bool earlyCertified=false;int lastTriedRank=-1;auto tryCertificate=[&]()->bool{vector<int>free;for(int i=0;i<nu0;i++)if(basis[i].empty())free.push_back(i);assert((int)free.size()==nu0-(int)chosen.size());vector<vector<P>>S;P tmp,sum;fmpz_t modulus,residue;fmpz_init_set_ui(modulus,prime);fmpz_init(residue);
int targetVariable=perm[nu0-1];if(find(free.begin(),free.end(),targetVariable)==free.end()){fmpz_clear(modulus);fmpz_clear(residue);return false;}for(int f:vector<int>{targetVariable}){vector<ulong>values(nu0);values[f]=1;for(int i=nu0-1;i>=0;i--)if(!basis[i].empty()){ulong sumv=0;for(int j=1;j<(int)basis[i].size();j++){auto term=basis[i][j];sumv=n_addmod(sumv,mul(term.second,values[term.first]),prime);}values[i]=sumv?prime-sumv:0;}vector<P>relation(n);for(int i=0;i<nu0;i++)if(values[i]){fmpz_set_ui(residue,values[i]);if(!fmpq_reconstruct_fmpz(q,residue,modulus)){fmpz_clear(modulus);fmpz_clear(residue);return false;}int orig=order[i];auto mu=mus[unknownMu[orig]];ulong x[3]={mu[0],mu[1],mu[2]};fmpq_mpoly_set_coeff_fmpq_ui(relation[unknownRow[orig]].p,q,x,ctx);}for(int j=0;j<m;j++){fmpq_mpoly_zero(sum.p,ctx);for(int i=0;i<n;i++){fmpq_mpoly_mul(tmp.p,relation[i].p,A[i][j].p,ctx);fmpq_mpoly_add(sum.p,sum.p,tmp.p,ctx);}if(!fmpq_mpoly_is_zero(sum.p,ctx)){fmpz_clear(modulus);fmpz_clear(residue);return false;}}S.push_back(move(relation));}
fmpz_clear(modulus);fmpz_clear(residue);certifiedS=move(S);return true;};
for(int e:eqorder){Row row=equations[e];while(!row.empty()){int pivot=row[0].first;if(basis[pivot].empty()){ulong inv=n_invmod(row[0].second,prime);for(auto&term:row)term.second=mul(term.second,inv);basisTerms+=row.size();assert(basisTerms<25000000);basis[pivot]=move(row);chosen.push_back(e);if(chosen.size()%500==0)cerr<<"SPARSE_PIVOTS "<<chosen.size()<<" TERMS "<<basisTerms<<" SEC "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';break;}auto&b=basis[pivot];ulong mult=row[0].second;Row out;out.reserve(row.size()+b.size());int x=1,y=1;while(x<(int)row.size()||y<(int)b.size()){if(y==(int)b.size()||(x<(int)row.size()&&row[x].first<b[y].first)){out.push_back(row[x++]);}else if(x==(int)row.size()||b[y].first<row[x].first){ulong c=mul(mult,b[y].second);if(c)out.push_back({b[y].first,prime-c});y++;}else{ulong c=sub(row[x].second,mul(mult,b[y].second));if(c)out.push_back({row[x].first,c});x++;y++;}operations++;}row=move(out);}processed++;if((int)chosen.size()==nu0)break;if((int)chosen.size()>=1&&chosen.size()%25==0&&lastTriedRank!=(int)chosen.size()){lastTriedRank=chosen.size();if(tryCertificate()){earlyCertified=true;cerr<<"EARLY_EXACT_QQ_KERNEL_CERTIFICATE RANK "<<chosen.size()<<" PROCESSED "<<processed<<'\n';break;}}if(processed%10000==0)cerr<<"SPARSE_PROCESSED "<<processed<<" RANK "<<chosen.size()<<" TERMS "<<basisTerms<<" OPERATIONS "<<operations<<'\n';}
int rank=chosen.size(),dim=nu0-rank;cerr<<"TARGET_FREE "<<basis[perm[nu0-1]].empty()<<"\n";cerr<<"SPARSE_COMPLETE_RANK "<<rank<<" NULLITY "<<dim<<" SEC "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';ofstream selected(string(argv[3])+".selected-equations");for(int e:chosen)selected<<e<<'\n';
if(!earlyCertified&&!tryCertificate()){
int targetVariable=perm[nu0-1];if(!basis[targetVariable].empty()){cerr<<"NO_CERTIFICATE_AT_THIS_PRIME_AND_BOUND\n";return 4;}
vector<int>rowIndex(eq.size(),-1),columnIndex(nu0,-1),pivotVariables;for(int i=0;i<rank;i++)rowIndex[chosen[i]]=i;for(int v=0;v<nu0;v++)if(!basis[v].empty()){columnIndex[v]=pivotVariables.size();pivotVariables.push_back(v);}assert((int)pivotVariables.size()==rank);
fmpq_mat_t BB,bb,xx;fmpq_mat_init(BB,rank,rank);fmpq_mat_init(bb,rank,1);fmpq_mat_init(xx,rank,1);
for(int i=0;i<n;i++)for(int j=0;j<m;j++)for(int t=0;t<fmpq_mpoly_length(A[i][j].p,ctx);t++){fmpq_mpoly_get_term_exp_ui(ex2,A[i][j].p,t,ctx);fmpq_mpoly_get_term_coeff_fmpq(q,A[i][j].p,t,ctx);for(int u=0;u<rowCount[i];u++){auto mu=mus[u];int er=eq.at({(ulong)j,ex2[0]+mu[0],ex2[1]+mu[1],ex2[2]+mu[2]});int rr=rowIndex[er];if(rr<0)continue;int v=perm[offset[i]+u];if(v==targetVariable)fmpq_sub(fmpq_mat_entry(bb,rr,0),fmpq_mat_entry(bb,rr,0),q);else if(columnIndex[v]>=0)fmpq_add(fmpq_mat_entry(BB,rr,columnIndex[v]),fmpq_mat_entry(BB,rr,columnIndex[v]),q);}}
cerr<<"EXACT_QQ_DIXON_SOLVE_BEGIN "<<rank<<"\n";assert(fmpq_mat_solve_dixon(xx,BB,bb));cerr<<"EXACT_QQ_DIXON_SOLVE_DONE\n";
vector<P>relation(n);fmpq_mpoly_one(relation[n-1].p,ctx);for(int k=0;k<rank;k++)if(!fmpq_is_zero(fmpq_mat_entry(xx,k,0))){int orig=order[pivotVariables[k]];auto mu=mus[unknownMu[orig]];ulong x[3]={mu[0],mu[1],mu[2]};fmpq_mpoly_set_coeff_fmpq_ui(relation[unknownRow[orig]].p,fmpq_mat_entry(xx,k,0),x,ctx);}
P sum,tmp;for(int j=0;j<m;j++){fmpq_mpoly_zero(sum.p,ctx);for(int i=0;i<n;i++){fmpq_mpoly_mul(tmp.p,relation[i].p,A[i][j].p,ctx);fmpq_mpoly_add(sum.p,sum.p,tmp.p,ctx);}if(!fmpq_mpoly_is_zero(sum.p,ctx)){cerr<<"EXACT_QQ_IDENTITY_REFUSAL\n";return 4;}}
certifiedS.push_back(move(relation));fmpq_mat_clear(BB);fmpq_mat_clear(bb);fmpq_mat_clear(xx);
}assert(certifiedS.size()==1&&fmpq_mpoly_is_one(certifiedS[0][n-1].p,ctx));auto&S=certifiedS;
ofstream out(argv[3]);out<<S.size()<<' '<<n<<" 0\n\n";for(int i=0;i<(int)S.size();i++)out<<i<<' ';out<<'\n';for(auto&row:S)for(auto&p:row)wr(out,p);cerr<<"EXACT_FULL_QQ_UNIT_LIFT_PASS "<<S.size()<<" DEGREE "<<deg<<'\n';fmpq_clear(q);}
fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
