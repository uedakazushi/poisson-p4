// Exact quotient-algebra certificate for the zero12 first base-adic quotient.
// args: modular-input, QQ-coefficients, filtered-trace, selected-targets,
//       selected-output, observable-input, observable-matrix-output,
//       rational-queries, normal-form-answers.
// The guide trace is replayed over QQ. It proves only row membership until the
// quotient multiplication, monomial evaluation and original equations below
// have ALL been verified. Base coordinates are fixed to [0,1,3].
#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>
#include <map>
#include <gmpxx.h>
using namespace std;using Q=mpq_class;using Row=vector<pair<int,Q>>;
struct Trace{int pivot,kind,source,var;vector<int>red;};
int main(int argc,char**argv){if(argc!=10)return 1;ifstream in(argv[1]);int p,n,N,F;in>>p>>n>>N>>F;vector<int>deg(N);vector<vector<int>>shift(n,vector<int>(N));for(int&i:deg)in>>i;for(auto&r:shift)for(int&i:r)in>>i;int offset=0;while(deg[offset]<2)offset++;ifstream cf(argv[2]);vector<Row>init(F);for(int j=0;j<F;j++)for(int k=offset;k<N;k++){Q a;cf>>a;if(a!=0)init[j].emplace_back(k,std::move(a));}if(!cf){cerr<<"bad coefficients\n";return 2;}ifstream tr(argv[3]);vector<Trace>trace;vector<int>bid(N,-1);int id;while(tr>>id){Trace t;int c;tr>>t.pivot>>t.kind>>t.source>>t.var>>c;t.red.resize(c);for(int&i:t.red)tr>>i;bid[t.pivot]=id;trace.push_back(std::move(t));}ifstream ns(argv[4]);int nn;ns>>nn;vector<int>need(nn);for(int&i:need)ns>>i;int nt;ns>>nt;vector<int>targets(nt);for(int&i:targets)ns>>i;vector<Row>B(trace.size());vector<int>state(trace.size());int done=0;long long reductions=0;
function<void(int)> replay=[&](int id){if(state[id]==2)return;if(state[id]){cerr<<"cycle "<<id<<"\n";exit(5);}state[id]=1;auto&t=trace[id];vector<Q>r(N);if(t.kind==0)for(auto&[k,c]:init[t.source])r[k]=c;else{replay(t.source);for(auto&[k,c]:B[t.source])if(shift[t.var][k]>=0)r[shift[t.var][k]]=c;}Q product;for(int k=0;k<N;k++)if(r[k]!=0&&bid[k]>=0&&bid[k]<id){int j=bid[k];replay(j);Q a=r[k];for(auto&[l,c]:B[j]){mpq_mul(product.get_mpq_t(),a.get_mpq_t(),c.get_mpq_t());mpq_sub(r[l].get_mpq_t(),r[l].get_mpq_t(),product.get_mpq_t());}reductions++;}int first=0;while(first<N&&r[first]==0)first++;if(first!=t.pivot){cerr<<"pivot mismatch "<<id<<" expected "<<t.pivot<<" actual "<<first<<"\n";exit(3);}Q a=r[t.pivot];for(int k=first;k<N;k++)if(r[k]!=0)B[id].emplace_back(k,r[k]/a);state[id]=2;done++;if(done%100==0)cerr<<"replayed "<<done<<" node "<<id<<" terms "<<B[id].size()<<" reductions "<<reductions<<"\n";};
for(int id:targets)replay(id);

for(int i=0;i<(int)trace.size();i++)replay(i);
vector<int>freecols,fi(N,-1);for(int k=0;k<N;k++)if(bid[k]<0){fi[k]=freecols.size();freecols.push_back(k);}
int L=freecols.size();vector<Row>NF(N);long long nfterms=0;
for(int k=N-1;k>=0;k--){
 if(bid[k]<0){NF[k].emplace_back(fi[k],Q(1));continue;}
 vector<Q>v(L);Q product;
 for(auto&[j,c]:B[bid[k]])if(j!=k){if(j<k){cerr<<"bad triangular row\n";return 20;}
  for(auto&[l,a]:NF[j]){mpq_mul(product.get_mpq_t(),c.get_mpq_t(),a.get_mpq_t());mpq_sub(v[l].get_mpq_t(),v[l].get_mpq_t(),product.get_mpq_t());}}
 for(int i=0;i<L;i++)if(v[i]!=0)NF[k].emplace_back(i,std::move(v[i]));nfterms+=NF[k].size();
}
cerr<<"NF READY length "<<L<<" terms "<<nfterms<<"\n";
vector<vector<Row>>M(n,vector<Row>(L));for(int j=0;j<n;j++)for(int a=0;a<L;a++)if(shift[j][freecols[a]]>=0)M[j][a]=NF[shift[j][freecols[a]]];
auto action=[&](int j,const Row&v){vector<Q>out(L);Q product;for(auto&[a,c]:v)for(auto&[b,d]:M[j][a]){mpq_mul(product.get_mpq_t(),c.get_mpq_t(),d.get_mpq_t());mpq_add(out[b].get_mpq_t(),out[b].get_mpq_t(),product.get_mpq_t());}Row r;for(int i=0;i<L;i++)if(out[i]!=0)r.emplace_back(i,std::move(out[i]));return r;};
long long commutators=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)for(int a=0;a<L;a++){
 if(action(i,M[j][a])!=action(j,M[i][a])){cerr<<"commutator failed "<<i<<" "<<j<<" "<<a<<"\n";return 21;}commutators++;}
cerr<<"COMMUTATORS PASS "<<commutators<<"\n";
vector<int>parent(N,-1),which(N,-1);vector<unsigned long long>key(N);for(int k=0;k<N;k++)for(int j=0;j<n;j++)if(shift[j][k]>=0&&parent[shift[j][k]]<0){int h=shift[j][k];parent[h]=k;which[h]=j;key[h]=key[k]+(1ULL<<(4*j));}
for(int k=1;k<N;k++){if(parent[k]<0||action(which[k],NF[parent[k]])!=NF[k]){cerr<<"monomial evaluation failed "<<k<<"\n";return 22;}}
for(int j=0;j<n;j++)for(int a=0;a<L;a++)for(auto&[b,c]:M[j][a])if(deg[freecols[b]]<deg[freecols[a]]+1){cerr<<"degree decrease\n";return 23;}
for(int i:{0,1,3})for(int j:{0,1,3})if(!action(i,M[j][fi[0]]).empty()){cerr<<"base square not zero\n";return 24;}
for(int i=0;i<F;i++){vector<Q>v(L);Q product;for(auto&[k,c]:init[i])for(auto&[a,d]:NF[k]){mpq_mul(product.get_mpq_t(),c.get_mpq_t(),d.get_mpq_t());mpq_add(v[a].get_mpq_t(),v[a].get_mpq_t(),product.get_mpq_t());}for(auto&c:v)if(c!=0){cerr<<"input equation failed "<<i<<"\n";return 25;}}
cerr<<"QUOTIENT ALGEBRA PASS "<<L<<" original "<<F<<" all monomials "<<N<<"\n";
ifstream ob(argv[6]);vector<Q>tv(N);for(auto&c:tv)ob>>c;if(!ob)return 26;
vector<Q>tnf(L);for(int k=0;k<N;k++)if(tv[k]!=0)for(auto&[a,c]:NF[k])tnf[a]+=tv[k]*c;
cout<<"length "<<L<<"\nt";for(int a=0;a<L;a++)if(tnf[a]!=0)cout<<" "<<a<<":"<<tnf[a];cout<<"\n";
std::map<unsigned long long,int>ki;for(int k=0;k<N;k++)ki[key[k]]=k;
vector<Row>Mt(L);for(int a=0;a<L;a++){vector<Q>v(L);for(int k=0;k<N;k++)if(tv[k]!=0){auto it=ki.find(key[k]+key[freecols[a]]);if(it!=ki.end())for(auto&[b,c]:NF[it->second])v[b]+=tv[k]*c;}
 for(int b=0;b<L;b++)if(v[b]!=0)Mt[a].emplace_back(b,std::move(v[b]));}
vector<Q>tt(L);for(int a=0;a<L;a++)if(tnf[a]!=0)for(auto&[b,c]:Mt[a])tt[b]+=tnf[a]*c;
for(auto&c:tt)if(c!=0){cerr<<"t squared not zero\n";return 27;}
cout<<"t_squared_zero\n";
ofstream tm(argv[7]);tm<<L<<"\n";for(int a=0;a<L;a++){tm<<Mt[a].size();for(auto&[b,c]:Mt[a])tm<<" "<<b<<" "<<c;tm<<"\n";}
ifstream queries(argv[8]);ofstream answers(argv[9]);int qc;queries>>qc;answers<<qc<<"\n";
for(int qi=0;qi<qc;qi++){int length;queries>>length;vector<Q>v(L);for(int a=0;a<length;a++){int k;Q c;queries>>k>>c;if(k<0||k>=N)return 30;for(auto&[b,d]:NF[k])v[b]+=c*d;}int nz=0;for(auto&c:v)if(c!=0)nz++;answers<<nz;for(int b=0;b<L;b++)if(v[b]!=0)answers<<" "<<b<<" "<<v[b];answers<<"\n";}
answers.flush();if(!queries||!answers)return 31;
ofstream out(argv[5]);out<<nt<<"\n";for(int id:targets){auto&t=trace[id];Row r;for(auto&[k,c]:B[id]){if(deg[k]<deg[t.pivot]){cerr<<"unexpected low term "<<id<<" "<<k<<"\n";return 4;}if(deg[k]==deg[t.pivot])r.emplace_back(k,c);}out<<r.size();for(auto&[k,c]:r)out<<" "<<k<<" "<<c;out<<"\n";}cerr<<"EXACT ADAPTIVE REPLAY PASS "<<nt<<" nodes "<<done<<" reductions "<<reductions<<"\n";}
