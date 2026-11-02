// QQ replay with only degree6 projected by a previously proved homogeneous ideal.
#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>
#include <map>
#include <gmpxx.h>
using namespace std;using Q=mpq_class;using Row=vector<pair<int,Q>>;
struct Trace{int pivot,kind,source,var;};
Row readrow(istream&in){int l;in>>l;Row r;for(int k=0;k<l;k++){int j;Q c;in>>j>>c;r.emplace_back(j,c);}return r;}
int main(int argc,char**argv){if(argc!=7)return 1;
 ifstream ma(argv[1]);int n,N,lowN,F;ma>>n>>N>>lowN>>F;vector<int>deg(N);for(int&x:deg)ma>>x;
 vector<vector<Row>>shift(n,vector<Row>(N));for(auto&v:shift)for(auto&r:v)r=readrow(ma);if(!ma)return 2;
 ifstream cf(argv[2]);vector<Row>init(F);int off=0;while(deg[off]<2)off++;
 for(auto&r:init)for(int j=off;j<N;j++){Q c;cf>>c;if(c!=0)r.emplace_back(j,c);}if(!cf)return 3;
 ifstream tr(argv[3]);vector<Trace>trace;vector<int>bid(N,-1);int id;
 while(tr>>id){Trace t;int l,x;tr>>t.pivot>>t.kind>>t.source>>t.var>>l;while(l--)tr>>x;bid[t.pivot]=id;trace.push_back(t);}
 vector<Row>B(trace.size());vector<int>state(trace.size());int done=0;
 auto shifted=[&](const Row&r,int j){vector<Q>v(N);for(auto&[k,c]:r)for(auto&[l,d]:shift[j][k])v[l]+=c*d;Row out;for(int k=0;k<N;k++)if(v[k]!=0)out.emplace_back(k,std::move(v[k]));return out;};
 function<void(int)>replay=[&](int id){if(state[id]==2)return;if(state[id]){cerr<<"cycle\n";exit(4);}state[id]=1;auto&t=trace[id];Row source;if(t.kind==0)source=init[t.source];else{replay(t.source);source=shifted(B[t.source],t.var);}vector<Q>v(N);for(auto&[k,c]:source)v[k]=c;Q prod;
 for(int k=0;k<lowN;k++)if(v[k]!=0&&bid[k]>=0&&bid[k]<id){int b=bid[k];replay(b);Q c=v[k];for(auto&[l,d]:B[b]){mpq_mul(prod.get_mpq_t(),c.get_mpq_t(),d.get_mpq_t());mpq_sub(v[l].get_mpq_t(),v[l].get_mpq_t(),prod.get_mpq_t());}}
 int p=0;while(p<N&&v[p]==0)p++;if(p!=t.pivot){cerr<<"pivot mismatch "<<id<<" "<<p<<" "<<t.pivot<<"\n";exit(5);}Q c=v[p];for(int k=p;k<N;k++)if(v[k]!=0)B[id].emplace_back(k,v[k]/c);state[id]=2;done++;if(done%100==0)cerr<<"replay "<<done<<" terms "<<B[id].size()<<"\n";};
 for(int i=0;i<(int)trace.size();i++)replay(i);
 vector<int>freecols,fi(N,-1);for(int k=0;k<N;k++)if(bid[k]<0){fi[k]=freecols.size();freecols.push_back(k);}int L=freecols.size();vector<Row>NF(N);
 for(int k=N-1;k>=0;k--){if(bid[k]<0){NF[k].emplace_back(fi[k],Q(1));continue;}vector<Q>v(L);for(auto&[j,c]:B[bid[k]])if(j!=k)for(auto&[l,d]:NF[j])v[l]-=c*d;for(int j=0;j<L;j++)if(v[j]!=0)NF[k].emplace_back(j,std::move(v[j]));}
 auto normal=[&](const Row&r){vector<Q>v(L);for(auto&[j,c]:r)for(auto&[k,d]:NF[j])v[k]+=c*d;Row out;for(int j=0;j<L;j++)if(v[j]!=0)out.emplace_back(j,std::move(v[j]));return out;};
 vector<vector<Row>>M(n,vector<Row>(L));for(int j=0;j<n;j++)for(int a=0;a<L;a++)M[j][a]=normal(shift[j][freecols[a]]);
 auto action=[&](int j,const Row&r){vector<Q>v(L);for(auto&[a,c]:r)for(auto&[b,d]:M[j][a])v[b]+=c*d;Row out;for(int a=0;a<L;a++)if(v[a]!=0)out.emplace_back(a,std::move(v[a]));return out;};
 for(int j=0;j<n;j++)for(int k=j+1;k<n;k++)for(int a=0;a<L;a++)if(action(j,M[k][a])!=action(k,M[j][a])){cerr<<"commutator\n";return 6;}
 for(auto&r:init)if(!normal(r).empty()){cerr<<"input\n";return 7;}
 // Every original monomial, including every projected highest monomial,
 // is evaluated in the algebra rather than assuming projection compatibility.
 ifstream ev(argv[4]);int count;ev>>count;vector<Row>image(count);
 for(int i=0;i<count;i++){int parent,j;ev>>parent>>j;Row expected=normal(readrow(ev));Row got;if(i==0)got={{fi[0],Q(1)}};else got=action(j,image[parent]);if(got!=expected){cerr<<"evaluation "<<i<<"\n";return 8;}image[i]=std::move(got);}
 if(!ev)return 9;
 for(int j=0;j<n;j++)for(int a=0;a<L;a++)for(auto&[b,c]:M[j][a])if(deg[freecols[b]]<=deg[freecols[a]])return 10;
 ifstream query(argv[5]);ofstream answer(argv[6]);int qc;query>>qc;answer<<L<<" "<<qc<<"\n";
 for(int i=0;i<qc;i++){Row r=normal(readrow(query));answer<<r.size();for(auto&[j,c]:r)answer<<" "<<j<<" "<<c;answer<<"\n";}if(!query||!answer)return 11;
 ofstream exported(string(argv[6])+".full");exported<<L<<" "<<count<<" "<<n<<"\n";
 for(int a=0;a<L;a++)exported<<deg[freecols[a]]<<" ";exported<<"\n";
 for(auto&r:image){exported<<r.size();for(auto&[j,c]:r)exported<<" "<<j<<" "<<c;exported<<"\n";}
 for(auto&variable:M)for(auto&r:variable){exported<<r.size();for(auto&[j,c]:r)exported<<" "<<j<<" "<<c;exported<<"\n";}
 if(!exported)return 12;
 cerr<<"EXACT PROJECTED ALGEBRA PASS rows "<<done<<" length "<<L<<" monomials "<<count<<"\n";
}
