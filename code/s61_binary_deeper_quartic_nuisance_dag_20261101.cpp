// Exact quartic certificate: nuisance dag.
#include "s61_binary_deeper_quartic_field_20261101.hpp"
#include <array>
#include <set>
using Mon=array<int,5>;using OP=map<Mon,Rat>;Mon zero{};
Rat neg(const Rat&a){Rat c=a;fmpq_mpoly_neg(c.n,c.n,ctx);return c;}
void addto(OP&x,const Mon&m,const Rat&v){Rat z=sub(x[m],neg(v));if(z.zero())x.erase(m);else x[m]=std::move(z);}
OP times(const OP&a,const OP&b){OP c;for(auto&[m,v]:a)for(auto&[n,w]:b){Mon q;for(int j=0;j<5;j++)q[j]=m[j]+n[j];addto(c,q,mul(v,w));}return c;}
void subtractfrom(OP&a,const OP&b){for(auto&[m,v]:b)addto(a,m,neg(v));}
int main(int argc,char**argv){fmpq_mpoly_ctx_init(ctx,2,ORD_LEX);flint_set_num_threads(2);{ifstream in(argv[1]),inv(argv[2]);int n,nn,nc;in>>n;inv>>nn>>nc;assert(n==24&&nn==n&&nc==n);vector<vector<OP>> B(n,vector<OP>(n)),N=B,X=B;vector<Row> BI(n);for(int i=0;i<n;i++)for(int j=0;j<n;j++){int l;in>>l;for(int k=0;k<l;k++){Mon m;for(int&v:m)in>>v;Rat c;readpoly(in,c.n);readpoly(in,c.d);c.reduce();B[i][j].emplace(m,std::move(c));}}for(int i=0;i<n;i++){int l;inv>>l;for(int k=0;k<l;k++){int j;inv>>j;Rat c;readpoly(inv,c.n);readpoly(inv,c.d);c.reduce();BI[i].emplace(j,std::move(c));}}
Rat one;fmpq_mpoly_one(one.n,ctx);for(int i=0;i<n;i++)for(auto&[k,v]:BI[i])for(int j=0;j<n;j++)for(auto&[m,w]:B[k][j])addto(N[i][j],m,mul(v,w));for(int i=0;i<n;i++)addto(N[i][i],zero,neg(one));for(int i=0;i<n;i++)for(int j=0;j<n;j++)assert(!N[i][j].count(zero));cerr<<"CORE INVERSE VERIFIED\n";
vector<int>order;set<int>remaining;for(int i=0;i<n;i++)remaining.insert(i);while(!remaining.empty()){vector<int>ready;for(int i:remaining){bool good=true;for(int j:remaining)if(!N[i][j].empty())good=false;if(good)ready.push_back(i);}assert(!ready.empty());for(int i:ready){order.push_back(i);remaining.erase(i);}}
int edges=0;for(auto&row:N)for(auto&p:row)edges+=!p.empty();cerr<<"DAG EDGES "<<edges<<"\n";assert(edges==144);
for(int i:order){for(auto&[j,v]:BI[i])X[i][j][zero]=v;for(int k=0;k<n;k++)if(!N[i][k].empty())for(int j=0;j<n;j++)if(!X[k][j].empty())subtractfrom(X[i][j],times(N[i][k],X[k][j]));cerr<<"INVERSE ROW "<<i<<"\n";}
for(int i=0;i<n;i++){for(int j=0;j<n;j++){OP value;for(int k=0;k<n;k++)if(!X[i][k].empty()&&!B[k][j].empty())subtractfrom(value,times(X[i][k],B[k][j]));if(i==j)addto(value,zero,one);assert(value.empty());}cerr<<"CHECK ROW "<<i<<"\n";}
ofstream out(argv[3]);assert(out);out<<n<<'\n';for(auto&row:X)for(auto&p:row){out<<p.size()<<'\n';for(auto&[m,v]:p){for(int x:m)out<<x<<' ';out<<'\n';writepoly(out,v.n);writepoly(out,v.d);}}ofstream info(string(argv[3])+".info");assert(info);info<<edges<<'\n';for(int i:order)info<<i<<' ';info<<'\n';}fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
