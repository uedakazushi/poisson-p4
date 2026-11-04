// Independent fixed QQ[z] source/identity checker. Denominators use only five
// named linear units. No rational-function gcd, basis solve or Gaussian search.
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include <set>
long LIVE=0,PEAK=0;const long LIVE_CAP=150000;
struct P{
 fmpq_poly_t v;mutable long counted=0;
 void sync()const{long n=fmpq_poly_length(v);LIVE+=n-counted;counted=n;PEAK=std::max(PEAK,LIVE);if(LIVE>LIVE_CAP)throw std::runtime_error("live_RF_polynomial_slots_150000");}
 P(){fmpq_poly_init(v);}P(const P&p):P(){fmpq_poly_set(v,p.v);sync();}
 P(P&&p)noexcept:P(){fmpq_poly_swap(v,p.v);counted=p.counted;p.counted=0;}
 P&operator=(const P&p){fmpq_poly_set(v,p.v);sync();return *this;}
 P&operator=(P&&p)noexcept{fmpq_poly_swap(v,p.v);std::swap(counted,p.counted);return *this;}
 ~P(){LIVE-=counted;fmpq_poly_clear(v);}
};
void check(bool x,const char*s){if(!x)throw std::runtime_error(s);}
P constant(const std::string&s){P p;fmpq_t q;fmpq_init(q);check(!fmpq_set_str(q,s.c_str(),10),"QQ input");fmpq_canonicalise(q);fmpq_poly_set_fmpq(p.v,q);fmpq_clear(q);p.sync();return p;}
P pm(const P&a,const P&b){P p;fmpq_poly_mul(p.v,a.v,b.v);p.sync();return p;}
P pa(const P&a,const P&b){P p;fmpq_poly_add(p.v,a.v,b.v);p.sync();return p;}
P pn(const P&a){P p;fmpq_poly_neg(p.v,a.v);p.sync();return p;}
P divide_scalar(const P&a,const fmpq_t c){P p;fmpq_poly_scalar_div_fmpq(p.v,a.v,c);p.sync();return p;}
P readP(){int n;check(bool(std::cin>>n),"EOF");P p;fmpq_t q;fmpq_init(q);for(int j=0;j<n;j++){int e;std::string v;std::cin>>e>>v;check(e>=0&&!fmpq_set_str(q,v.c_str(),10),"poly input");fmpq_canonicalise(q);fmpq_poly_set_coeff_fmpq(p.v,e,q);}fmpq_clear(q);p.sync();return p;}
std::array<P,5> factors;
std::array<std::vector<P>,5> powers;
const P& fp(int k,int e){check(e>=0&&e<120,"unit exponent bound");while((int)powers[k].size()<=e)powers[k].push_back(pm(powers[k].back(),factors[k]));return powers[k][e];}
std::array<int,5> unit_factor(P p,fmpq_t lc){
 check(!fmpq_poly_is_zero(p.v),"nonzero unit factor");fmpq_poly_get_coeff_fmpq(lc,p.v,fmpq_poly_degree(p.v));p=divide_scalar(p,lc);std::array<int,5> e{};
 for(int k=0;k<5;k++)for(;;){P q,r;fmpq_poly_divrem(q.v,r.v,p.v,factors[k].v);q.sync();r.sync();if(!fmpq_poly_is_zero(r.v))break;p=std::move(q);e[k]++;}
 check(fmpq_poly_is_one(p.v),"denominator is not product of named whole-curve units");return e;
}
struct R{
 P n;std::array<int,5> e{};
 R(){}R(P p):n(std::move(p)){}R(P p,std::array<int,5> a):n(std::move(p)),e(a){normalize();}
 void normalize(){if(fmpq_poly_is_zero(n.v)){e={};return;}for(int k=0;k<5;k++)while(e[k]){P q,r;fmpq_poly_divrem(q.v,r.v,n.v,factors[k].v);q.sync();r.sync();if(!fmpq_poly_is_zero(r.v))break;n=std::move(q);e[k]--;}}
};
R rc(const std::string&s){return R(constant(s));}
bool zero(const R&r){return fmpq_poly_is_zero(r.n.v);}
R rn(const R&a){return R(pn(a.n),a.e);}
R rm(const R&a,const R&b){std::array<int,5> e;for(int k=0;k<5;k++)e[k]=a.e[k]+b.e[k];return R(pm(a.n,b.n),e);}
R ra(const R&a,const R&b){if(zero(a))return b;if(zero(b))return a;std::array<int,5> e;P x=a.n,y=b.n;for(int k=0;k<5;k++){e[k]=std::max(a.e[k],b.e[k]);x=pm(x,fp(k,e[k]-a.e[k]));y=pm(y,fp(k,e[k]-b.e[k]));}return R(pa(x,y),e);}
R rs(const R&a,const R&b){return ra(a,rn(b));}
bool eq(const R&a,const R&b){return a.e==b.e&&fmpq_poly_equal(a.n.v,b.n.v);}
R inv_unit(const R&a){check(!zero(a),"nonzero unit");fmpq_t lc;fmpq_init(lc);auto e=unit_factor(a.n,lc);P n=constant("1");for(int k=0;k<5;k++)n=pm(n,fp(k,a.e[k]));n=divide_scalar(n,lc);fmpq_clear(lc);return R(n,e);}
R readR(){P n=readP(),d=readP();fmpq_t lc;fmpq_init(lc);auto e=unit_factor(d,lc);n=divide_scalar(n,lc);fmpq_clear(lc);return R(n,e);}
using V=std::vector<R>;
V va(const V&a,const V&b){V c=a;for(size_t k=0;k<c.size();k++)c[k]=ra(c[k],b[k]);return c;}
V vs(const V&a,const R&b){V c=a;for(auto&v:c)v=rm(v,b);return c;}
bool veq(const V&a,const V&b){if(a.size()!=b.size())return false;for(size_t i=0;i<a.size();i++)if(!eq(a[i],b[i]))return false;return true;}
struct T{int a,e,d;R c;};using Source=std::vector<T>;
Source readS(){int n;std::cin>>n;Source s;for(int j=0;j<n;j++){int a,e,d;std::string c;std::cin>>a>>e>>d>>c;check(a>=0&&e>=0&&d>=0,"source exponent");s.push_back({a,e,d,rc(c)});}return s;}
std::vector<R> ap,ep,dp;std::map<std::tuple<int,int,int>,R> monomials;
R eval(const Source&s){R z;for(const auto&t:s){while((int)ap.size()<=t.a)ap.push_back(rm(ap.back(),ap[1]));while((int)ep.size()<=t.e)ep.push_back(rm(ep.back(),ep[1]));while((int)dp.size()<=t.d)dp.push_back(rm(dp.back(),dp[1]));auto key=std::make_tuple(t.a,t.e,t.d);auto it=monomials.find(key);if(it==monomials.end())it=monomials.emplace(key,rm(rm(ap[t.a],ep[t.e]),dp[t.d])).first;z=ra(z,rm(t.c,it->second));}return z;}
std::vector<std::array<int,2>> pairs;
int mi(int i,int j){if(i>j)std::swap(i,j);for(int k=0;k<28;k++)if(pairs[k][0]==i&&pairs[k][1]==j)return k;throw std::runtime_error("pair");}
using Mons=std::array<int,8>;
std::vector<Mons> readMons(int n){std::vector<Mons>a(n);for(auto&m:a)for(int&v:m)std::cin>>v;return a;}
std::array<int,2> pair_of(const Mons&m,const std::vector<int>&physical){std::vector<int> ix;for(int j=0;j<7;j++)for(int k=0;k<m[physical[j]];k++)ix.push_back(j);check(ix.size()==2,"two light indices");return {ix[0],ix[1]};}
V push(const V&q,const std::vector<V>&lambda){V v(28);for(int m=0;m<28;m++)if(!zero(q[m])){int i=pairs[m][0],j=pairs[m][1];for(int a=0;a<7;a++)for(int b=0;b<7;b++)v[mi(a,b)]=ra(v[mi(a,b)],rm(q[m],rm(lambda[i][a],lambda[j][b])));}return v;}

// New independent fixed pure-row verification follows.
using Mon=std::array<int,7>;
struct Node{long id;char op;long a=0,b=0;std::vector<std::pair<long,long>> terms;};
std::vector<Node> readNodes(){long count;std::cin>>count;std::vector<Node> out;for(long k=0;k<count;k++){
 Node n;std::cin>>n.id>>n.op;
 if(n.op=='c'||n.op=='v')std::cin>>n.a;
 else if(n.op=='m'||n.op=='a'||n.op=='s')std::cin>>n.a>>n.b;
 else if(n.op=='l'){long t;std::cin>>t;for(long i=0;i<t;i++){long c,j;std::cin>>c>>j;n.terms.emplace_back(c,j);}}
 else check(n.op=='z',"typed node opcode");out.push_back(std::move(n));
}return out;}
std::vector<long> children(const Node&n){if(n.op=='a'||n.op=='s'||n.op=='m')return{n.a,n.b};std::vector<long> out;if(n.op=='l')for(auto p:n.terms)out.push_back(p.second);return out;}
R value(const Node&n,const std::map<long,R>& v,const std::map<long,R>& c){
 if(n.op=='z')return R();if(n.op=='c')return c.at(n.a);
 if(n.op=='a')return ra(v.at(n.a),v.at(n.b));if(n.op=='s')return rs(v.at(n.a),v.at(n.b));if(n.op=='m')return rm(v.at(n.a),v.at(n.b));
 check(n.op=='l',"scalar/linear source node");R r;for(auto p:n.terms)r=ra(r,rm(c.at(p.first),v.at(p.second)));return r;
}
std::map<long,R> readAtoms(){long count;std::cin>>count;std::map<long,R> out;for(long k=0;k<count;k++){long i;std::cin>>i;Source n=readS(),d=readS();R nn=eval(n),dd=eval(d);out.emplace(i,rm(nn,inv_unit(dd)));}monomials.clear();return out;}
std::vector<std::vector<R>> readRMatrix(int n,int m){std::vector<std::vector<R>> a(n,std::vector<R>(m));for(auto&row:a)for(auto&v:row)v=readR();return a;}
std::vector<Mon> readMon(int count){std::vector<Mon> m(count);for(auto&e:m)for(int&v:e)std::cin>>v;return m;}
P denominator(const R&r){P d=constant("1");for(int k=0;k<5;k++)d=pm(d,fp(k,r.e[k]));return d;}
bool named_unit(const P&input){if(fmpq_poly_is_zero(input.v))return false;P p=input;for(int k=0;k<5;k++)for(;;){P q,r;fmpq_poly_divrem(q.v,r.v,p.v,factors[k].v);q.sync();r.sync();if(!fmpq_poly_is_zero(r.v))break;p=std::move(q);}return fmpq_poly_degree(p.v)==0;}
int main(){try{
 auto start=std::chrono::steady_clock::now();
 for(int k=0;k<5;k++){factors[k]=readP();powers[k]={constant("1")};}
 P z;fmpq_poly_set_coeff_si(z.v,1,1);z.sync();P zp20=pa(z,constant("20"));
 R a(pm(constant("-10/9"),zp20),{2,0,0,0,0}),e(pm(constant("1/9"),pa(z,constant("100"))),{1,0,0,0,0});
 R d(pm(constant("-3/100"),pm(z,zp20)));ap={rc("1"),a};ep={rc("1"),e};dp={rc("1"),d};
 auto ca=readAtoms();auto nodes=readNodes();check(nodes.size()==19960,"pinned coefficient DAG count");
 std::vector<std::vector<long>> roots(6,std::vector<long>(84));for(auto&row:roots)for(long&i:row)std::cin>>i;
 auto oldmon=readMon(84);std::map<long,long> uses;
 for(auto&n:nodes){uses[n.id]+=0;for(long c:children(n))uses[c]++;}for(auto&row:roots)for(long r:row)uses[r]++;
 std::map<long,R> vals;
 for(auto&n:nodes){vals.emplace(n.id,value(n,vals,ca));for(long c:children(n))if(--uses[c]==0)vals.erase(c);}
 std::vector<std::vector<R>> oldpi(6,std::vector<R>(84));
 for(int i=0;i<6;i++)for(int j=0;j<84;j++){oldpi[i][j]=vals.at(roots[i][j]);check(eq(oldpi[i][j],readR()),"all504 normalized coefficient values");if(--uses[roots[i][j]]==0)vals.erase(roots[i][j]);}
 vals.clear();ca.clear();nodes.clear();
 std::vector<R>b(23);for(auto&v:b){Source n=readS(),q=readS();v=rm(eval(n),inv_unit(eval(q)));check(eq(v,readR()),"actual q5H cache0 values23");}monomials.clear();
 long active;std::cin>>active;check(active==42,"exact42 source rows");std::map<int,R> knew,kold;
 for(long k=0;k<active;k++){int i;std::cin>>i;R x,y;std::vector<R>cofs(23);for(auto&v:cofs)v=readR();for(int j=0;j<23;j++)x=ra(x,rm(cofs[j],b[j]));for(int j=22;j>=0;j--)y=ra(y,rm(cofs[j],b[j]));check(eq(x,y),"42 cofactor sums in both orders");knew[i]=x;check(eq(x,readR()),"saved new source cofactor");}
 R cc=rs(rs(ra(rm(rc("100"),a),rm(rc("18"),rm(e,e))),rm(rc("9"),e)),rc("2"));R delta=rm(rc("-64000"),rm(cc,cc));R idelta=inv_unit(delta);kold=knew;
 for(int j=0;j<7;j++){kold[42+j]=rm(delta,knew.at(42+j));kold[7+j]=rs(knew.at(7+j),kold.at(42+j));check(eq(rm(kold.at(42+j),idelta),knew.at(42+j)),"delta inverse source coefficient");check(eq(ra(kold.at(7+j),kold.at(42+j)),knew.at(7+j)),"delta old/new source sum");}
 for(auto&pair:kold){int i;std::cin>>i;check(i==pair.first&&eq(pair.second,readR()),"all42 old source coefficients");}knew.clear();b.clear();
 std::vector<V> backward(7,V(7)),forward(7,V(7));
 R jac=rm(rc("3200/81"),rm(rm(a,cc),rm(d,d)));R ijac=inv_unit(jac);
 std::vector<V> n49(7,V(7));for(auto&row:n49)for(auto&v:row)v=rm(eval(readS()),ijac);
 int order[7]={0,1,2,6,3,5,4};for(int i=0;i<7;i++)for(int j=0;j<7;j++)backward[i][j]=i==5?rs(ra(n49[5][j],n49[1][j]),rm(rc("3"),n49[0][j])):n49[order[i]][j];n49.clear();
 std::vector<V> ff(7,V(7));for(auto&row:ff)for(auto&v:row)v=eval(readS());
 std::vector<std::vector<int>> Hi(7,std::vector<int>(7));for(auto&row:Hi)for(int&v:row)std::cin>>v;
 for(int i=0;i<7;i++)for(int j=0;j<7;j++)for(int k=0;k<7;k++)forward[i][j]=ra(forward[i][j],rm(ff[i][k],rc(std::to_string(Hi[j][k]))));ff.clear();
 for(int i=0;i<7;i++)for(int j=0;j<7;j++){R x,y;for(int k=0;k<7;k++){x=ra(x,rm(forward[i][k],backward[k][j]));y=ra(y,rm(backward[i][k],forward[k][j]));}check(eq(x,rc(i==j?"1":"0"))&&eq(y,rc(i==j?"1":"0")),"both49 physical/joint inverse products");}forward.clear();monomials.clear();
 std::vector<Mon> mon3=oldmon;std::sort(mon3.begin(),mon3.end());std::map<Mon,int> index3;for(int j=0;j<84;j++)index3[mon3[j]]=j;
 std::vector<std::vector<R>> joint(6,std::vector<R>(84));
 for(int m=0;m<84;m++){
   std::map<Mon,R> image;Mon zeroM{};image.emplace(zeroM,rc("1"));
   for(int j=0;j<7;j++)for(int count=0;count<oldmon[m][j];count++){
     std::map<Mon,R> next;for(auto&pair:image)for(int k=0;k<7;k++)if(!zero(backward[j][k])){Mon ex=pair.first;ex[k]++;R v=rm(pair.second,backward[j][k]);auto it=next.find(ex);if(it==next.end())next.emplace(ex,std::move(v));else it->second=ra(it->second,v);}image=std::move(next);
   }
   for(auto&pair:image)for(int s=0;s<6;s++)joint[s][index3.at(pair.first)]=ra(joint[s][index3.at(pair.first)],rm(oldpi[s][m],pair.second));
 }
 oldpi.clear();for(auto&row:joint)for(auto&v:row)check(eq(v,readR()),"all504 joint cubic coefficients");
 int genuine[6]={1,6,8,9,11,14};std::vector<V> li(6,V(7));for(int s=0;s<6;s++)for(int j=0;j<7;j++)for(int k=0;k<7;k++)li[s][k]=ra(li[s][k],rm(kold.at(7*genuine[s]+j),backward[j][k]));
 backward.clear();kold.clear();auto mon4=readMon(210);std::map<Mon,int> index4;for(int m=0;m<210;m++)index4[mon4[m]]=m;
 std::vector<R> psi(210);
 for(int m=0;m<210;m++){
   R v;for(int j=0;j<7;j++)if(mon4[m][j]){Mon sub=mon4[m];sub[j]--;int c=index3.at(sub);for(int s=0;s<6;s++)v=rs(v,rm(joint[s][c],li[s][j]));}
   psi[m]=v;check(eq(v,readR()),"all210 actual negative full-tail combinations");
 }
 joint.clear();li.clear();auto nfAtoms=readAtoms();auto scalarNodes=readNodes();check(scalarNodes.size()==244,"pure98 scalar closure");std::map<long,R> scalar;
 for(auto&n:scalarNodes)scalar.emplace(n.id,value(n,scalar,nfAtoms));nfAtoms.clear();scalarNodes.clear();
 auto linearNodes=readNodes();check(linearNodes.size()==547,"pure98 linear closure");std::vector<long> outroots(98),labels(98),inputvars(239);for(auto&i:outroots)std::cin>>i;for(auto&i:labels)std::cin>>i;for(auto&i:inputvars)std::cin>>i;
 std::map<long,int> variablecol;for(int i=0;i<239;i++)variablecol[inputvars[i]]=i;
 std::vector<std::array<int,8>> allmon(239);for(auto&m:allmon)for(int&i:m)std::cin>>i;
 std::map<long,long> lu;for(auto&n:linearNodes){lu[n.id]+=0;for(auto i:children(n))lu[i]++;}for(long i:outroots)lu[i]++;
 std::map<long,R> linear;
 for(auto&n:linearNodes){R v;if(n.op=='v'){int i=variablecol.at(n.a);check(allmon[i][7]==0,"pure NF input H degree");Mon m;std::copy_n(allmon[i].begin(),7,m.begin());v=psi[index4.at(m)];}else v=value(n,linear,scalar);linear.emplace(n.id,std::move(v));for(auto i:children(n))if(--lu[i]==0)linear.erase(i);}
 scalar.clear();psi.clear();std::vector<R> out(98);long individualUnits=0;
 for(int i=0;i<98;i++){out[i]=linear.at(outroots[i]);check(eq(out[i],readR()),"all98 actual NF outputs");check(!zero(out[i]),"all98 source outputs nonzero");individualUnits+=named_unit(out[i].n);if(--lu[outroots[i]]==0)linear.erase(outroots[i]);}
 check(individualUnits==0,"individual-unit metadata");linear.clear();
 P g=readP();check(fmpq_poly_is_one(g.v),"saved common numerator gcd1");
 for(int i=0;i<98;i++){P q=readP();check(fmpq_poly_equal(pm(q,g).v,out[i].n.v),"all98 gcd quotient identities");}
 int contributors;std::cin>>contributors;check(contributors==2,"two Bezout contributors");P bez;R dualsum;
 for(int k=0;k<2;k++){
   int pos;std::cin>>pos;check(pos==k,"fixed contributor positions0/1");P A=readP();bez=pa(bez,pm(A,out[pos].n));R dual=readR();P actual=pm(A,denominator(out[pos]));check(eq(dual,R(actual)),"regular row dual AiDi");dualsum=ra(dualsum,rm(dual,out[pos]));
 }
 check(fmpq_poly_is_one(bez.v)&&eq(dualsum,rc("1")),"exact polynomial Bezout and regular pure-row dual");
 std::string extra;check(!(std::cin>>extra),"trailing input");
 std::cout<<"{\"all_passed\":true,\"source504_values_checked\":true,\"cofactor42_source_transport_checked\":true,\"all504_joint_coefficients_checked\":true,\"all210_actual_combinations_checked\":true,\"all98_NF_values_checked\":true,\"both49_coordinate_inverse_products_checked\":true,\"literal_Bezout_positions0_1_equals1\":true,\"regular_row_dual_equals1\":true,\"individual_unit_entries\":0,\"live_polynomial_slots_peak\":"<<PEAK<<",\"no_GCD_rank_search_new_source_or_field_query\":true,\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
