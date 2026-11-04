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
struct P{fmpq_poly_t v;P(){fmpq_poly_init(v);}P(const P&p):P(){fmpq_poly_set(v,p.v);}P(P&&p)noexcept:P(){fmpq_poly_swap(v,p.v);}P&operator=(const P&p){fmpq_poly_set(v,p.v);return *this;}P&operator=(P&&p)noexcept{fmpq_poly_swap(v,p.v);return *this;}~P(){fmpq_poly_clear(v);}};
void check(bool x,const char*s){if(!x)throw std::runtime_error(s);}
P constant(const std::string&s){P p;fmpq_t q;fmpq_init(q);check(!fmpq_set_str(q,s.c_str(),10),"QQ input");fmpq_canonicalise(q);fmpq_poly_set_fmpq(p.v,q);fmpq_clear(q);return p;}
P pm(const P&a,const P&b){P p;fmpq_poly_mul(p.v,a.v,b.v);return p;}
P pa(const P&a,const P&b){P p;fmpq_poly_add(p.v,a.v,b.v);return p;}
P pn(const P&a){P p;fmpq_poly_neg(p.v,a.v);return p;}
P divide_scalar(const P&a,const fmpq_t c){P p;fmpq_poly_scalar_div_fmpq(p.v,a.v,c);return p;}
P readP(){int n;check(bool(std::cin>>n),"EOF");P p;fmpq_t q;fmpq_init(q);for(int j=0;j<n;j++){int e;std::string v;std::cin>>e>>v;check(e>=0&&!fmpq_set_str(q,v.c_str(),10),"poly input");fmpq_canonicalise(q);fmpq_poly_set_coeff_fmpq(p.v,e,q);}fmpq_clear(q);return p;}
std::array<P,5> factors;
std::array<std::vector<P>,5> powers;
const P& fp(int k,int e){check(e>=0&&e<120,"unit exponent bound");while((int)powers[k].size()<=e)powers[k].push_back(pm(powers[k].back(),factors[k]));return powers[k][e];}
std::array<int,5> unit_factor(P p,fmpq_t lc){
 check(!fmpq_poly_is_zero(p.v),"nonzero unit factor");fmpq_poly_get_coeff_fmpq(lc,p.v,fmpq_poly_degree(p.v));p=divide_scalar(p,lc);std::array<int,5> e{};
 for(int k=0;k<5;k++)for(;;){P q,r;fmpq_poly_divrem(q.v,r.v,p.v,factors[k].v);if(!fmpq_poly_is_zero(r.v))break;p=std::move(q);e[k]++;}
 check(fmpq_poly_is_one(p.v),"denominator is not product of named whole-curve units");return e;
}
struct R{
 P n;std::array<int,5> e{};
 R(){}R(P p):n(std::move(p)){}R(P p,std::array<int,5> a):n(std::move(p)),e(a){normalize();}
 void normalize(){if(fmpq_poly_is_zero(n.v)){e={};return;}for(int k=0;k<5;k++)while(e[k]){P q,r;fmpq_poly_divrem(q.v,r.v,n.v,factors[k].v);if(!fmpq_poly_is_zero(r.v))break;n=std::move(q);e[k]--;}}
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
int main(){try{
 auto start=std::chrono::steady_clock::now();
 for(int i=0;i<7;i++)for(int j=i;j<7;j++)pairs.push_back({i,j});
 for(int k=0;k<5;k++){factors[k]=readP();powers[k]={constant("1")};}
 P z;fmpq_poly_set_coeff_si(z.v,1,1);P zp20=pa(z,constant("20"));
 R a(pm(constant("-10/9"),zp20),{2,0,0,0,0}),e(pm(constant("1/9"),pa(z,constant("100"))),{1,0,0,0,0});
 R d(pm(constant("-3/100"),pm(z,zp20)));ap={rc("1"),a};ep={rc("1"),e};dp={rc("1"),d};
 std::vector<V> raw(16,V(7));for(auto&row:raw)for(auto&v:row)v=eval(readS());
 auto oldmons=readMons(29);std::vector<V> qold(6,V(28));std::vector<int> phys={0,1,2,3,4,6,7};
 for(int i=0;i<6;i++)for(int k=0;k<29;k++){Source p=readS();if(oldmons[k][5]){check(p.empty(),"original quadratic H coefficient literally zero");continue;}auto ij=pair_of(oldmons[k],phys);qold[i][mi(ij[0],ij[1])]=eval(p);}
 std::vector<V> forward(7,V(7));for(auto&row:forward)for(auto&v:row)v=eval(readS());
 auto newmons=readMons(29);std::vector<std::vector<Source>> primenum(6,std::vector<Source>(29));for(auto&row:primenum)for(auto&v:row)v=readS();
 std::vector<R> den(6);for(auto&v:den)v=inv_unit(eval(readS()));
 std::vector<V> prime(6,V(28));for(int i=0;i<6;i++)for(int k=0;k<29;k++){if(newmons[k][7]){check(primenum[i][k].empty(),"corrected prime H coefficient literally zero");continue;}auto ij=pair_of(newmons[k],{0,1,2,3,4,5,6});prime[i][mi(ij[0],ij[1])]=rm(eval(primenum[i][k]),den[i]);}
 std::vector<std::vector<int>> H(7,std::vector<int>(7)),Hi=H;for(auto&row:H)for(int&v:row)std::cin>>v;for(auto&row:Hi)for(int&v:row)std::cin>>v;
 for(int i=0;i<7;i++)for(int j=0;j<7;j++){int xy=0,yx=0;for(int k=0;k<7;k++){xy+=H[i][k]*Hi[k][j];yx+=Hi[i][k]*H[k][j];}check(xy==(i==j)&&yx==(i==j),"constant GL7 inverse");}
 auto hmons=readMons(23);V eta(7);for(auto&v:eta)v=readR();
 check(eq(eta[3],rc("1")),"physical W normalization");
 for(const auto&row:raw){R x;for(int j=0;j<7;j++)x=ra(x,rm(row[j],eta[j]));check(zero(x),"all16 literal raw dual annihilations");}
 V qvalues(6);for(int i=0;i<6;i++){R x;for(int k=0;k<28;k++)x=ra(x,rm(qold[i][k],rm(eta[pairs[k][0]],eta[pairs[k][1]])));qvalues[i]=readR();check(eq(x,qvalues[i]),"all6 original q eta values");}
 R u=R(constant("137500/243"),{2,0,0,0,0});check(eq(qvalues[1],u),"unit q1 eta value");R iu=inv_unit(u);
 std::vector<V> sigma(7,V(16));std::vector<int> basis={1,6,8,9,11,14};
 for(int j=0;j<7;j++)for(int s=0;s<16;s++){sigma[j][s]=readR();if(!zero(sigma[j][s]))check(std::find(basis.begin(),basis.end(),s)!=basis.end(),"fixed six source basis");}
 for(int i=0;i<7;i++)for(int k=0;k<7;k++){R x;for(int s=0;s<16;s++)x=ra(x,rm(sigma[i][s],raw[s][k]));R target=rc(i==k?"1":"0");if(k==3)target=rs(target,eta[i]);check(eq(x,target),"all49 source lift slots");}
 std::vector<V> lambda(7,V(7));for(int i=0;i<7;i++)for(int p=0;p<7;p++)for(int k=0;k<7;k++)lambda[i][p]=ra(lambda[i][p],rm(forward[i][k],rc(std::to_string(Hi[p][k]))));
 V joint_eta(7);for(int i=0;i<7;i++)for(int p=0;p<7;p++)joint_eta[i]=ra(joint_eta[i],rm(lambda[i][p],eta[p]));check(eq(joint_eta[6],rc("1")),"joint W normalization");
 std::vector<int> piv={mi(0,0),mi(0,2),mi(1,5),mi(1,6),mi(1,1)};std::vector<V> five(prime.begin(),prime.begin()+5);
 for(int r=0;r<5;r++){if(r<4)check(eq(five[r][piv[r]],rc("1")),"fixed first four monic pivots");five[r]=vs(five[r],inv_unit(five[r][piv[r]]));V expected(28);for(auto&v:expected)v=readR();check(veq(five[r],expected),"all5 corrected normalized NF rows");for(int t=r+1;t<5;t++)five[t]=va(five[t],vs(five[r],rn(five[t][piv[r]])));}
 auto nf=[&](V v){for(int r=0;r<5;r++)v=va(v,vs(five[r],rn(v[piv[r]])));return v;};
 std::vector<V> orig(5,V(28));orig[0]=prime[4];orig[1]=vs(prime[0],rc("1/8100"));orig[3]=va(vs(prime[0],rm(e,rc("1/8100"))),vs(prime[1],rm(a,rc("-400/8100"))));
 orig[2]=vs(va(va(prime[3],vs(orig[1],ra(rm(d,e),rc("8")))),vs(orig[3],rn(d))),rc("1/8"));
 orig[4]=vs(va(va(vs(orig[1],rs(ra(rm(rc("10"),a),rm(rc("3"),rm(e,e))),rm(rc("3"),e))),vs(orig[3],rn(rm(rc("3"),rs(e,rc("1")))))),vs(prime[2],rn(a))),inv_unit(rm(rc("10"),a)));
 for(int i=0;i<5;i++){check(veq(push(orig[i],lambda),qold[i]),"full original-to-corrected-prime source transformation");for(auto&v:nf(orig[i]))check(zero(v),"all5 original qH normal forms zero");}
 std::vector<int> retained;for(int i=0;i<28;i++)if(std::find(piv.begin(),piv.end(),i)==piv.end())retained.push_back(i);check(retained.size()==23,"retained H1 count");
 std::vector<V> L(112,V(23));for(auto&row:L)for(auto&v:row)v=readR();V rho(23);for(auto&v:rho)v=readR();
 R delta=eval(readS());R cc=rs(rs(ra(rm(rc("100"),a),rm(rc("18"),rm(e,e))),rm(rc("9"),e)),rc("2"));check(eq(delta,rm(rc("-64000"),rm(cc,cc))),"actual delta7 unit binding");R idelta=inv_unit(delta);
 std::vector<V> newraw=raw;for(int k=0;k<7;k++)newraw[6][k]=rm(delta,rs(raw[6][k],raw[1][k]));
 long checks=0;
 for(int t=0;t<23;t++){
   check(hmons[t][7]==1,"H degree1 target");auto ij=pair_of(hmons[t],{0,1,2,3,4,5,6});check(mi(ij[0],ij[1])==retained[t],"actual retained target order");
   V target(28);for(int j=0;j<7;j++)for(int k=0;k<7;k++)target[mi(j,k)]=ra(target[mi(j,k)],rm(lambda[ij[0]][j],lambda[ij[1]][k]));
   check(eq(rho[t],rm(rm(joint_eta[ij[0]],joint_eta[ij[1]]),iu)),"q1 correction rho source");
   V cn(112),co(112);for(int s=0;s<112;s++)co[s]=cn[s]=L[s][t];
   for(int j=0;j<7;j++){co[6*7+j]=rm(delta,cn[6*7+j]);co[1*7+j]=rs(cn[1*7+j],co[6*7+j]);check(eq(rm(co[6*7+j],idelta),cn[6*7+j]),"old/new source unit transport");}
   for(int frame=0;frame<2;frame++){
     const V&coef=frame?cn:co;const auto& rr=frame?newraw:raw;V actual=vs(qold[1],rho[t]);
     for(int s=0;s<16;s++)for(int j=0;j<7;j++)if(!zero(coef[7*s+j]))for(int k=0;k<7;k++)actual[mi(j,k)]=ra(actual[mi(j,k)],rm(coef[7*s+j],rr[s][k]));
     for(int k=0;k<28;k++){check(eq(actual[k],target[k]),"all23 full source quadratic identities in both frames");checks++;}
   }
 }
 std::string extra;check(!(std::cin>>extra),"trailing input");
 std::cout<<"{\"all_passed\":true,\"raw16_dual_annihilations\":16,\"original_q_eta_values\":6,\"linear_lift_slots\":49,\"all5_corrected_Hquad_NF_relations\":true,\"all5_original_q_source_transforms\":true,\"H1_targets\":23,\"physical_coefficient_checks_both_frames\":"<<checks<<",\"all_serialized_denominators_named_curve_units\":true,\"no_GCD_Gaussian_unknown_inverse_or_optional_row_query\":true,\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
