// ONE all-H15 finite algebra test, strict coefficient/rank lifetime separation.
#include <flint/flint.h>
#include <flint/nmod_poly.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
constexpr uint32_t p=32003;int D=1;array<uint32_t,4>modc{};long live=0,peak=0,cap=107743,constructs=0;double limit=50;string stage="start",folder;chrono::steady_clock::time_point began;
ofstream inverse_log;long inverse_checks=0;
void guard(){peak=max(peak,live*D);if(live*D>cap)throw runtime_error("field_slot_guard_"+stage);if((++constructs%100000)==0&&chrono::duration<double>(chrono::steady_clock::now()-began).count()>limit)throw runtime_error("time_guard_"+stage);}
struct F{array<uint32_t,4>x{};F(){++live;guard();}F(uint32_t n):F(){x[0]=n%p;}F(const F&a):F(){x=a.x;}F(F&&a):F(){x=a.x;}~F(){--live;}F&operator=(const F&a){x=a.x;return*this;}F&operator=(F&&a){x=a.x;return*this;}};
bool zero(const F&a){for(int i=0;i<D;i++)if(a.x[i])return false;return true;}bool eq(const F&a,const F&b){return a.x==b.x;}
F add(const F&a,const F&b){F c;for(int i=0;i<D;i++)c.x[i]=(a.x[i]+b.x[i])%p;return c;}F neg(const F&a){F c;for(int i=0;i<D;i++)c.x[i]=a.x[i]?p-a.x[i]:0;return c;}F sub(const F&a,const F&b){return add(a,neg(b));}
F mul(const F&a,const F&b){uint64_t c[7]={};for(int i=0;i<D;i++)for(int j=0;j<D;j++)c[i+j]=(c[i+j]+uint64_t(a.x[i])*b.x[j])%p;for(int k=2*D-2;k>=D;k--)for(int j=0;j<D;j++)c[k-D+j]=(c[k-D+j]+p-c[k]*modc[j]%p)%p;F r;for(int i=0;i<D;i++)r.x[i]=c[i];return r;}
F inv(const F&a){if(zero(a))throw runtime_error("NONUNIT_DENOMINATOR_"+stage);nmod_poly_t aa,g,ii;nmod_poly_init(aa,p);nmod_poly_init(g,p);nmod_poly_init(ii,p);for(int i=0;i<D;i++){nmod_poly_set_coeff_ui(aa,i,a.x[i]);nmod_poly_set_coeff_ui(g,i,modc[i]);}nmod_poly_set_coeff_ui(g,D,1);int ok=nmod_poly_invmod(ii,aa,g);if(!ok)throw runtime_error("NONUNIT_PROPER_FACTOR_"+stage);F b;for(int i=0;i<D;i++)b.x[i]=nmod_poly_get_coeff_ui(ii,i);nmod_poly_clear(aa);nmod_poly_clear(g);nmod_poly_clear(ii);assert(eq(mul(a,b),F(1)));if(inverse_log.is_open()){uint8_t phase=stage=="rank";inverse_log.write((char*)&phase,1);for(const F*t:{&a,static_cast<const F*>(&b)})for(int i=0;i<D;i++){uint16_t c=t->x[i];inverse_log.write((char*)&c,2);}inverse_checks++;}return b;}
F divf(const F&a,const F&b){return mul(a,inv(b));}F powf(F a,int n){assert(n>=0);F r(1);while(n){if(n&1)r=mul(r,a);n>>=1;if(n)a=mul(a,a);}return r;}
void jf(ostream&o,const F&a){o<<'[';for(int i=0;i<D;i++){if(i)o<<',';o<<a.x[i];}o<<']';}
using V=vector<F>;using Poly=vector<pair<int,F>>;
vector<int>codes,weights;map<int,int>monidx;vector<int>joint239,cubic91,pure98,H1,H2;vector<vector<int>>cubeMon;vector<int>purepos;
int mid(int code){auto it=monidx.find(code);if(it==monidx.end())throw runtime_error("monomial_out_of_range");return it->second;}
void accumulate(map<int,F>&z,int i,const F&t){auto it=z.find(i);if(it==z.end())z.emplace(i,t);else it->second=add(it->second,t);}
Poly finish(map<int,F>&z){Poly r;for(auto&t:z)if(!zero(t.second))r.emplace_back(t.first,move(t.second));return r;}
Poly padd(const Poly&a,const Poly&b){map<int,F>z;for(auto&t:a)accumulate(z,t.first,t.second);for(auto&t:b)accumulate(z,t.first,t.second);return finish(z);}
Poly scale(const F&c,const Poly&a){Poly r;if(zero(c))return r;for(auto&t:a){F q=mul(c,t.second);if(!zero(q))r.emplace_back(t.first,move(q));}return r;}
Poly pmul(const Poly&a,const Poly&b){map<int,F>z;for(auto&t:a)for(auto&u:b){int i=mid(codes[t.first]+codes[u.first]);accumulate(z,i,mul(t.second,u.second));}return finish(z);}
vector<Poly>linear8;
Poly monimage(int code){Poly r;r.emplace_back(mid(0),F(1));for(int j=0;j<8;j++)for(int k=0;k<((code>>(3*j))&7);k++)r=pmul(r,linear8[j]);return r;}
Poly joint(const Poly&a){map<int,F>z;for(auto&t:a){Poly u=monimage(codes[t.first]);for(auto&v:u)accumulate(z,v.first,mul(t.second,v.second));}return finish(z);}
struct Node{char op='z';int x=0,y=0;vector<pair<int,int>>l;};vector<Node>ON,SN,LN;vector<int>oroots,olabels,nfroots;V OA,SC;vector<vector<int>>closure,usecount;
V corepowers[3];
F corepoly(istream&in){int n;in>>n;assert(n>=0&&n<=150000);F r;for(int k=0;k<n;k++){uint32_t c;int ex[3];in>>c>>ex[0]>>ex[1]>>ex[2];F t(c);for(int j=0;j<3;j++){assert(ex[j]>=0);while(int(corepowers[j].size())<=ex[j])corepowers[j].push_back(mul(corepowers[j].back(),corepowers[j][1]));t=mul(t,corepowers[j][ex[j]]);}r=add(r,t);}return r;}
F corerf(istream&in){F a=corepoly(in),b=corepoly(in);return divf(a,b);}F unipoly(istream&in,const F&z){int n;in>>n;F r;for(int k=0;k<n;k++){int e;uint32_t c;in>>e>>c;r=add(r,mul(F(c),powf(z,e)));}return r;}F unirf(istream&in,const F&z){F a=unipoly(in,z),b=unipoly(in,z);return divf(a,b);}
Node readnode(istream&in){Node n;in>>n.op;if(n.op=='m'||n.op=='a'||n.op=='s')in>>n.x>>n.y;else if(n.op=='v'||n.op=='c')in>>n.x;else if(n.op=='l'){int k;in>>k;for(int j=0;j<k;j++){int a,b;in>>a>>b;n.l.emplace_back(a,b);}}else assert(n.op=='z');return n;}
V nf(const Poly&a){V input(239);for(auto&t:a){auto it=find(joint239.begin(),joint239.end(),t.first);assert(it!=joint239.end());input[it-joint239.begin()]=t.second;}V v(LN.size());for(int i=0;i<int(LN.size());i++){auto&n=LN[i];if(n.op=='v')v[i]=input[n.x];else if(n.op=='l')for(auto&c:n.l)v[i]=add(v[i],mul(SC[c.first],v[c.second]));else assert(n.op=='z');}V r;for(int i:nfroots)r.push_back(v[i]);return r;}
Poly oracle(int ri,const vector<F>&Dv){vector<unique_ptr<Poly>>value(ON.size());vector<int>uses=usecount[ri];for(int i:closure[ri]){auto&n=ON[i];Poly q;if(n.op=='v'){q.emplace_back(mid(1<<(3*n.x)),Dv[n.x]);}else if(n.op=='m')q=pmul(*value[n.x],*value[n.y]);else if(n.op=='l'){map<int,F>z;for(auto&c:n.l)for(auto&t:*value[c.second])accumulate(z,t.first,mul(OA[c.first],t.second));q=finish(z);}else assert(n.op=='z');value[i]=make_unique<Poly>(move(q));vector<int>children;if(n.op=='m')children={n.x,n.y};else if(n.op=='l')for(auto&c:n.l)children.push_back(c.second);for(int j:children)if(!--uses[j])value[j].reset();}return move(*value[oroots[ri]]);}
void prepare_closure(){closure.clear();usecount.clear();for(int r:oroots){set<int>need;vector<int>todo={r};vector<int>use(ON.size());while(todo.size()){int i=todo.back();todo.pop_back();if(!need.insert(i).second)continue;auto&n=ON[i];vector<int>children;if(n.op=='m')children={n.x,n.y};else if(n.op=='l')for(auto&c:n.l)children.push_back(c.second);for(int j:children){use[j]++;todo.push_back(j);}}use[r]++;closure.emplace_back(need.begin(),need.end());usecount.push_back(move(use));}}
struct Src{int type=0,index=0,mon=0;};vector<Src>catalogue;
V gamma(const V&f,const V&chi,const F&U,const vector<V>&section){F h=divf(f[H2[0]],U);V out;for(int k:pure98)out.push_back(sub(f[k],mul(h,chi[k])));for(int j=0;j<23;j++){F c=sub(f[H1[j]],mul(h,chi[H1[j]]));if(!zero(c))for(int k=0;k<98;k++)out[k]=sub(out[k],mul(c,section[j][k]));}return out;}
void writeGamma(ofstream&o,int row,const V&v){uint16_t r=row;o.write((char*)&r,2);for(auto&a:v)for(int k=0;k<D;k++){uint16_t c=a.x[k];o.write((char*)&c,2);}}
void phase(int f,const string&s,int count=0){ofstream o(folder+"/native-phase.json");o<<"{\"factor_index\":"<<f<<",\"stage\":\""<<s<<"\",\"rows\":"<<count<<",\"live_slots\":"<<live*D<<",\"peak_slots\":"<<peak<<"}\n";}
void coeff(int fi,const string&shared,const string&cfile,const string&gfile){cap=107743;stage="coefficient";ifstream in(shared);int k,n;in>>n;assert(n==int(codes.size()));for(int i=0;i<n;i++){int a,b;in>>a>>b;assert(a==codes[i]&&b==weights[i]);}
 F z;if(D==1)z=F(modc[0]?p-modc[0]:0);else z.x[1]=1;vector<pair<int,int>>poles={{1,0},{1,20},{1,60},{203,2700},{1,12}};for(auto[a,b]:poles){F t=add(mul(F(a),z),F(b));assert(!zero(t));(void)inv(t);}
 F a=divf(mul(F(p-10),add(z,F(20))),mul(F(9),mul(z,z))),e=divf(add(z,F(100)),mul(F(9),z)),d=mul(F(p-3),divf(mul(z,add(z,F(20))),F(100)));
 corepowers[0]={F(1),a};corepowers[1]={F(1),e};corepowers[2]={F(1),d};F C=sub(sub(add(mul(F(100),a),mul(F(18),mul(e,e))),mul(F(9),e)),F(2));F B=sub(add(mul(sub(sub(add(mul(F(50),a),mul(F(18),mul(e,e))),mul(F(9),e)),F(2)),d),mul(F(60),e)),F(40));F delta=mul(F(p-64000%p),mul(C,C));assert(!zero(delta));F DH=divf(B,F(50));vector<F>Dv={a,a,a,F(1),a,DH,a,a};F invDH2=inv(mul(DH,DH));
 in>>n;OA.clear();OA.reserve(n);for(int i=0;i<n;i++)OA.push_back(corerf(in));in>>n;ON.clear();for(int i=0;i<n;i++)ON.push_back(readnode(in));in>>n;oroots.clear();olabels.clear();for(int i=0;i<n;i++){int l,r;in>>l>>r;olabels.push_back(l);oroots.push_back(r);}prepare_closure();
 in>>n;V nfatom;for(int i=0;i<n;i++)nfatom.push_back(corerf(in));in>>n;SN.clear();SC.clear();for(int i=0;i<n;i++){Node q=readnode(in);SN.push_back(q);F t;if(q.op=='c')t=nfatom[q.x];else if(q.op=='m')t=mul(SC[q.x],SC[q.y]);else if(q.op=='a')t=add(SC[q.x],SC[q.y]);else if(q.op=='s')t=sub(SC[q.x],SC[q.y]);else assert(q.op=='z');SC.push_back(move(t));}nfatom.clear();
 in>>n;LN.clear();for(int i=0;i<n;i++)LN.push_back(readnode(in));in>>n;assert(n==122);nfroots.resize(n);for(auto&x:nfroots)in>>x;
 linear8.clear();vector<int>light={0,1,2,3,4,6,7};for(int i=0;i<8;i++){Poly q;if(i==5)q.emplace_back(mid(1<<(3*7)),F(1));else{int row=find(light.begin(),light.end(),i)-light.begin();for(int j=0;j<7;j++){F t=unirf(in,z);if(!zero(t))q.emplace_back(mid(1<<(3*j)),move(t));}}linear8.push_back(move(q));}
 F u=corerf(in),v=corerf(in),U=corerf(in);assert(!zero(U));assert(eq(U,mul(powf(mul(mul(F(50),a),mul(B,C)),2),F(1))));
 vector<Poly>quad(6);for(int i=0;i<6;i++){in>>n;Poly raw;for(int j=0;j<n;j++){int code;in>>code;F t=corepoly(in);if(!zero(t))raw.emplace_back(mid(code),move(t));}quad[i]=joint(raw);}
 int active;in>>active;vector<int>lr(active);vector<V>L(active,V(23));for(int i=0;i<active;i++){in>>lr[i];for(auto&t:L[i])t=unirf(in,z);}string extra;assert(!(in>>extra));
 // H2 source first, with only one oracle polynomial root resident at a time.
 phase(fi,"chi_source");int r77=find(olabels.begin(),olabels.end(),77)-olabels.begin(),r131=find(olabels.begin(),olabels.end(),131)-olabels.begin();V chi(122);
 {Poly raw=scale(invDH2,oracle(r77,Dv));V n=nf(joint(raw));for(int j=0;j<122;j++)chi[j]=mul(u,n[j]);}
 {Poly raw=scale(invDH2,oracle(r131,Dv));V n=nf(joint(raw));for(int j=0;j<122;j++)chi[j]=add(chi[j],mul(v,n[j]));}assert(eq(chi[H2[0]],U));
 // All16 full cubics come from the streamed original normalized57x91 source.
 phase(fi,"cubic_original16");ifstream cf(cfile);int cc;cf>>cc;assert(cc==16);vector<Poly>cub(16);for(int i=0;i<16;i++){Poly raw;for(int code:cubic91){F t=corepoly(cf);if(!zero(t))raw.emplace_back(mid(code),move(t));}cub[i]=joint(raw);}assert(!(cf>>extra));
 vector<V>f112;f112.reserve(112);for(int i=0;i<16;i++)for(int l:light)f112.push_back(nf(pmul(cub[i],linear8[l])));cub.clear();
 vector<V>section(23,V(98));vector<vector<F>>identity(23,V(23));for(int rr=0;rr<active;rr++){int r=lr[rr];V original=f112[r];if(r/7==6){for(int k=0;k<122;k++)original[k]=mul(delta,sub(f112[r][k],f112[7+r%7][k]));}for(int j=0;j<23;j++)if(!zero(L[rr][j])){for(int k=0;k<98;k++)section[j][k]=add(section[j][k],mul(L[rr][j],original[pure98[k]]));for(int k=0;k<23;k++)identity[j][k]=add(identity[j][k],mul(L[rr][j],original[H1[k]]));}}
 for(int i=0;i<23;i++)for(int j=0;j<23;j++)assert(eq(identity[i][j],F(i==j)));identity.clear();L.clear();
 ofstream go(gfile,ios::binary);int count=0;for(int r=0;r<112;r++){assert(zero(f112[r][H2[0]]));writeGamma(go,174+r,gamma(f112[r],chi,U,section));count++;}f112.clear();phase(fi,"112_cubics_written_and_lift_cache_released",count);
 for(int r=0;r<174;r++){auto&s=catalogue[r];assert(s.type==0);Poly q=pmul(quad[s.index],monimage(s.mon));V ff=nf(q);writeGamma(go,r,gamma(ff,chi,U,section));count++;}quad.clear();
 for(int j=0;j<33;j++){phase(fi,"one_raw_quartic",count);Poly raw=scale(invDH2,oracle(j,Dv));V ff=nf(joint(raw));writeGamma(go,286+j,gamma(ff,chi,U,section));count++;}assert(count==319);go.close();phase(fi,"all319_written",count);
 OA.clear();SC.clear();linear8.clear();for(auto&x:corepowers)x.clear();
}
void rankfactor(int fi,const string&gfile,const string&outfile){assert(live==0);cap=81712;stage="rank";phase(fi,"rank_after_all_coefficient_cache_release");ifstream in(gfile,ios::binary);vector<V>basis(98),selected;vector<int>pivots,rows;V pivotvalues;int rnk=0;for(int n=0;n<319;n++){uint16_t row;in.read((char*)&row,2);assert(in);V orig(98);for(auto&t:orig)for(int k=0;k<D;k++){uint16_t c;in.read((char*)&c,2);assert(in&&c<p);t.x[k]=c;}V r=orig;for(int j=0;j<98;j++)if(basis[j].size()&&!zero(r[j])){F c=r[j];for(int k=j;k<98;k++)r[k]=sub(r[k],mul(c,basis[j][k]));}int j=0;while(j<98&&zero(r[j]))j++;if(j==98)continue;F pv=r[j],ip=inv(pv);for(auto&t:r)t=mul(t,ip);basis[j]=move(r);selected.push_back(move(orig));rows.push_back(row);pivots.push_back(j);pivotvalues.push_back(pv);rnk++;if(rnk==98)break;}
 ofstream o(outfile);o<<"{\"factor_index\":"<<fi<<",\"degree\":"<<D<<",\"rank\":"<<rnk<<",\"selected_source319_rows\":[";for(int i=0;i<rnk;i++){if(i)o<<',';o<<rows[i];}o<<"],\"pivot_columns\":[";for(int i=0;i<rnk;i++){if(i)o<<',';o<<pivots[i];}o<<"],\"pivot_values\":[";for(int i=0;i<rnk;i++){if(i)o<<',';jf(o,pivotvalues[i]);}F det(1);int swaps=0;for(int i=0;i<rnk;i++){det=mul(det,pivotvalues[i]);for(int j=i+1;j<rnk;j++)swaps+=pivots[i]>pivots[j];}if(swaps%2)det=neg(det);o<<"],\"selected98_minor_determinant\":";jf(o,det);o<<",\"selected_original_rows\":[";for(int i=0;i<rnk;i++){if(i)o<<',';o<<'[';for(int j=0;j<98;j++){if(j)o<<',';jf(o,selected[i][j]);}o<<']';}o<<"],\"rank_phase_peak_slots\":"<<peak<<",\"all319_source_coefficients_computed_before_rank\":true,\"positive_full98\":"<<(rnk==98?"true":"false")<<"}\n";if(rnk!=98)throw runtime_error("INCONCLUSIVE_DEFICIENT_FACTOR_"+to_string(fi));assert(!zero(det));}
int main(int argc,char**argv){assert(argc==5);began=chrono::steady_clock::now();limit=stod(argv[4]);folder=argv[3];try{
 ifstream cfg(string(argv[1])+"/config.txt");int n;cfg>>n;for(int i=0;i<n;i++){int c,w;cfg>>c>>w;monidx[c]=codes.size();codes.push_back(c);weights.push_back(w);}cfg>>n;for(int i=0;i<n;i++){int c;cfg>>c;joint239.push_back(mid(c));}cfg>>n;for(int i=0;i<n;i++){int c;cfg>>c;cubic91.push_back(c);}for(auto*arr:{&pure98,&H1,&H2}){cfg>>n;arr->resize(n);for(auto&x:*arr)cfg>>x;}cfg>>n;assert(n==319);for(int i=0;i<n;i++){Src s;cfg>>s.type>>s.index>>s.mon;catalogue.push_back(s);}cfg>>n;vector<vector<uint32_t>>factors(n);for(auto&f:factors){int d;cfg>>d;f.resize(d+1);for(auto&c:f)cfg>>c;assert(f.back()==1);}cfg.close();int done=0;long globalpeak=0;for(int i=0;i<n;i++){assert(live==0);D=factors[i].size()-1;for(int j=0;j<D;j++)modc[j]=factors[i][j];peak=0;string gf=string(argv[2])+"/Gamma-factor"+to_string(i)+".bin";inverse_log.open(string(argv[2])+"/inverse-factor"+to_string(i)+".bin",ios::binary);{coeff(i,string(argv[1])+"/shared.txt",string(argv[1])+"/cubics.txt",gf);}assert(live==0);long cp=peak;peak=0;rankfactor(i,gf,folder+"/factor"+to_string(i)+".json");assert(live==0);inverse_log.close();globalpeak=max({globalpeak,cp,peak});done++;ofstream sum(folder+"/native-summary.json");sum<<"{\"all_passed\":"<<(done==n?"true":"false")<<",\"completed_full98_factors\":"<<done<<",\"factor_count\":"<<n<<",\"max_live_prime_field_slots\":"<<globalpeak<<",\"all_inverse_products_checked_and_binary_logged\":"<<inverse_checks<<",\"no_new_prime_or_source_point\":true,\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";}
 }catch(const exception&e){ofstream o(folder+"/frontier.json");o<<"{\"all_passed\":false,\"stage\":\""<<stage<<"\",\"reason\":\""<<e.what()<<"\",\"live_slots\":"<<live*D<<",\"peak_slots\":"<<peak<<",\"no_auto_retry\":true}\n";cerr<<e.what()<<'\n';return 2;}flint_cleanup();}
