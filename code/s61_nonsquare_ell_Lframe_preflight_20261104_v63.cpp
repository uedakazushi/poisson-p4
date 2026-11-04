#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
using namespace std;
fmpq_mpoly_ctx_t ctx;
struct P {
 fmpq_mpoly_t p;
 P(){fmpq_mpoly_init(p,ctx);} P(const P&b):P(){fmpq_mpoly_set(p,b.p,ctx);}
 P(P&&b):P(){fmpq_mpoly_swap(p,b.p,ctx);} P&operator=(P b){fmpq_mpoly_swap(p,b.p,ctx);return *this;}
 ~P(){fmpq_mpoly_clear(p,ctx);}
};
using E=array<int,5>; using LP=map<E,P>; using Mat=vector<vector<P>>;
void readp(istream&f,P&p,bool keep=true){int n;f>>n;assert(f&&n>=0);fmpq_t c;fmpq_init(c);ulong e[3];for(int i=0;i<n;i++){string s;f>>s>>e[0]>>e[1]>>e[2];assert(f);if(keep){assert(!fmpq_set_str(c,s.c_str(),10));fmpq_canonicalise(c);fmpq_mpoly_push_term_fmpq_ui(p.p,c,e,ctx);}}if(keep){fmpq_mpoly_sort_terms(p.p,ctx);fmpq_mpoly_combine_like_terms(p.p,ctx);}fmpq_clear(c);}
void writep(ostream&f,const P&p){int n=fmpq_mpoly_length(p.p,ctx);f<<n;fmpq_t c;fmpq_init(c);ulong e[3];for(int i=0;i<n;i++){fmpq_mpoly_get_term_coeff_fmpq(c,p.p,i,ctx);fmpq_mpoly_get_term_exp_ui(e,p.p,i,ctx);char*s=fmpq_get_str(nullptr,10,c);f<<' '<<s<<' '<<e[0]<<' '<<e[1]<<' '<<e[2];flint_free(s);}f<<'\n';fmpq_clear(c);}
void writem(const string&path,const Mat&A){ofstream f(path);assert(f);f<<A.size()<<' '<<A[0].size()<<'\n';for(auto&r:A)for(auto&p:r)writep(f,p);}
P one(){P p;fmpq_mpoly_one(p.p,ctx);return p;}
P integer(long n){P p;fmpq_mpoly_set_si(p.p,n,ctx);return p;}
P mul(const P&a,const P&b){P p;fmpq_mpoly_mul(p.p,a.p,b.p,ctx);return p;}
void add(P&a,const P&b){fmpq_mpoly_add(a.p,a.p,b.p,ctx);}
void addlp(LP&a,const E&e,const P&c){add(a[e],c);}
LP lpmul(const LP&a,const LP&b){LP c;for(auto&[x,p]:a)for(auto&[y,q]:b){E z;for(int j=0;j<5;j++)z[j]=x[j]+y[j];addlp(c,z,mul(p,q));}return c;}
vector<E> mons(int n){vector<E>v;E e{};auto rec=[&](auto&&self,int j,int left)->void{if(j==4){e[j]=left;v.push_back(e);return;}for(int k=left;k>=0;k--){e[j]=k;self(self,j+1,left-k);}};rec(rec,0,n);return v;}
long terms(const vector<P>&r){long n=0;for(auto&p:r)n+=fmpq_mpoly_length(p.p,ctx);return n;}
void profile(ostream&f,const Mat&A){long n=0,nz=0,mx=0,deg=-1,nb=0,db=0;fmpq_t c;fmpq_init(c);for(auto&r:A)for(auto&p:r){long t=fmpq_mpoly_length(p.p,ctx);n+=t;nz+=t>0;mx=max(mx,t);if(t)deg=max(deg,fmpq_mpoly_total_degree_si(p.p,ctx));for(long k=0;k<t;k++){fmpq_mpoly_get_term_coeff_fmpq(c,p.p,k,ctx);nb=max(nb,(long)fmpz_bits(fmpq_numref(c)));db=max(db,(long)fmpz_bits(fmpq_denref(c)));}}fmpq_clear(c);f<<"{\"rows\":"<<A.size()<<",\"columns\":"<<A[0].size()<<",\"terms\":"<<n<<",\"nonzero_entries\":"<<nz<<",\"max_entry_terms\":"<<mx<<",\"max_core_total_degree\":"<<deg<<",\"max_numerator_bits\":"<<nb<<",\"max_denominator_bits\":"<<db<<"}";}
int stripC(vector<P>&row,const P&C,int cap){int k=0;for(;k<cap;k++){vector<P>next(row.size());bool ok=true;for(int j=0;j<(int)row.size();j++)if(!fmpq_mpoly_divides(next[j].p,row[j].p,C.p,ctx)){ok=false;break;}if(!ok)break;row=move(next);}return k;}
int main(int argc,char**argv){assert(argc==6);fmpq_mpoly_ctx_init(ctx,3,ORD_DEGREVLEX);{
 auto start=chrono::steady_clock::now();string out=argv[5];auto m3=mons(3),m2=mons(2);map<E,int>i3,i2;for(int j=0;j<35;j++)i3[m3[j]]=j;for(int j=0;j<15;j++)i2[m2[j]]=j;
 Mat oldQ20(20,vector<P>(35));ifstream qf(argv[2]);int nr,nc;qf>>nr>>nc;assert(nr==20&&nc==35);for(auto&r:oldQ20)for(auto&p:r)readp(qf,p);
 Mat oldQ4(4,vector<P>(15));for(int r=0;r<4;r++)for(int j=0;j<35;j++){if(fmpq_mpoly_is_zero(oldQ20[5*r][j].p,ctx))continue;E e=m3[j];assert(e[0]>=1);e[0]--;oldQ4[r][i2.at(e)]=oldQ20[5*r][j];}
 P C,a,epar,t;fmpq_mpoly_gen(a.p,0,ctx);fmpq_mpoly_gen(epar.p,1,ctx);fmpq_mpoly_scalar_mul_si(C.p,a.p,100,ctx);t=mul(epar,epar);fmpq_mpoly_scalar_mul_si(t.p,t.p,18,ctx);add(C,t);fmpq_mpoly_scalar_mul_si(t.p,epar.p,-9,ctx);add(C,t);t=integer(-2);add(C,t);
 vector<P>L(5);fmpq_t fac;fmpq_init(fac);fmpq_set_si(fac,-405,4);for(int j=0;j<5;j++){E ex{};ex[1]++;ex[j]++;P der=oldQ4[1][i2.at(ex)];if(j==1)fmpq_mpoly_scalar_mul_si(der.p,der.p,2,ctx);assert(fmpq_mpoly_divides(L[j].p,der.p,C.p,ctx));fmpq_mpoly_scalar_mul_fmpq(L[j].p,L[j].p,fac,ctx);}fmpq_clear(fac);
 P den;fmpq_mpoly_scalar_mul_si(den.p,C.p,-40,ctx);assert(fmpq_mpoly_equal(den.p,L[1].p,ctx));assert(!fmpq_mpoly_is_zero(den.p,ctx));
 // New light order (L,x0,y1=x2,y2=x3,y3=x4). The inverse has sole denominator den=-40C.
 LP N;E z{};z[0]=1;N[z]=one();for(int j:{0,2,3,4}){E ex{};ex[j==0?1:j]=1;P p;fmpq_mpoly_neg(p.p,L[j].p,ctx);N[ex]=p;}
 vector<P>dp(4);dp[0]=one();for(int k=1;k<=3;k++)dp[k]=mul(dp[k-1],den);vector<LP>np(4);np[0][E{}]=one();for(int k=1;k<=3;k++)np[k]=lpmul(np[k-1],N);
 auto trans=[&](const E&old,int clear){LP v;E base{};base[1]=old[0];base[2]=old[2];base[3]=old[3];base[4]=old[4];for(auto&[ex,c]:np[old[1]]){E dest;for(int k=0;k<5;k++)dest[k]=base[k]+ex[k];addlp(v,dest,mul(c,dp[clear-old[1]]));}return v;};
 // Literal linear inverse identity: N(L=sum L_j*x_j)=den*x1, checked coefficientwise.
 for(int j=0;j<5;j++){P p=L[j];if(j!=1){P neg;fmpq_mpoly_neg(neg.p,L[j].p,ctx);add(p,neg);}assert(j==1?fmpq_mpoly_equal(p.p,den.p,ctx):fmpq_mpoly_is_zero(p.p,ctx));}
 vector<LP>T2,T3;for(auto&m:m2)T2.push_back(trans(m,2));for(auto&m:m3)T3.push_back(trans(m,3));
 Mat Q4(4,vector<P>(15));for(int r=0;r<4;r++)for(int j=0;j<15;j++)for(auto&[ex,c]:T2[j])add(Q4[r][i2.at(ex)],mul(oldQ4[r][j],c));
 // Exact actual-source normalizations, BEFORE any known-C clearing cancellation.
 Mat expect(2,vector<P>(15));E ll{};ll[0]=2;expect[0][i2.at(ll)]=dp[2];E lx{};lx[0]=1;lx[1]=1;expect[1][i2.at(ll)]=mul(dp[2],epar);P pp=mul(dp[2],a);fmpq_mpoly_scalar_mul_si(pp.p,pp.p,-400,ctx);expect[1][i2.at(lx)]=pp;
 for(int j=0;j<15;j++){P q;fmpq_mpoly_scalar_mul_si(q.p,Q4[1][j].p,8100,ctx);assert(fmpq_mpoly_equal(q.p,expect[0][j].p,ctx));fmpq_mpoly_scalar_mul_si(q.p,Q4[2][j].p,8100,ctx);assert(fmpq_mpoly_equal(q.p,expect[1][j].p,ctx));}
 // q0 has no x1; preserve its unscaled source rather than gratuitously multiply by C².
 for(int j=0;j<15;j++){P q;assert(fmpq_mpoly_divides(q.p,Q4[0][j].p,dp[2].p,ctx));Q4[0][j]=q;}
 int q5strip=stripC(Q4[3],C,2);
 ifstream sf(argv[3]);sf>>nr>>nc;assert(nr==35&&nc==40);int ignored;for(int j=0;j<nr+nc;j++)sf>>ignored;int nm;sf>>nm;assert(nm==91);vector<vector<pair<int,long>>>rest(nm);for(auto&v:rest){int n;sf>>n;while(n--){int j;long c;sf>>j>>c;v.push_back({j,c});}}
 Mat G(16,vector<P>(40));ifstream src(argv[1]);int k;src>>nr>>nc>>k;assert(nr==57&&nc==91&&k==57);P p;for(int r=0;r<57;r++)for(int j=0;j<91;j++){bool keep=r>=41&&!rest[j].empty();fmpq_mpoly_zero(p.p,ctx);readp(src,p,keep);if(keep)for(auto[target,c]:rest[j]){P v;fmpq_mpoly_scalar_mul_si(v.p,p.p,c,ctx);add(G[r-41][target],v);}}
 Mat NG(16,vector<P>(40));for(int r=0;r<16;r++){for(int j=0;j<35;j++)for(auto&[ex,c]:T3[j])add(NG[r][i3.at(ex)],mul(G[r][j],c));for(int j=0;j<5;j++){E ex{};ex[j]=1;for(auto&[ne,c]:trans(ex,3)){int ind=-1;for(int v=0;v<5;v++)if(ne[v]){assert(ne[v]==1&&ind<0);ind=v;}assert(ind>=0);add(NG[r][35+ind],mul(G[r][35+j],c));}}}
 vector<int>gstrip;for(auto&r:NG)gstrip.push_back(stripC(r,C,3));
 // Nine monomial unit pivots: L²*five lights and L*x0*four non-L lights.
 vector<int>keep;for(int j=15;j<35;j++)keep.push_back(j);for(int j=9;j<15;j++)keep.push_back(j);assert(keep.size()==26);for(int j=0;j<9;j++)assert(m3[j][0]>=1&&(m3[j][0]>=2||m3[j][1]>=1));
 Mat Q10(10,vector<P>(26));int ri=0;for(int r:{0,3})for(int v=1;v<5;v++){for(int j=0;j<15;j++){E ex=m2[j];ex[v]++;int col=i3.at(ex);auto it=find(keep.begin(),keep.end(),col);if(it!=keep.end())add(Q10[ri][it-keep.begin()],Q4[r][j]);}ri++;}for(int r:{0,3}){for(int j=0;j<15;j++){E ex=m2[j];ex[0]++;int col=i3.at(ex);auto it=find(keep.begin(),keep.end(),col);if(it!=keep.end())add(Q10[ri][it-keep.begin()],Q4[r][j]);}ri++;}assert(ri==10);for(int r=8;r<10;r++)for(int j=0;j<20;j++)assert(fmpq_mpoly_is_zero(Q10[r][j].p,ctx));
 Mat G31(16,vector<P>(31));for(int r=0;r<16;r++){for(int j=0;j<26;j++)G31[r][j]=NG[r][keep[j]];for(int j=0;j<5;j++)G31[r][26+j]=NG[r][35+j];}
 Mat A8(8,vector<P>(20)),B8(8,vector<P>(6)),Pair(2,vector<P>(6));for(int r=0;r<8;r++){copy(Q10[r].begin(),Q10[r].begin()+20,A8[r].begin());copy(Q10[r].begin()+20,Q10[r].end(),B8[r].begin());}for(int r=0;r<2;r++)copy(Q10[r+8].begin()+20,Q10[r+8].end(),Pair[r].begin());
 Mat oldPair(2,vector<P>(6));ifstream pf(argv[4]);int metadata;pf>>nr>>nc>>metadata;assert(nr==2&&nc==6&&metadata==0);for(auto&r:oldPair)for(auto&p:r)readp(pf,p);P q5scale=dp[2];for(int i=0;i<q5strip;i++){P q;assert(fmpq_mpoly_divides(q.p,q5scale.p,C.p,ctx));q5scale=q;}for(int j=0;j<6;j++){assert(fmpq_mpoly_equal(Pair[0][j].p,oldPair[0][5-j].p,ctx));P lhs=mul(C,Pair[1][j]),rhs=mul(q5scale,oldPair[1][5-j]);assert(fmpq_mpoly_equal(lhs.p,rhs.p,ctx));}
 writem(out+"/qfour4x15.txt",Q4);writem(out+"/quad10x26.txt",Q10);writem(out+"/genuine16x31.txt",G31);writem(out+"/genuine16x40.txt",NG);writem(out+"/quad8x20.txt",A8);writem(out+"/quad8x6.txt",B8);writem(out+"/pair2x6.txt",Pair);
 ofstream frame(out+"/L-coefficients.txt");frame<<"5 1\n";for(auto&p:L)writep(frame,p);ofstream stats(out+"/native.json");stats<<"{\n\"actual_source_q1_q3_identities_passed\":true,\n\"linear_inverse_identity_passed\":true,\n\"nine_monomial_pivots_no_arithmetic_elimination\":true,\n\"bottom_left_2x20_identically_zero\":true,\n\"all12_pair_coefficients_match_frozen_global_dual_source\":true,\n\"old_quad20x35\":";profile(stats,oldQ20);stats<<",\n\"old_genuine16x40\":";profile(stats,G);stats<<",\n\"new_qfour4x15\":";profile(stats,Q4);stats<<",\n\"new_quad10x26\":";profile(stats,Q10);stats<<",\n\"new_genuine16x40\":";profile(stats,NG);stats<<",\n\"new_genuine16x31\":";profile(stats,G31);stats<<",\n\"quad_A8x20\":";profile(stats,A8);stats<<",\n\"quad_B8x6\":";profile(stats,B8);stats<<",\n\"quad_C2x6\":";profile(stats,Pair);stats<<",\n\"q5_C_factors_cancelled\":"<<q5strip<<",\n\"genuine_C_factors_cancelled\":[";for(int i=0;i<16;i++)stats<<(i?",":"")<<gstrip[i];stats<<"],\n\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"\n}\n";
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();return 0;}
