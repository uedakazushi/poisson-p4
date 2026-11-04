// Independent scalar-DAG evaluation, COLUMN-BY-COLUMN linear evaluation,
// and direct physical-quad dual transports. No search or quartic matrix.
#include "s61_binary_hyperplane_cubic_qq3_20261102.hpp"
#include <array>
using Mat=vector<vector<Rat>>;
Rat rd(istream&f){Rat x;readpoly(f,x.n);readpoly(f,x.d);return x;}
Rat one(){Rat x;fmpq_mpoly_one(x.n,ctx);return x;}
Rat number(long x){Rat a;fmpq_mpoly_set_si(a.n,x,ctx);return a;}
Rat sumrat(const Rat&a,const Rat&b){return sub(a,sub(Rat(),b));}
bool eq(const Rat&a,const Rat&b){return sub(a,b).zero();}
vector<vector<int>> mons(int n,int d){if(n==1)return{{d}};vector<vector<int>>v;for(int i=0;i<=d;i++)for(auto z:mons(n-1,d-i)){z.insert(z.begin(),i);v.push_back(z);}return v;}
vector<Rat> restriction(const vector<Rat>&Q,const Mat&T){
 int n=T[0].size();auto old=mons(7,2),dest=mons(n,2);map<vector<int>,int>idx;for(int i=0;i<(int)dest.size();i++)idx[dest[i]]=i;
 vector<Rat>v(dest.size());for(int h=0;h<28;h++)if(!Q[h].zero()){
  vector<int>is;for(int i=0;i<7;i++)for(int k=0;k<old[h][i];k++)is.push_back(i);assert(is.size()==2);
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(!T[is[0]][i].zero()&&!T[is[1]][j].zero()){
   vector<int>e(n);e[i]++;e[j]++;int k=idx[e];v[k]=sumrat(v[k],mul(Q[h],mul(T[is[0]][i],T[is[1]][j])));
  }
 }return v;
}
int main(int argc,char**argv){assert(argc==2);fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);{
 ifstream f(argv[1]);int na,ns,nl,nc;f>>na>>ns>>nl>>nc;assert(na==378&&ns==2491&&nl==1685&&nc==45);
 vector<Rat>atoms(na);const char*names[]={"a","e","d"};vector<Rat>u;
 for(const char*s:{"a","d","500*a+90*e^2+3*e-24","50*a*d+18*d*e^2-9*d*e-2*d+60*e-40","100*a+18*e^2-9*e-2","5*a*d-6*e-4"}){Rat x;assert(!fmpq_mpoly_set_str_pretty(x.n,s,names,ctx));u.push_back(x);}
 fmpq_mpoly_t den,tmp;fmpq_mpoly_init(den,ctx);fmpq_mpoly_init(tmp,ctx);
 for(auto&x:atoms){x=rd(f);assert(!fmpq_mpoly_is_zero(x.d,ctx));fmpq_mpoly_set(den,x.d,ctx);for(auto&v:u)while(fmpq_mpoly_divides(tmp,den,v.n,ctx))fmpq_mpoly_swap(den,tmp,ctx);assert(fmpq_mpoly_is_fmpq(den,ctx));}
 fmpq_mpoly_clear(den,ctx);fmpq_mpoly_clear(tmp,ctx);
 vector<Rat>s(ns);for(int i=0;i<ns;i++){char t;f>>t;if(t=='z')continue;int a,b;f>>a;if(t=='c'){assert(a<na);s[i]=atoms[a];continue;}f>>b;assert(a<i&&b<i);s[i]=t=='a'?sumrat(s[a],s[b]):t=='s'?sub(s[a],s[b]):mul(s[a],s[b]);assert(t=='a'||t=='s'||t=='m');}
 struct Node{char t;int v=0;vector<pair<int,int>>terms;};vector<Node>ln(nl);
 for(int i=0;i<nl;i++){auto&n=ln[i];f>>n.t;if(n.t=='v'){f>>n.v;assert(n.v>=0&&n.v<174);}else if(n.t=='l'){int k;f>>k;while(k--){int c,j;f>>c>>j;assert(c<ns&&j<i);n.terms.push_back({c,j});}}else assert(n.t=='z');}
 for(int i=0;i<nc;i++){int root,want;f>>root>>want;assert(eq(s[root],number(want)));}
 vector<int>outputs(22);for(int&i:outputs)f>>i;Mat K(174,vector<Rat>(22));for(auto&r:K)for(auto&x:r)x=rd(f);
 // Independent evaluation order: scalar operator on one frozen column at a time.
 for(int j=0;j<22;j++){vector<Rat>v(nl);for(int i=0;i<nl;i++){auto&n=ln[i];if(n.t=='v')v[i]=K[n.v][j];else for(auto[c,k]:n.terms)if(!s[c].zero()&&!v[k].zero())v[i]=sumrat(v[i],mul(s[c],v[k]));}
  for(int i=0;i<22;i++)assert(eq(v[outputs[i]],number(i==j)));
 }
 cout<<"ALL484_LITERAL_P22_K_EQUALS_I22_PASS\n";
 Mat Q(6,vector<Rat>(28));for(auto&r:Q)for(auto&x:r)x=rd(f);vector<Rat>L(7),ell(7);for(auto&x:L)x=rd(f);for(auto&x:ell)x=rd(f);
 Mat A(6,vector<Rat>(6)),inv=A;for(auto&r:A)for(auto&x:r)x=rd(f);for(auto&r:inv)for(auto&x:r)x=rd(f);
 for(int i=0;i<6;i++)for(int j=0;j<6;j++){Rat z;for(int k=0;k<6;k++)z=sumrat(z,mul(A[i][k],inv[k][j]));assert(eq(z,number(i==j)));}
 assert(eq(L[1],mul(number(-40),u[4]))&&eq(ell[1],divi(mul(number(-16),mul(u[4],u[5])),number(243))));
 Rat g=divi(mul(number(-80),mul(u[0],mul(u[1],u[1]))),number(81));
 vector<Rat>D(5);D[1]=divi(sub(mul(number(-3),L[2]),mul(number(-7),divi(L[4],number(2)))),L[1]);D[2]=number(3);D[4]=divi(number(-7),number(2));
 Rat dl,de;for(int i=0;i<5;i++){dl=sumrat(dl,mul(D[i],L[i]));de=sumrat(de,mul(D[i],ell[i]));}assert(dl.zero()&&eq(de,g));
 Mat pairM(2,vector<Rat>(6));for(auto&r:pairM)for(auto&x:r)x=rd(f);vector<Rat>fs(5),bs(5);for(auto&x:fs)x=rd(f);for(auto&x:bs)x=rd(f);
 Rat bz;for(int i=0;i<5;i++)bz=sumrat(bz,mul(fs[i],bs[i]));assert(eq(bz,one()));
 Mat pairT(7,vector<Rat>(3));for(int i=0;i<3;i++){pairT[i+2][i]=one();pairT[1][i]=divi(sub(Rat(),L[i+2]),L[1]);}
 Mat pairq(2);pairq[0]=restriction(Q[0],pairT);pairq[1]=restriction(Q[3],pairT);
 for(int i=0;i<2;i++)for(int j=0;j<6;j++)assert(eq(mul(pairq[i][j],i?u[4]:one()),pairM[i][j]));
 assert(eq(pairM[0][2],mul(number(5),u[0])));vector<int>rem={0,1,3,4,5};Rat alpha=sub(Rat(),mul(u[0],u[4]));
 for(int k=0;k<5;k++){int j=rem[k];assert(eq(sub(pairM[1][j],divi(mul(pairM[1][2],pairM[0][j]),pairM[0][2])),mul(alpha,fs[k])));}
 Mat raw4(21,vector<Rat>(4)),saved4(4,vector<Rat>(21));for(auto&r:raw4)for(auto&x:r)x=rd(f);for(auto&r:saved4)for(auto&x:r)x=rd(f);
 Mat T(7,vector<Rat>(6));vector<int>ix={0,2,3,4,5,6};for(int j=0;j<6;j++){T[ix[j]][j]=one();T[1][j]=divi(sub(Rat(),ell[ix[j]]),ell[1]);}
 Mat qbar(4);for(int i=0;i<4;i++)qbar[i]=restriction(Q[i],T);
 vector<Rat>mult={one(),divi(mul(u[5],u[5]),number(25)),divi(mul(u[5],u[5]),number(25)),divi(u[5],number(5))};
 vector<int>sh={0,2,2,1};Rat g4=sub(Rat(),mul(mul(mul(u[0],u[0]),mul(u[0],u[0])),mul(mul(u[1],u[1]),u[5])));
 for(int i=0;i<4;i++)for(int j=0;j<21;j++)assert(eq(mul(qbar[i][j],mult[i]),saved4[i][j]));
 Mat V=raw4;for(int j=0;j<21;j++)for(int i=0;i<4;i++){Rat aa=one();for(int k=0;k<sh[i];k++)aa=mul(aa,u[0]);V[j][i]=divi(mul(raw4[j][i],mul(aa,mult[i])),g4);}
 for(int i=0;i<4;i++)for(int j=0;j<4;j++){Rat z;for(int k=0;k<21;k++)z=sumrat(z,mul(qbar[i][k],V[k][j]));assert(eq(z,number(i==j)));}
 vector<int>lam(15);for(int&i:lam)f>>i;vector<vector<int>>pd(6,vector<int>(2)),qd(6,vector<int>(28));for(auto&r:pd)for(int&i:r)f>>i;for(auto&r:qd)for(int&i:r)f>>i;
 for(int i=0;i<2;i++)for(int j=0;j<2;j++){Rat z;for(int k=0;k<6;k++)z=sumrat(z,mul(pairq[i][k],s[pd[k][j]]));assert(eq(z,number(i==j)));}
 for(int i=0;i<6;i++)for(int j=0;j<6;j++){Rat z;for(int k=0;k<28;k++)z=sumrat(z,mul(Q[i][k],s[qd[j][k]]));assert(eq(z,number(i==j)));}
 auto m5=mons(5,2),m7=mons(7,2);map<vector<int>,int>at7;for(int i=0;i<28;i++)at7[m7[i]]=i;
 vector<Rat>beta(15);for(int k=0;k<15;k++){auto&m=m5[k];Rat z;for(int j=0;j<5;j++)if(m[1]&&(m[j]-(j==1)))z=sumrat(z,mul(number(m[1]*(m[j]-(j==1))),D[j]));beta[k]=divi(z,mul(g,L[1]));}
 vector<Rat>bq(2);for(int t=0;t<2;t++)for(int k=0;k<15;k++){auto ex=m5[k];ex.push_back(0);ex.push_back(0);bq[t]=sumrat(bq[t],mul(beta[k],Q[t?3:0][at7[ex]]));}
 for(int k=0;k<15;k++){
  vector<Rat>basis(28);auto ex=m5[k];ex.push_back(0);ex.push_back(0);basis[at7[ex]]=one();auto cut=restriction(basis,pairT);Rat value=beta[k];
  for(int t=0;t<2;t++){Rat z;for(int h=0;h<6;h++)z=sumrat(z,mul(cut[h],s[pd[h][t]]));value=sub(value,mul(bq[t],z));}assert(eq(value,s[lam[k]]));
 }
 for(int i=0;i<4;i++){Rat z;for(int k=0;k<15;k++){auto ex=m5[k];ex.push_back(0);ex.push_back(0);z=sumrat(z,mul(s[lam[k]],Q[i][at7[ex]]));}assert(z.zero());}
 Rat ellL;for(int k=0;k<15;k++){vector<int>is;for(int j=0;j<5;j++)for(int z=0;z<m5[k][j];z++)is.push_back(j);Rat c=mul(ell[is[0]],L[is[1]]);if(is[0]!=is[1])c=sumrat(c,mul(ell[is[1]],L[is[0]]));ellL=sumrat(ellL,mul(c,s[lam[k]]));}assert(eq(ellL,one()));
 string extra;assert(!(f>>extra));cout<<"ACTUAL_PAIR12_SCHUR_BEZOUT116_AND_84_ELL_SOURCE_TRANSPORT_PASS\nACTUAL_4DUAL16_6DUAL36_LAMBDA15_AND_ANCHORS_PASS\nALL_KNOWN_FULLU_DENOMINATORS_AND_LINEARITY_PASS\n";
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();
}
