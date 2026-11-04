// Independent inverse-coordinate audit of the frozen L-frame input.
// No rank, determinant, forward-coordinate elimination or producer imports.
#include <flint/fmpq_mpoly.h>
#include <array>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
using namespace std;
fmpq_mpoly_ctx_t qc;
struct Poly {
 fmpq_mpoly_t x;
 Poly(){fmpq_mpoly_init(x,qc);} Poly(const Poly&b):Poly(){fmpq_mpoly_set(x,b.x,qc);}
 Poly(Poly&&b):Poly(){fmpq_mpoly_swap(x,b.x,qc);} Poly&operator=(Poly b){fmpq_mpoly_swap(x,b.x,qc);return *this;}
 ~Poly(){fmpq_mpoly_clear(x,qc);}
};
using Exp=array<int,5>;using Light=map<Exp,Poly>;using Matrix=vector<vector<Poly>>;
Poly constant(long a){Poly p;fmpq_mpoly_set_si(p.x,a,qc);return p;}
Poly multiply(const Poly&a,const Poly&b){Poly p;fmpq_mpoly_mul(p.x,a.x,b.x,qc);return p;}
Poly scale(const Poly&a,long n){Poly p;fmpq_mpoly_scalar_mul_si(p.x,a.x,n,qc);return p;}
void add(Poly&a,const Poly&b){fmpq_mpoly_add(a.x,a.x,b.x,qc);}
void accum(Light&a,const Exp&e,const Poly&p){add(a[e],p);}
bool same(const Poly&a,const Poly&b){return fmpq_mpoly_equal(a.x,b.x,qc);}
bool zero(const Poly&a){return fmpq_mpoly_is_zero(a.x,qc);}
void load(istream&f,Poly&p){long n;f>>n;assert(f&&n>=0);fmpq_t c;fmpq_init(c);ulong e[3];
 for(long i=0;i<n;i++){string v;f>>v>>e[0]>>e[1]>>e[2];assert(f&&!fmpq_set_str(c,v.c_str(),10));fmpq_canonicalise(c);fmpq_mpoly_push_term_fmpq_ui(p.x,c,e,qc);}
 fmpq_mpoly_sort_terms(p.x,qc);fmpq_mpoly_combine_like_terms(p.x,qc);fmpq_clear(c);
}
Matrix matrix(const string&p,int nr,int nc,int extra=-1){ifstream f(p);int r,c,k;f>>r>>c;assert(f&&r==nr&&c==nc);if(extra>=0){f>>k;assert(k==extra);}Matrix a(r,vector<Poly>(c));for(auto&row:a)for(auto&v:row)load(f,v);string tail;assert(!(f>>tail));return a;}
vector<Exp> monomials(int degree){vector<Exp>v;Exp e{};auto rec=[&](auto&&self,int j,int left)->void{
 if(j==4){e[j]=left;v.push_back(e);return;}for(int k=left;k>=0;k--){e[j]=k;self(self,j+1,left-k);}};rec(rec,0,degree);return v;
}
Light product(const Light&a,const Light&b){Light z;for(auto&[e,c]:a)for(auto&[f,d]:b){Exp g;for(int j=0;j<5;j++)g[j]=e[j]+f[j];accum(z,g,multiply(c,d));}return z;}
void equalmat(const Matrix&a,const Matrix&b){assert(a.size()==b.size());for(int i=0;i<(int)a.size();i++){assert(a[i].size()==b[i].size());for(int j=0;j<(int)a[i].size();j++)assert(same(a[i][j],b[i][j]));}}
long count_terms(const Matrix&a){long n=0;for(auto&r:a)for(auto&p:r)n+=fmpq_mpoly_length(p.x,qc);return n;}
long max_degree(const Matrix&a){long n=-1;for(auto&r:a)for(auto&p:r)if(!zero(p))n=max(n,fmpq_mpoly_total_degree_si(p.x,qc));return n;}
int main(int argc,char**argv){assert(argc==6);fmpq_mpoly_ctx_init(qc,3,ORD_LEX);{
 auto start=chrono::steady_clock::now();string in=argv[1],out=argv[5];ifstream meta(in+"/independent-meta.txt");
 int n;meta>>n;assert(n==91);vector<vector<pair<int,long>>> restriction(n);
 for(auto&row:restriction){int k;meta>>k;while(k--){int j;long c;meta>>j>>c;row.push_back({j,c});}}
 vector<int>strip(16);int qstrip;meta>>qstrip;assert(qstrip==2);for(int&i:strip){meta>>i;assert(i>=0&&i<=3);}assert(meta);
 auto m2=monomials(2),m3=monomials(3);map<Exp,int>i2,i3;for(int i=0;i<15;i++)i2[m2[i]]=i;for(int i=0;i<35;i++)i3[m3[i]]=i;
 Matrix oldQ=matrix(in+"/original-qfour4x15.txt",4,15);
 Matrix Q4=matrix(out+"/qfour4x15.txt",4,15),G40=matrix(out+"/genuine16x40.txt",16,40);
 Matrix Q20=matrix(argv[3],20,35),Gold(16,vector<Poly>(40));
 // Rebuild all original quad products from independently restricted source.
 for(int r=0;r<4;r++)for(int v=0;v<5;v++){vector<Poly>p(35);for(int j=0;j<15;j++){Exp e=m2[j];e[v]++;p[i3.at(e)]=oldQ[r][j];}
  for(int j=0;j<35;j++)assert(same(Q20[5*r+v][j],p[j]));}
 // Read the original actual57 source; restriction is independently rebuilt
 // from its public8-variable monomial metadata by the Python reviewer.
 ifstream raw(argv[2]);int r,c,k;raw>>r>>c>>k;assert(r==57&&c==91&&k==57);
 for(int i=0;i<57;i++)for(int j=0;j<91;j++){Poly p;load(raw,p);if(i>=41)for(auto[col,m]:restriction[j])add(Gold[i-41][col],scale(p,m));}
 string tail;assert(!(raw>>tail));
 Matrix L=matrix(out+"/L-coefficients.txt",5,1);Poly C,a,epar;fmpq_mpoly_gen(a.x,0,qc);fmpq_mpoly_gen(epar.x,1,qc);
 C=scale(a,100);add(C,scale(multiply(epar,epar),18));add(C,scale(epar,-9));add(C,constant(-2));Poly den=scale(C,-40);assert(same(L[1][0],den));
 Light linear;for(int j=0;j<5;j++){Exp e{};e[j]=1;linear[e]=L[j][0];}
 // The original source identities q1=L²/8100 and
 // q3=L(eL-400a*x0)/8100, BEFORE row rescaling.
 Light lsq=product(linear,linear),lx;for(auto&[e,p]:linear){Exp f=e;f[0]++;lx[f]=p;}
 for(int j=0;j<15;j++){
  assert(same(scale(oldQ[1][j],8100),lsq[m2[j]]));
  Poly expected=multiply(epar,lsq[m2[j]]);add(expected,scale(multiply(a,lx[m2[j]]),-400));assert(same(scale(oldQ[2][j],8100),expected));
 }
 // Recover L coefficients from the actual q1 derivative; no forward producer helper.
 for(int v=0;v<5;v++){Exp ex{};ex[1]++;ex[v]++;Poly derivative=scale(oldQ[1][i2.at(ex)],v==1?2:1);
  assert(same(scale(multiply(C,L[v][0]),4),scale(derivative,-405)));}
 // Substitute L=sum_j L_j*old_xj BACK into every new row. No polynomial
 // division is performed: multiply by the saved C powers on the left.
 vector<Light>powers(4);powers[0][Exp{}]=constant(1);for(int j=1;j<=3;j++)powers[j]=product(powers[j-1],linear);
 auto back=[&](const Exp&ne){Light z;Exp base{};base[0]=ne[1];base[2]=ne[2];base[3]=ne[3];base[4]=ne[4];
  for(auto&[e,p]:powers[ne[0]]){Exp g;for(int j=0;j<5;j++)g[j]=base[j]+e[j];accum(z,g,p);}return z;};
 vector<Light>back2,back3,back1;for(auto&x:m2)back2.push_back(back(x));for(auto&x:m3)back3.push_back(back(x));for(int j=0;j<5;j++){Exp e{};e[j]=1;back1.push_back(back(e));}
 vector<Poly>dp(4),cp(4);dp[0]=cp[0]=constant(1);for(int j=1;j<=3;j++){dp[j]=multiply(dp[j-1],den);cp[j]=multiply(cp[j-1],C);}
 Matrix Qback(4,vector<Poly>(15));for(int r=0;r<4;r++)for(int j=0;j<15;j++)for(auto&[e,p]:back2[j])add(Qback[r][i2.at(e)],multiply(p,Q4[r][j]));
 for(int r=0;r<4;r++)for(int j=0;j<15;j++){
  if(r==0){assert(same(Qback[r][j],oldQ[r][j]));if(m2[j][1])assert(zero(oldQ[r][j]));}
  else assert(same(multiply(Qback[r][j],cp[r==3?qstrip:0]),multiply(oldQ[r][j],dp[2])));
 }
 for(int r=0;r<16;r++){
  vector<Poly>rec(40);for(int j=0;j<35;j++)for(auto&[e,p]:back3[j])add(rec[i3.at(e)],multiply(p,G40[r][j]));
  for(int j=0;j<5;j++)for(auto&[e,p]:back1[j]){int v=-1;for(int h=0;h<5;h++)if(e[h]){assert(e[h]==1&&v<0);v=h;}assert(v>=0);add(rec[35+v],multiply(p,G40[r][35+j]));}
  for(int j=0;j<40;j++)assert(same(multiply(rec[j],cp[strip[r]]),multiply(Gold[r][j],dp[3])));
 }
 // The normalized L² and L*x0 products span exactly the first nine
 // cubic monomial columns; this is a fixed unit deletion, not rank search.
 vector<int>covered(35);for(int v=0;v<5;v++){Exp e{};e[0]=2;e[v]++;covered[i3.at(e)]=1;}for(int v=1;v<5;v++){Exp e{};e[0]=e[1]=1;e[v]++;covered[i3.at(e)]=1;}
 for(int j=0;j<35;j++)assert(covered[j]==int(j<9));
 vector<int>keep;for(int j=15;j<35;j++)keep.push_back(j);for(int j=9;j<15;j++)keep.push_back(j);
 Matrix Q10(10,vector<Poly>(26));int row=0;for(int q:{0,3})for(int v=1;v<5;v++){
  for(int j=0;j<15;j++){Exp ex=m2[j];ex[v]++;for(int col=0;col<26;col++)if(m3[keep[col]]==ex)add(Q10[row][col],Q4[q][j]);}row++;}
 for(int q:{0,3}){for(int j=0;j<15;j++){Exp ex=m2[j];ex[0]++;for(int col=0;col<26;col++)if(m3[keep[col]]==ex)add(Q10[row][col],Q4[q][j]);}row++;}assert(row==10);
 equalmat(Q10,matrix(out+"/quad10x26.txt",10,26));
 Matrix A(8,vector<Poly>(20)),B(8,vector<Poly>(6)),Pair(2,vector<Poly>(6)),G31(16,vector<Poly>(31));
 for(int i=0;i<8;i++){for(int j=0;j<20;j++)A[i][j]=Q10[i][j];for(int j=0;j<6;j++)B[i][j]=Q10[i][20+j];}
 for(int i=0;i<2;i++){for(int j=0;j<20;j++)assert(zero(Q10[8+i][j]));for(int j=0;j<6;j++)Pair[i][j]=Q10[8+i][20+j];}
 for(int i=0;i<16;i++){for(int j=0;j<26;j++)G31[i][j]=G40[i][keep[j]];for(int j=0;j<5;j++)G31[i][26+j]=G40[i][35+j];}
 equalmat(A,matrix(out+"/quad8x20.txt",8,20));equalmat(B,matrix(out+"/quad8x6.txt",8,6));equalmat(Pair,matrix(out+"/pair2x6.txt",2,6));equalmat(G31,matrix(out+"/genuine16x31.txt",16,31));
 Matrix frozenPair=matrix(argv[4],2,6,0);
 for(int i=0;i<6;i++){
  assert(same(Pair[0][i],frozenPair[0][5-i]));
  assert(same(multiply(multiply(C,cp[qstrip]),Pair[1][i]),multiply(dp[2],frozenPair[1][5-i])));
 }
 assert(count_terms(Q10)==317&&count_terms(Pair)==40&&count_terms(G31)==50085&&max_degree(G31)==19);
 cout<<"{\"all_passed\":true,\"original_qfour_and_16genuine_source_rebuilt\":true,\"all60_quad_and640_genuine_inverse_polynomial_identities\":true,\"all700_original_quadproduct_entries\":true,\"L_derivative_and_q1_q3_identities\":true,\"nine_monomial_unit_deletions\":true,\"all10x26_and16x31_entries\":true,\"all12_pair_coefficients\":true,\"quad10_terms\":"<<count_terms(Q10)<<",\"genuine16x31_terms\":"<<count_terms(G31)<<",\"pair_terms\":"<<count_terms(Pair)<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}"<<endl;
 }fmpq_mpoly_ctx_clear(qc);flint_cleanup();return 0;}
