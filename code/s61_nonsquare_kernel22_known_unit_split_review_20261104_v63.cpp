// Independent fixed-step replay. No pivot candidates are searched.
#include "s61_binary_hyperplane_cubic_qq3_20261102.hpp"
#include <sstream>
using Matrix=vector<vector<Rat>>;
Rat getrat(istream&in){Rat c;readpoly(in,c.n);readpoly(in,c.d);assert(!fmpq_mpoly_is_zero(c.d,ctx));return c;}
Rat one(){Rat x;fmpq_mpoly_one(x.n,ctx);return x;}
bool equal(const Rat&a,const Rat&b){return sub(a,b).zero();}
int main(int argc,char**argv){
 assert(argc==4);fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);{
  ifstream input(argv[1]),state(argv[2]),log(argv[3]);int n,m,r,nn,mm;
  input>>nn>>mm;state>>n>>m>>r;assert(nn==n&&mm==m&&n==174&&m==22&&r==20);
  vector<int>rows(r),cols(r);for(int&i:rows)state>>i;for(int&j:cols)state>>j;
  string status;int verified;double seconds;long numterms;state>>status>>verified>>seconds>>numterms;
  assert(status=="no_individual_known_unit_entry"&&verified==1);
  Matrix O(n,vector<Rat>(m));for(auto&row:O)for(auto&c:row)c=getrat(input);
  Matrix A(r,vector<Rat>(r)),L=A,U=A,S(n,vector<Rat>(m));
  for(auto&row:A)for(auto&c:row)c=getrat(state);for(auto&row:L)for(auto&c:row)c=getrat(state);
  for(auto&row:U)for(auto&c:row)c=getrat(state);for(auto&row:S)for(auto&c:row)c=getrat(state);
  string extra;assert(!(input>>extra)&&!(state>>extra));
  const char*names[]={"a","e","d"};vector<Rat>known;
  for(const char*s:{"a","d","500*a+90*e^2+3*e-24","50*a*d+18*d*e^2-9*d*e-2*d+60*e-40","100*a+18*e^2-9*e-2","5*a*d-6*e-4"}){
   Rat v;assert(!fmpq_mpoly_set_str_pretty(v.n,s,names,ctx));known.push_back(v);
  }
  Matrix B=O,low(n,vector<Rat>(r));vector<bool>usedrow(n),usedcol(m);string line;
  for(int k=0;k<r;k++){
   getline(log,line);istringstream t(line);string a,b,c,d,e,f,g,h,scalar;int kk,ii,jj;long terms;double sec;
   t>>a>>kk>>b>>ii>>c>>jj>>d>>terms>>e>>sec>>f>>scalar>>g;
   assert(a=="PIVOT"&&b=="ROW"&&c=="COL"&&d=="TERMS"&&e=="SEC"&&f=="SCALAR"&&g=="EXPONENTS");
   assert(kk==k&&ii==rows[k]&&jj==cols[k]&&!usedrow[ii]&&!usedcol[jj]);
   Rat unit;assert(!fmpq_mpoly_set_str_pretty(unit.n,scalar.c_str(),names,ctx));
   for(int j=0;j<6;j++){int p;t>>p;for(int l=0;l<abs(p);l++)unit=p>0?mul(unit,known[j]):divi(unit,known[j]);}
   assert(!(t>>extra));Rat pivot=B[ii][jj];assert(equal(pivot,unit));
   assert(equal(pivot,U[k][k]));
   // Eliminate only the SPECIFIED pivot; all unselected rows/columns retained.
   for(int i=0;i<n;i++)if(!usedrow[i]&&i!=ii&&!B[i][jj].zero()){
    Rat multiplier=divi(B[i][jj],pivot);low[i][k]=multiplier;
    for(int j=0;j<m;j++)if(!usedcol[j]&&j!=jj)B[i][j]=sub(B[i][j],mul(multiplier,B[ii][j]));
    B[i][jj]=Rat();
   }
   usedrow[ii]=usedcol[jj]=true;
  }
  for(int i=0;i<n;i++)for(int j=0;j<m;j++)assert(equal(B[i][j],S[i][j]));
  for(int i=0;i<r;i++)for(int j=0;j<r;j++){
   assert(equal(A[i][j],O[rows[i]][cols[j]]));
   Rat expectedL=i==j?one():(i>j?low[rows[i]][j]:Rat());assert(equal(expectedL,L[i][j]));
   Rat expectedU=i<=j?B[rows[i]][cols[j]]:Rat();assert(equal(expectedU,U[i][j]));
   Rat sum;for(int k=0;k<r;k++)sum=sub(sum,sub(Rat(),mul(L[i][k],U[k][j])));
   assert(equal(sum,A[i][j]));
  }
  cout<<"FIXED20_UNIT_LOG_ALL3828_SCHUR_ENTRIES_AND_A_L_U_PASS\n";
  cout<<"NO_ADDITIONAL_PIVOT_MINOR_OR_BEZOUT_SEARCH\n";
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();
}
