// Same exact integer recursion as integer_implicit_jets_20261101.cpp,
// in a homogeneous monomial quotient: multiplication-map entry -1 means zero.
// After the F rows, append the observable coordinate44 row in the same scaling.
// Exact implicit jets for quadratic systems with an invertible linear block.
// Input: F,n,D,r; den,td; degree sizes; pivot rows; integer J,T,U; quadratic terms; monomial product maps.
// U/den is a selected right inverse of J, and T/td is the tangent frame.
// Each output row contains rational coefficients in degree order 2..D.
// See S61_FLAG_BOUNDARY_EXACT_JETS_20261101.md for the normalization and validation.
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>
#include <gmpxx.h>
using namespace std;
using Z=mpz_class;using Poly=vector<pair<int,Z>>;
int main(int argc,char**argv){
 if(argc!=3)return 1;
 auto start=chrono::steady_clock::now();ifstream in(argv[1]);int F,n,D,r;in>>F>>n>>D>>r;Z den,td;in>>den>>td;
 vector<int>sizes(D+1),rows(r);for(int&i:sizes)in>>i;for(int&i:rows)in>>i;
 vector<Poly>J(F),U(45);vector<vector<Poly>>W(D+1,vector<Poly>(45));
 for(int i=0;i<F;i++)for(int j=0;j<45;j++){Z a;in>>a;if(a!=0)J[i].emplace_back(j,a);}
 for(int i=0;i<45;i++)for(int j=0;j<n;j++){Z a;in>>a;if(a!=0)W[1][i].emplace_back(j,a);}
 for(int i=0;i<45;i++)for(int j=0;j<r;j++){Z a;in>>a;if(a!=0)U[i].emplace_back(j,a);}
 map<pair<int,int>,vector<pair<int,Z>>>H;
 for(int i=0;i<F;i++){int count;in>>count;for(int k=0;k<count;k++){int a,b;Z c;in>>a>>b>>c;if(a>b)swap(a,b);H[{a,b}].emplace_back(i,c);}}
 map<pair<int,int>,vector<int>>maps;
 for(int a=1;a<D;a++)for(int b=1;a+b<=D;b++){auto&v=maps[{a,b}];v.resize(sizes[a]*sizes[b]);for(int&i:v)in>>i;}
 if(!in){cerr<<"bad input\n";return 2;}
 vector<vector<vector<Z>>>A(D+1);
 for(int d=2;d<=D;d++){
  vector<vector<Z>>Q(F,vector<Z>(sizes[d]));
  for(auto&[ij,entries]:H){
   int i=ij.first,j=ij.second;vector<Z>product(sizes[d]);
   for(int a=1;a<d;a++){
    int b=d-a;auto&mp=maps[{a,b}];
    for(auto&[k,x]:W[a][i])for(auto&[l,y]:W[b][j])
     if(mp[k*sizes[b]+l]>=0)mpz_addmul(product[mp[k*sizes[b]+l]].get_mpz_t(),x.get_mpz_t(),y.get_mpz_t());
   }
   for(auto&[row,c]:entries)for(int k=0;k<sizes[d];k++)if(product[k]!=0)
    mpz_addmul(Q[row][k].get_mpz_t(),c.get_mpz_t(),product[k].get_mpz_t());
  }
  for(int i=0;i<45;i++){
   vector<Z>v(sizes[d]);for(auto&[j,c]:U[i])for(int k=0;k<sizes[d];k++)if(Q[rows[j]][k]!=0)
    mpz_submul(v[k].get_mpz_t(),c.get_mpz_t(),Q[rows[j]][k].get_mpz_t());
   for(int k=0;k<sizes[d];k++)if(v[k]!=0)W[d][i].emplace_back(k,std::move(v[k]));
  }
  for(int i=0;i<F;i++){
   for(auto&c:Q[i])c*=den;
   for(auto&[j,c]:J[i])for(auto&[k,v]:W[d][j])mpz_addmul(Q[i][k].get_mpz_t(),c.get_mpz_t(),v.get_mpz_t());
  }
  A[d]=std::move(Q);
  for(int row:rows)for(auto&v:A[d][row])if(v!=0){cerr<<"pivot residual not zero\n";return 3;}
  cerr<<"integer degree "<<d<<" ready seconds "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"\n";
 }
 vector<Z>scale(D+1);for(int d=2;d<=D;d++){Z a,b;mpz_pow_ui(a.get_mpz_t(),den.get_mpz_t(),d-1);mpz_pow_ui(b.get_mpz_t(),td.get_mpz_t(),d);scale[d]=a*b;}
 ofstream out(argv[2]);for(int i=0;i<F;i++){
  for(int d=2;d<=D;d++)for(auto&a:A[d][i]){mpq_class q(a,scale[d]);q.canonicalize();out<<q.get_num()<<"/"<<q.get_den()<<" ";}out<<"\n";
 }
 for(int d=2;d<=D;d++){vector<Z>v(sizes[d]);for(auto&[k,c]:W[d][44])v[k]=c;
  for(auto&a:v){mpq_class q(a,scale[d]);q.canonicalize();out<<q.get_num()<<"/"<<q.get_den()<<" ";}}out<<"\n";
 cerr<<"EXACT INTEGER JETS PASS seconds "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"\n";
}
