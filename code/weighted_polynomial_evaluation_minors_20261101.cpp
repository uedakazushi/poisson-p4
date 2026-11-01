#include <flint/nmod_mat.h>
#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
using U=unsigned long;
struct Term {U c; std::array<int,10> e;};
int main(int argc,char**argv){
 if(argc!=3)return 2; std::ifstream in(argv[1]);std::ofstream out(argv[2]);
 U p,seed;int np,npoints,ng;in>>p>>seed>>np>>npoints>>ng;if(!in || !out || p!=32003 || np!=45 || npoints<1 || ng<1)return 6;
 std::vector<std::vector<Term>> polys(np);
 for(auto &f:polys){int n;in>>n;f.resize(n);for(auto&t:f){in>>t.c;for(int&v:t.e)in>>v;}}
 std::mt19937_64 rng(seed);std::vector<std::array<U,10>> points(npoints);
 std::vector<std::vector<U>> values(npoints,std::vector<U>(np));
 for(int i=0;i<npoints;i++){
  for(auto&v:points[i])v=rng()%p;
  for(int j=0;j<np;j++)for(auto&t:polys[j]){
   U a=t.c;for(int k=0;k<10;k++)for(int l=0;l<t.e[k];l++)a=a*points[i][k]%p;
   values[i][j]=(values[i][j]+a)%p;
  }
 }
 auto st=std::chrono::steady_clock::now();long total=0;
 out<<"{\"prime\":"<<p<<",\"seed\":"<<seed<<",\"points\":[";
 for(int i=0;i<npoints;i++){if(i)out<<',';out<<'[';for(int j=0;j<10;j++){if(j)out<<',';out<<points[i][j];}out<<']';}out<<"],\"groups\":[";
 for(int g=0;g<ng;g++){
  int a,b,r,n;in>>a>>b>>r>>n;std::vector<std::array<int,4>> products(n);
  for(auto&v:products)for(int&x:v)in>>x;
  nmod_mat_t E,B,S;nmod_mat_init(E,r,n,p);nmod_mat_init(B,r,n,p);
  for(int i=0;i<r;i++)for(int j=0;j<n;j++){
   U co=1;for(int k:products[j])co=co*values[i][k]%p;nmod_mat_entry(E,i,j)=co;
  }
  nmod_mat_set(B,E);long rank=nmod_mat_rref(B);
  std::cerr<<"weight "<<a<<' '<<b<<" rank "<<rank<<'/'<<r<<" cols "<<n<<std::endl;
  if(rank!=r)return 3;
  std::vector<int> piv;
  for(int i=0;i<r;i++){int j=0;while(j<n && !nmod_mat_entry(B,i,j))j++;if(j==n)return 4;piv.push_back(j);}
  nmod_mat_init(S,r,r,p);for(int i=0;i<r;i++)for(int j=0;j<r;j++)nmod_mat_entry(S,i,j)=nmod_mat_entry(E,i,piv[j]);
  U det=nmod_mat_det(S);if(!det)return 5;
  if(g)out<<',';out<<"{\"weight\":["<<a<<','<<b<<"],\"rank\":"<<r<<",\"columns\":"<<n<<",\"determinant\":"<<det<<",\"selected_products\":[";
  for(int j=0;j<r;j++){if(j)out<<',';out<<'[';for(int k=0;k<4;k++){if(k)out<<',';out<<products[piv[j]][k];}out<<']';}out<<"]}";
  total+=r;nmod_mat_clear(E);nmod_mat_clear(B);nmod_mat_clear(S);
 }
 double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-st).count();
 out<<"],\"total_rank\":"<<total<<",\"seconds\":"<<seconds<<"}\n";
 std::cerr<<"ALL EVALUATION MINORS PASS "<<total<<" seconds "<<seconds<<std::endl;
}
