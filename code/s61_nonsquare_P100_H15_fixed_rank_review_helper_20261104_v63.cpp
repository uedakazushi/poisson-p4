// Review helper: fixed selected square matrices over the seven frozen fields.
// No coefficient generation, source queries, pivot search, factor or prime choice.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
constexpr long p=32003;using E=array<long,4>;
long d;array<long,4> f;long field_muls=0;
E add(E a,const E&b){for(long i=0;i<d;i++)a[i]=(a[i]+b[i])%p;return a;}
E sub(E a,const E&b){for(long i=0;i<d;i++)a[i]=(a[i]+p-b[i])%p;return a;}
bool zero(const E&a){for(long i=0;i<d;i++)if(a[i])return false;return true;}
E one(){E a{};a[0]=1;return a;}
E mul(const E&a,const E&b){
 field_muls++;array<long,7> c{};
 for(long i=0;i<d;i++)for(long j=0;j<d;j++)c[i+j]=(c[i+j]+a[i]*b[j])%p;
 for(long k=2*d-2;k>=d;k--)for(long i=0;i<d;i++)c[k-d+i]=(c[k-d+i]+p-(c[k]*f[i])%p)%p;
 E r{};copy_n(c.begin(),d,r.begin());return r;
}
E inv(E a){assert(!zero(a));uint64_t q=1;for(long i=0;i<d;i++)q*=p;uint64_t e=q-2;E r=one(),original=a;
 while(e){if(e&1)r=mul(r,a);e>>=1;if(e)a=mul(a,a);}assert(mul(original,r)==one());return r;
}
E readE(istream&in){E a{};for(long i=0;i<d;i++){in>>a[i];assert(in&&0<=a[i]&&a[i]<p);}return a;}
void writeE(ostream&out,const E&a){out<<'[';for(long i=0;i<d;i++){if(i)out<<',';out<<a[i];}out<<']';}
int main(int argc,char**argv){assert(argc==3);auto start=chrono::steady_clock::now();string folder=argv[2];try{
 ifstream in(argv[1]);long nf,prime;in>>nf>>prime;assert(nf==7&&prime==p);ofstream out(folder+"/fixed-rank-review-native.json");
 out<<"{\"independent_field_array_arithmetic\":true,\"no_pivot_search\":true,\"factors\":[";
 for(long factor=0;factor<nf;factor++){
  long fi,n;in>>fi>>d>>n;assert(fi==factor&&1<=d&&d<=4&&n==98);
  f={};for(long i=0;i<d;i++){in>>f[i];assert(in&&0<=f[i]&&f[i]<p);}long leading;in>>leading;assert(leading==1);
  vector<long> source(n),cols(n);array<bool,319> seen{};array<bool,98> used{};
  for(long&i:source){in>>i;assert(0<=i&&i<319&&!seen[i]);seen[i]=true;}
  for(long&i:cols){in>>i;assert(0<=i&&i<98&&!used[i]);used[i]=true;}
  vector<vector<E>> A(n,vector<E>(n)),basis;for(auto&row:A)for(auto&x:row)x=readE(in);
  assert(2*n*n*d+2*n*d+256<150000);vector<E> pivots;E determinant=one();
  for(long i=0;i<n;i++){
   vector<E> row=A[i];
   for(long j=0;j<i;j++){E c=row[cols[j]];if(!zero(c))for(long k=0;k<n;k++)row[k]=sub(row[k],mul(c,basis[j][k]));}
   for(long j=0;j<i;j++)assert(zero(row[cols[j]]));E pivot=row[cols[i]];assert(!zero(pivot));pivots.push_back(pivot);determinant=mul(determinant,pivot);E ip=inv(pivot);
   for(E&x:row)x=mul(x,ip);assert(row[cols[i]]==one());basis.push_back(move(row));
   if(chrono::duration<double>(chrono::steady_clock::now()-start).count()>49)throw runtime_error("49s_fixed_verifier_guard");
  }
  long parity=0;for(long i=0;i<n;i++)for(long j=0;j<i;j++)parity^=cols[j]>cols[i];if(parity)determinant=sub(E{},determinant);assert(!zero(determinant));
  if(factor)out<<',';out<<"{\"factor_index\":"<<factor<<",\"degree\":"<<d<<",\"fixed98_pivots_all_nonzero\":true,\"determinant_in_original_column_order\":";writeE(out,determinant);out<<",\"fixed_pivot_values\":[";
  for(long i=0;i<n;i++){if(i)out<<',';writeE(out,pivots[i]);}out<<"]}";
 }
 string extra;assert(!(in>>extra));out<<"],\"all_passed\":true,\"CRT_full_Hbar_row_surjectivity_gate\":true,\"field_multiplications\":"<<field_muls<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
 }catch(const exception&e){ofstream out(folder+"/frontier.json");out<<"{\"all_passed\":false,\"reason\":\""<<e.what()<<"\"}\n";return 2;}}
