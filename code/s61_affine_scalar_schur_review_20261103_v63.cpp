// Independent integer modular determinant of the original224x224 coefficient minor.
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
constexpr uint32_t review_p=32003;
uint32_t inv(uint32_t a){uint32_t r=1,n=review_p-2;while(n){if(n&1)r=(uint64_t)r*a%review_p;a=(uint64_t)a*a%review_p;n>>=1;}return r;}
int main(int argc,char**argv){
 assert(argc==2);ifstream f(argv[1]);int n;uint32_t p;f>>n>>p;assert(n==224&&p==review_p);
 vector<uint32_t>A((size_t)n*n);for(auto&x:A){f>>x;assert(f&&x<p);}string extra;assert(!(f>>extra));uint32_t det=1;
 for(int k=0;k<n;k++){
  int q=k;while(q<n&&!A[(size_t)q*n+k])q++;assert(q<n);
  if(q!=k){for(int j=k;j<n;j++)swap(A[(size_t)q*n+j],A[(size_t)k*n+j]);det=p-det;}
  uint32_t v=A[(size_t)k*n+k];det=(uint64_t)det*v%p;uint32_t vi=inv(v);
  for(int i=k+1;i<n;i++)if(A[(size_t)i*n+k]){
   uint32_t a=(uint64_t)A[(size_t)i*n+k]*vi%p;A[(size_t)i*n+k]=0;
   for(int j=k+1;j<n;j++){uint32_t b=(uint64_t)a*A[(size_t)k*n+j]%p;auto&x=A[(size_t)i*n+j];x=x>=b?x-b:x+p-b;}
  }
 }
 assert(det);cout<<"INDEPENDENT_INTEGER224_DETERMINANT "<<det<<'\n';
}
