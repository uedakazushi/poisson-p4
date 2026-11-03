// Independent integer Gaussian determinant of the pinned670x670 minor.
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
static constexpr uint32_t p_review=32003;
uint32_t power(uint32_t a,uint32_t n){uint32_t r=1;while(n){if(n&1)r=(uint64_t)r*a%p_review;a=(uint64_t)a*a%p_review;n>>=1;}return r;}
int main(int argc,char**argv){
 assert(argc==2);ifstream f(argv[1]);int n;uint32_t p;f>>n>>p;assert(n==670&&p==p_review);
 vector<uint32_t>A((size_t)n*n);for(auto&c:A){f>>c;assert(f&&c<p_review);}string extra;assert(!(f>>extra));
 uint32_t determinant=1;
 for(int k=0;k<n;k++){
  int pivot=k;while(pivot<n&&!A[(size_t)pivot*n+k])pivot++;assert(pivot<n);
  if(pivot!=k){for(int j=k;j<n;j++)swap(A[(size_t)pivot*n+j],A[(size_t)k*n+j]);determinant=determinant?p_review-determinant:0;}
  uint32_t diagonal=A[(size_t)k*n+k];determinant=(uint64_t)determinant*diagonal%p_review;uint32_t inverse=power(diagonal,p_review-2);
  for(int i=k+1;i<n;i++)if(A[(size_t)i*n+k]){
   uint32_t multiplier=(uint64_t)A[(size_t)i*n+k]*inverse%p_review;A[(size_t)i*n+k]=0;
   for(int j=k+1;j<n;j++){
    uint32_t term=(uint64_t)multiplier*A[(size_t)k*n+j]%p_review;uint32_t&entry=A[(size_t)i*n+j];entry=entry>=term?entry-term:entry+p_review-term;
   }
  }
 }
 assert(determinant);cout<<"INDEPENDENT_INTEGER_MINOR670_DETERMINANT "<<determinant<<'\n';
}
