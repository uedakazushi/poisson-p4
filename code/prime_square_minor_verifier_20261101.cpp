// Independent elementary Gaussian determinant verification over F_32003.
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
constexpr std::uint32_t P=32003;
std::uint32_t power(std::uint32_t a,std::uint32_t n){
 std::uint32_t r=1;while(n){if(n&1)r=r*a%P;a=a*a%P;n>>=1;}return r;
}
int main(int argc,char**argv){
 if(argc!=2)return 2;std::ifstream in(argv[1],std::ios::binary);
 auto read=[&](){std::uint32_t v=0;in.read(reinterpret_cast<char*>(&v),4);return v;};
 if(read()!=P)return 3;auto count=read();std::uint64_t total=0;
 for(std::uint32_t g=0;g<count;g++){
  auto n=read(),expected=read();std::vector<std::uint32_t>a(std::size_t(n)*n);
  in.read(reinterpret_cast<char*>(a.data()),a.size()*4);if(!in)return 4;
  std::uint32_t det=1;
  for(std::uint32_t k=0;k<n;k++){
   auto row=k;while(row<n && !a[std::size_t(row)*n+k])row++;if(row==n)return 5;
   if(row!=k){for(std::uint32_t j=k;j<n;j++)std::swap(a[std::size_t(k)*n+j],a[std::size_t(row)*n+j]);det=(P-det)%P;}
   auto pivot=a[std::size_t(k)*n+k];det=det*pivot%P;auto inverse=power(pivot,P-2);
   for(std::uint32_t i=k+1;i<n;i++){
    auto factor=a[std::size_t(i)*n+k]*inverse%P;if(!factor)continue;
    for(std::uint32_t j=k+1;j<n;j++){
     auto value=a[std::size_t(i)*n+j]+P-factor*a[std::size_t(k)*n+j]%P;
     a[std::size_t(i)*n+j]=value>=P?value-P:value;
    }
    a[std::size_t(i)*n+k]=0;
   }
  }
  if(det!=expected || !det){std::cerr<<"DETERMINANT MISMATCH "<<g<<' '<<det<<' '<<expected<<std::endl;return 6;}
  total+=n;
 }
 if(in.peek()!=std::ifstream::traits_type::eof())return 7;
 std::cout<<"{\"all_passed\":true,\"minors\":"<<count<<",\"total_rank\":"<<total<<"}\n";
}
