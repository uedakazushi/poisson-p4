#include <flint/nmod_mat.h>
#include <cassert>
#include <fstream>
#include <iostream>
int main(int argc,char**argv){
    assert(argc==2);std::ifstream f(argv[1]);long n,nnz;unsigned long p;f>>n>>p>>nnz;
    nmod_mat_t A;nmod_mat_init(A,n,n,p);
    for(long k=0;k<nnz;k++){long i,j;unsigned long c;f>>i>>j>>c;assert(i>=0&&i<n&&j>=0&&j<n);nmod_mat_entry(A,i,j)=(nmod_mat_entry(A,i,j)+c)%p;}
    assert(f.good());long rank=nmod_mat_rank(A);std::cout<<"SELECTED_SQUARE_RANK "<<rank<<" SIZE "<<n<<" PRIME "<<p<<"\n";
    nmod_mat_clear(A);flint_cleanup();return rank==n?0:2;
}
