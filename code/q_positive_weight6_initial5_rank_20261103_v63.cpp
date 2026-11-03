#include <flint/nmod_mat.h>
#include <fstream>
#include <iostream>
int main(int argc,char**argv){
 if(argc!=2)return 2;
 std::ifstream in(argv[1]);long n,m;unsigned long p;in>>n>>m>>p;
 nmod_mat_t A;nmod_mat_init(A,n,m,p);
 for(long i=0;i<n;i++)for(long j=0;j<m;j++){unsigned long v;in>>v;nmod_mat_entry(A,i,j)=v;}
 if(!in)return 3;
 long rank=nmod_mat_rref(A);std::cout<<rank<<'\n';
 for(long i=0;i<rank;i++){long j=0;while(j<m&&!nmod_mat_entry(A,i,j))j++;std::cout<<j<<' ';}
 std::cout<<'\n';nmod_mat_clear(A);return 0;
}
