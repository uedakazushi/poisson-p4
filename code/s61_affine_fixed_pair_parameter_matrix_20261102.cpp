#include <flint/nmod_mat.h>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
struct Mat{nmod_mat_t a;Mat(long r,long c){nmod_mat_init(a,r,c,32003);}~Mat(){nmod_mat_clear(a);}};
int terms(int d){return (d+1)*(d+2)/2;}
int ix(int i,int j,int degree){return i*(degree+1)-i*(i-1)/2+j;}
int main(int argc,char**argv) {
 assert(argc==4);ifstream f(argv[1]);ofstream o(argv[2]);int D=stoi(argv[3]),mrs=terms(D-3),rows=782*mrs,cols=142*terms(D-2)+23*terms(D-1)+terms(D);
 Mat M(rows,cols);int limits[166];for(int i=0;i<166;i++)limits[i]=i<142?1:i<165?2:3;
 vector<int>offset(166);for(int j=1;j<166;j++)offset[j]=offset[j-1]+terms(D-3+limits[j-1]);
 for(int z=0;z<10;z++){int i,j;f>>i>>j;for(int r=0;r<782;r++)for(int c=0;c<166;c++){
  ulong val;f>>val;assert(f&&val<32003);if(!val)continue;int b=limits[c];assert(i+j<=b);
  for(int a=0;a<=D-3;a++)for(int u=0;u<=D-3-a;u++)
   nmod_mat_entry(M.a,ix(a,u,D-3)*782+r,offset[c]+ix(i+a,j+u,D-3+b))=val;
 }}
 std::vector<slong> perm(rows);for(long i=0;i<rows;i++)perm[i]=i;
 nmod_mat_t E;nmod_mat_init_set(E,M.a);long rank=nmod_mat_lu(perm.data(),E,0);nmod_mat_clear(E);
 assert(rank==cols);Mat sq(cols,cols);for(long i=0;i<cols;i++)for(long j=0;j<cols;j++)nmod_mat_entry(sq.a,i,j)=nmod_mat_entry(M.a,perm[i],j);
 unsigned long det=nmod_mat_det(sq.a);assert(det!=0);o<<det<<"\n";for(long i=0;i<cols;i++)o<<perm[i]<<" ";o<<"\n";
 o<<D<<" "<<rows<<" "<<cols<<" "<<rank<<" "<<cols-rank<<"\n";o.flush();flint_cleanup();
}
