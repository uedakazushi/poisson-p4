#include <flint/nmod_mat.h>
#include <flint/nmod.h>
#include <array>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include <cassert>
using namespace std;
constexpr int eqs=134,n=25,n2=325,n3=2925,m3=3350,sy=425,p=32003,low=351,target=217;
struct Z {nmod_mat_t m; Z(int r,int c){nmod_mat_init(m,r,c,p);} ~Z(){nmod_mat_clear(m);} ulong& at(int i,int j){return nmod_mat_entry(m,i,j);}};
static vector<int> selectrows(Z &a,int rank){
    nmod_mat_t E;nmod_mat_init_set(E,a.m);vector<slong> perm(a.m->r);iota(perm.begin(),perm.end(),0);
    auto r=nmod_mat_lu(perm.data(),E,0);assert(r==rank);nmod_mat_clear(E);return {perm.begin(),perm.begin()+rank};
}
int main(int argc,char**argv){
 assert(argc==4);ifstream f(argv[1]);ofstream o(argv[2]);int np;f>>np;
 int pair[n][n],tri[n][n][n],ct=0;
 for(int i=0;i<n;i++)for(int j=i;j<n;j++)pair[i][j]=pair[j][i]=ct++;assert(ct==n2);ct=0;
 for(int i=0;i<n;i++)for(int j=i;j<n;j++)for(int k=j;k<n;k++){array<int,3> a{i,j,k};do{tri[a[0]][a[1]][a[2]]=ct;}while(next_permutation(a.begin(),a.end()));ct++;}assert(ct==n3);
 Z H(eqs,n2),H3(m3,n3);
 for(int r=0;r<eqs;r++)for(int i=0;i<n;i++)for(int j=i;j<n;j++){
  ulong v;f>>v;assert(f&&v<p);H.at(r,pair[i][j])=v;for(int k=0;k<n;k++)H3.at(r*n+k,tri[i][j][k])=v;
 }
 Z HT(n2,eqs);nmod_mat_transpose(HT.m,H.m);auto p2=selectrows(HT,eqs),p3=selectrows(H3,n3);
 vector<int>cf,rr;vector<bool>used2(n2),used3(m3);for(int j:p2)used2[j]=true;for(int j:p3)used3[j]=true;
 for(int j=0;j<n2;j++)if(!used2[j])cf.push_back(j);for(int j=0;j<m3;j++)if(!used3[j])rr.push_back(j);
 // C3 * H3[p3] = H3[rr]. Thus row-reduction multipliers are independent of all y's.
 Z c3(sy,n3),A3(n3,n3),B3(n3,sy),X3(n3,sy);
 for(int i=0;i<n3;i++){for(int j=0;j<n3;j++)A3.at(i,j)=H3.at(p3[j],i);for(int j=0;j<sy;j++)B3.at(i,j)=H3.at(rr[j],i);}
 ulong determinant3=nmod_mat_det(A3.m);assert(determinant3);
 assert(nmod_mat_solve(X3.m,A3.m,B3.m));nmod_mat_transpose(c3.m,X3.m);
 Z A2(eqs,eqs),A2inv(eqs,eqs);for(int i=0;i<eqs;i++)for(int j=0;j<eqs;j++)A2.at(i,j)=H.at(i,p2[j]);
 assert(nmod_mat_inv(A2inv.m,A2.m));o<<"fixed "<<eqs<<" "<<n3<<" "<<nmod_mat_det(A2.m)<<" "<<determinant3<<"\n";
 for(int col:p2)o<<col<<" ";o<<"\n";for(int row:p3)o<<row<<" ";o<<"\n";o.flush();
 for(int sample=0;sample<np;sample++){
  int label;f>>label;Z R(eqs,26);for(int i=0;i<eqs;i++)for(int j=0;j<26;j++){f>>R.at(i,j);assert(f&&R.at(i,j)<p);}
  Z L(m3,low);for(int i=0;i<eqs;i++)for(int j=0;j<n;j++){
    for(int k=0;k<n;k++)L.at(i*n+j,pair[j][k])=R.at(i,k);
    L.at(i*n+j,n2+j)=R.at(i,n);
  }
  Z highrows(n3,low),ex(sy,low),red(sy,low);
  for(int i=0;i<n3;i++)for(int j=0;j<low;j++)highrows.at(i,j)=L.at(p3[i],j);
  for(int i=0;i<sy;i++)for(int j=0;j<low;j++)ex.at(i,j)=L.at(rr[i],j);
  nmod_mat_mul(red.m,c3.m,highrows.m);nmod_mat_sub(red.m,ex.m,red.m);
  Z ph(sy,eqs),ret(sy,target),redf(eqs,target);
  for(int i=0;i<sy;i++){
    for(int j=0;j<eqs;j++)ph.at(i,j)=red.at(i,p2[j]);
    for(int j=0;j<int(cf.size());j++)ret.at(i,j)=red.at(i,cf[j]);
    for(int j=0;j<26;j++)ret.at(i,cf.size()+j)=red.at(i,n2+j);
  }
  for(int i=0;i<eqs;i++){
   for(int j=0;j<int(cf.size());j++)redf.at(i,j)=H.at(i,cf[j]);
   for(int j=0;j<26;j++)redf.at(i,cf.size()+j)=R.at(i,j);
  }
  Z xp(sy,eqs),r2(sy,target);nmod_mat_mul(xp.m,ph.m,A2inv.m);nmod_mat_mul(r2.m,xp.m,redf.m);nmod_mat_sub(r2.m,ret.m,r2.m);
  assert(np==1);
  std::ifstream probes(argv[3]);
  for (int which=0;which<5;which++){
    std::vector<ulong> point(n);for(auto & value:point){probes>>value;assert(probes&&value<p);}
    Z derivative(eqs,20),coeff(sy,eqs),der_out(sy,20);
    for(int a=0;a<eqs;a++)for(int b=0;b<20;b++){probes>>derivative.at(a,b);assert(probes&&derivative.at(a,b)<p);}
    for(int u=0;u<sy;u++){
      for(int k=0;k<n3;k++){int row=p3[k]/n,mul=p3[k]%n;
        coeff.at(u,row)=(coeff.at(u,row)+p-(c3.at(u,k)*point[mul])%p)%p;
      }
      coeff.at(u,rr[u]/n)=(coeff.at(u,rr[u]/n)+point[rr[u]%n])%p;
      for(int row=0;row<eqs;row++)coeff.at(u,row)=(coeff.at(u,row)+p-xp.at(u,row))%p;
    }
    nmod_mat_mul(der_out.m,coeff.m,derivative.m);
    for(int i=0;i<sy;i++){for(int b=0;b<20;b++)o<<der_out.at(i,b)<<" ";o<<"\n";}
  }
  int rank=nmod_mat_rank(r2.m);o<<label<<" "<<rank<<"\n";
  for(int i=0;i<sy;i++){for(int j=0;j<target;j++)o<<r2.at(i,j)<<" ";o<<"\n";}o.flush();
 }
}
