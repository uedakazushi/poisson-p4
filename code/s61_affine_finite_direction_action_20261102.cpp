#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod.h>
#include <array>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include <cassert>
using namespace std;
constexpr int eqs=134,n=24,n2=300,n3=2600,m3=3216,p=32003;
struct Z {nmod_mat_t m; Z(int r,int c){nmod_mat_init(m,r,c,p);} ~Z(){nmod_mat_clear(m);} ulong& at(int i,int j){return nmod_mat_entry(m,i,j);}};
static vector<int> selectrows(Z &a,int rank){
    nmod_mat_t E;nmod_mat_init_set(E,a.m);vector<slong> perm(a.m->r);iota(perm.begin(),perm.end(),0);
    auto r=nmod_mat_lu(perm.data(),E,0);assert(r==rank);nmod_mat_clear(E);return {perm.begin(),perm.begin()+rank};
}
#include <memory>
using Matrix = std::unique_ptr<Z>;

static Matrix identity(int size) {
 auto result=std::make_unique<Z>(size,size);
 nmod_mat_one(result->m);
 return result;
}

// Keep the combinations of the current columns on which the test vanishes.
static Matrix restrict_kernel(const Z& columns, const Z& test) {
 int width=columns.m->c;
 Z equations(test.m->r,width),kernel(width,width);
 nmod_mat_mul(equations.m,test.m,columns.m);
 int nullity=nmod_mat_nullspace(kernel.m,equations.m);
 Z check(equations.m->r,width);nmod_mat_mul(check.m,equations.m,kernel.m);assert(nmod_mat_is_zero(check.m));
 Z change(width,nullity);
 for(int i=0;i<width;i++)for(int j=0;j<nullity;j++)change.at(i,j)=kernel.at(i,j);
 auto result=std::make_unique<Z>(columns.m->r,nullity);
 nmod_mat_mul(result->m,columns.m,change.m);
 assert(nmod_mat_rank(result->m)==nullity);
 return result;
}

static Matrix commute_and_close(const std::vector<Matrix>& operators,std::vector<int>& trace) {
 int size=operators.front()->m->r;
 Matrix basis=identity(size);
 for(int i=0;i<int(operators.size()) && basis->m->c;i++) {
  for(int j=0;j<i && basis->m->c;j++) {
   int width=basis->m->c;
   Z a(size,width),b(size,width),ab(size,width),ba(size,width),test(size,width);
   nmod_mat_mul(a.m,operators[i]->m,basis->m);
   nmod_mat_mul(b.m,operators[j]->m,basis->m);
   nmod_mat_mul(ab.m,operators[i]->m,b.m);
   nmod_mat_mul(ba.m,operators[j]->m,a.m);
   nmod_mat_sub(test.m,ab.m,ba.m);
   if(nmod_mat_is_zero(test.m))continue;
   auto next=restrict_kernel(*identity(width),test);
   auto reduced=std::make_unique<Z>(size,next->m->c);
   nmod_mat_mul(reduced->m,basis->m,next->m);
   basis=std::move(reduced);
  }
 }
 trace.push_back(basis->m->c);
 while(basis->m->c!=0 && basis->m->c!=size) {
  int width=basis->m->c;
  auto rows=selectrows(*basis,width);
  Z pivot(width,width),inverse(width,width),projector(size,width);
  for(int i=0;i<width;i++)for(int j=0;j<width;j++)pivot.at(i,j)=basis->at(rows[i],j);
  assert(nmod_mat_inv(inverse.m,pivot.m));nmod_mat_mul(projector.m,basis->m,inverse.m);
  auto change=identity(width);
  for(const auto& op:operators){
   Z image(size,width),selected(width,width),res(size,width);
   nmod_mat_mul(image.m,op->m,basis->m);
   for(int i=0;i<width;i++)for(int j=0;j<width;j++)selected.at(i,j)=image.at(rows[i],j);
   nmod_mat_mul(res.m,projector.m,selected.m);nmod_mat_sub(res.m,image.m,res.m);
   if(!nmod_mat_is_zero(res.m))change=restrict_kernel(*change,res);
   if(change->m->c==0)break;
  }
  if(change->m->c==width)break;
  auto reduced=std::make_unique<Z>(size,change->m->c);
  nmod_mat_mul(reduced->m,basis->m,change->m);
  basis=std::move(reduced);trace.push_back(basis->m->c);
 }
 return basis;
}

int main(int argc,char**argv){
 assert(argc==3);ifstream f(argv[1]);ofstream o(argv[2]);int np;f>>np;
 int pair[n][n],tri[n][n][n],ct=0;
 for(int i=0;i<n;i++)for(int j=i;j<n;j++)pair[i][j]=pair[j][i]=ct++;assert(ct==n2);ct=0;
 for(int i=0;i<n;i++)for(int j=i;j<n;j++)for(int k=j;k<n;k++){array<int,3> a{i,j,k};do{tri[a[0]][a[1]][a[2]]=ct;}while(next_permutation(a.begin(),a.end()));ct++;}assert(ct==n3);
 Z H(eqs,n2),H3(m3,n3);
 for(int r=0;r<eqs;r++)for(int i=0;i<n;i++)for(int j=i;j<n;j++){
  ulong v;f>>v;assert(f&&v<p);H.at(r,pair[i][j])=v;for(int k=0;k<n;k++)H3.at(r*n+k,tri[i][j][k])=v;
 }
 Z HT(n2,eqs);nmod_mat_transpose(HT.m,H.m);auto p2=selectrows(HT,eqs),p3=selectrows(H3,n3);
 vector<int>cf;vector<bool>used2(n2);for(int j:p2)used2[j]=true;
 for(int j=0;j<n2;j++)if(!used2[j])cf.push_back(j);
 // The 191 symbols are the constant, 24 linears and 166 remaining quadrics.
 const int m=191;
 Z A2(eqs,eqs),inv2(eqs,eqs),mat3(n3,n3),rhs3(n3,m),red3(n3,m);
 for(int i=0;i<eqs;i++)for(int j=0;j<eqs;j++)A2.at(i,j)=H.at(i,p2[j]);
 assert(nmod_mat_inv(inv2.m,A2.m));
 for(int i=0;i<n3;i++)for(int j=0;j<n3;j++)mat3.at(i,j)=H3.at(p3[i],j);

 o<<"pivots "<<nmod_mat_det(A2.m)<<" "<<nmod_mat_det(mat3.m)<<"\n";
 for(auto i:p2)o<<i<<" ";o<<"\n";for(auto i:p3)o<<i<<" ";o<<"\n";
 for(int sample=0;sample<np;sample++){
  int label;f>>label;Z L(eqs,n+1);for(int i=0;i<eqs;i++)for(int j=0;j<=n;j++){f>>L.at(i,j);assert(f&&L.at(i,j)<p);}
  Z rf(eqs,m),reduced(eqs,m),red2(n2,m);
  for(int i=0;i<eqs;i++){rf.at(i,0)=L.at(i,n);for(int k=0;k<n;k++)rf.at(i,1+k)=L.at(i,k);for(int k=0;k<int(cf.size());k++)rf.at(i,1+n+k)=H.at(i,cf[k]);}
  nmod_mat_mul(reduced.m,inv2.m,rf.m);
  for(int i=0;i<eqs;i++)for(int j=0;j<m;j++)red2.at(p2[i],j)=(p-reduced.at(i,j))%p;
  for(int k=0;k<int(cf.size());k++)red2.at(cf[k],1+n+k)=1;
  // Reduce the low parts of the selected X_j F_r rows to the retained symbols.
  for(int i=0;i<n3;i++){int r=p3[i]/n,j=p3[i]%n;
   for(int k=0;k<m;k++){ulong v= (k==j+1?L.at(r,n):0);
    for(int t=0;t<n;t++)v+=L.at(r,t)*red2.at(pair[j][t],k);
    rhs3.at(i,k)=(p-v%p)%p;
   }
  }
  assert(nmod_mat_solve(red3.m,mat3.m,rhs3.m));

  std::vector<Matrix> operators;
  for(int i=0;i<n;i++) {
   auto op=std::make_unique<Z>(m,m);op->at(0,i+1)=1;
   for(int j=0;j<n;j++)for(int k=0;k<m;k++)op->at(j+1,k)=red2.at(pair[i][j],k);
   int u=0;
   for(int j=0;j<n;j++)for(int t=j;t<n;t++)if(!used2[pair[j][t]]){
    assert(cf[u]==pair[j][t]);
    for(int k=0;k<m;k++)op->at(n+1+u,k)=red3.at(tri[i][j][t],k);
    u++;
   }
   operators.push_back(std::move(op));
  }
  std::vector<int> dimensions;auto closed=commute_and_close(operators,dimensions);
  o<<"exact "<<label<<" trace";for(auto k:dimensions)o<<" "<<k;
  o<<" dimension "<<closed->m->c<<"\n";

  if(closed->m->c){
   // Determine the homogenizer on the invariant kernel and check every row.
   int k=closed->m->c;auto rows=selectrows(*closed,k);
   Z pre(k,k),rhs(k,k),action(k,k),full(m,k),re(m,k),power(k,k);
   nmod_mat_mul(full.m,operators[0]->m,closed->m);
   for(int i=0;i<k;i++)for(int j=0;j<k;j++){pre.at(i,j)=closed->at(rows[i],j);rhs.at(i,j)=full.at(rows[i],j);}
   assert(nmod_mat_solve(action.m,pre.m,rhs.m));nmod_mat_mul(re.m,closed->m,action.m);assert(nmod_mat_equal(re.m,full.m));
   nmod_mat_pow(power.m,action.m,m);
   o<<"t "<<nmod_mat_rank(action.m)<<" stable "<<nmod_mat_rank(power.m)<<" dim "<<k<<"\n";
  }else o<<"t 0 stable 0 dim 0\n";

  o.flush();

 }
 std::vector<Matrix> toy;toy.push_back(std::make_unique<Z>(3,3));toy.push_back(std::make_unique<Z>(3,3));
 for(int i=0;i<3;i++){toy[0]->at(i,i)=i;for(int j=0;j<3;j++)if(i!=j)toy[1]->at(i,j)=1;}
 std::vector<int> test_dimensions;auto zero=commute_and_close(toy,test_dimensions);
 assert(test_dimensions==std::vector<int>({1,0}));o<<"toy trace 1 0\n";

}
