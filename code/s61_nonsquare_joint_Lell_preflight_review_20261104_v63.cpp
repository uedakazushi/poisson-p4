// Inverse coordinate replay of actual quadrics, not forward producer transport.
#define main independent_Lframe_inverse_audit_main
#include "s61_nonsquare_ell_Lframe_review_20261104_v63.cpp"
#undef main
using E8=array<int,8>;using L8=map<E8,Poly>;
Poly qq(long n,long d=1){Poly p;fmpq_t q;fmpq_init(q);fmpq_set_si(q,n,d);fmpq_mpoly_set_fmpq(p.x,q,qc);fmpq_clear(q);return p;}
L8 prod8(const L8&a,const L8&b){L8 z;for(auto&[e,c]:a)for(auto&[f,d]:b){E8 m;for(int j=0;j<8;j++)m[j]=e[j]+f[j];add(z[m],multiply(c,d));}return z;}
int main(int argc,char**argv){assert(argc==4);fmpq_mpoly_ctx_init(qc,3,ORD_LEX);{
 auto started=chrono::steady_clock::now();string own=argv[1],frozen=argv[2],pre=argv[3];ifstream meta(own+"/independent-monomials.txt");
 int a2,a3,a4;meta>>a2>>a3>>a4;assert(a2==29&&a3==91&&a4==239);vector<E8>oldm2(29),m2(29),m3(91);for(auto*arr:{&oldm2,&m2,&m3})for(auto&m:*arr)for(int&v:m)meta>>v;
 map<E8,int>i2,i3;for(int j=0;j<29;j++)i2[m2[j]]=j;for(int j=0;j<91;j++)i3[m3[j]]=j;
 Matrix Q=matrix(own+"/original-six-quadrics.txt",6,29),F=matrix(frozen+"/coordinate-forward7x7.txt",7,7),N=matrix(frozen+"/coordinate-inverse-numerators7x7.txt",7,7);
 Matrix nums=matrix(frozen+"/six-prime-quadrics-numerators.txt",6,29),dens=matrix(frozen+"/six-prime-quadrics-denominators.txt",6,1),Lm=matrix(pre+"/L-coefficients.txt",5,1);
 Poly a,e,d;fmpq_mpoly_gen(a.x,0,qc);fmpq_mpoly_gen(e.x,1,qc);fmpq_mpoly_gen(d.x,2,qc);Poly C=scale(a,100);add(C,scale(multiply(e,e),18));add(C,scale(e,-9));add(C,constant(-2));
 Poly de8=multiply(d,e);add(de8,constant(8));Poly em=e;add(em,constant(-1));Poly T=scale(a,10);add(T,scale(multiply(e,e),3));add(T,scale(e,-3));
 // Independent physical-v8 -> oldX7,H substitution. v6=Z+x1-3x0.
 vector<L8>oldimages(8);vector<int>idx{0,1,2,6,3,7,-1,4};
 for(int v=0;v<8;v++)if(idx[v]>=0){E8 ex{};ex[idx[v]]=1;oldimages[v][ex]=constant(1);}
 for(auto[v,c]:vector<pair<int,int>>{{5,1},{1,1},{0,-3}}){E8 ex{};ex[v]=1;oldimages[6][ex]=constant(c);}
 vector<L8>oldmons;for(auto&m:oldm2){L8 p;p[E8{}]=constant(1);for(int v=0;v<8;v++)for(int k=0;k<m[v];k++)p=prod8(p,oldimages[v]);oldmons.push_back(move(p));}
 Matrix Qold(6,vector<Poly>(29));for(int r=0;r<6;r++)for(int j=0;j<29;j++)for(auto&[ex,p]:oldmons[j])add(Qold[r][i2.at(ex)],multiply(Q[r][j],p));
 // Recover full ell from the actual q6/W identity, including BOTH Z,W terms.
 vector<Poly>ell(7);for(int j=0;j<29;j++){
  Poly co=scale(Qold[2][j],8);add(co,multiply(d,Qold[3][j]));add(co,scale(multiply(de8,Qold[1][j]),-1));if(zero(co))continue;
  auto ex=m2[j];assert(ex[6]>0);ex[6]--;int v=-1;for(int k=0;k<8;k++)if(ex[k]){assert(ex[k]==1&&v<0);v=k;}assert(v>=0&&v<7);add(ell[v],co);
 }
 auto oldell=matrix(frozen+"/ell-x7.txt",7,1);for(int j=0;j<7;j++)assert(same(ell[j],oldell[j][0]));
 // Source-bound forward map and L restricted-source derivative identity.
 Matrix expectedF(7,vector<Poly>(7));for(int j=0;j<5;j++)expectedF[0][j]=Lm[j][0];expectedF[1]=ell;
 expectedF[2][0]=expectedF[3][3]=expectedF[5][5]=expectedF[6][6]=constant(1);expectedF[4][2]=qq(7,2);expectedF[4][4]=constant(3);equalmat(F,expectedF);
 for(int j=0;j<5;j++){E8 ex{};ex[1]++;ex[j]++;Poly der=scale(Qold[1][i2.at(ex)],j==1?2:1);assert(same(scale(multiply(C,Lm[j][0]),4),scale(der,-405)));}
 // Actual determinant reduces to this fixed3x3 expression because four rows
 // are coordinate rows; no parametric determinant solver is used.
 Poly g=scale(ell[2],3);add(g,multiply(qq(-7,2),ell[4]));Poly h=scale(F[0][2],3);add(h,multiply(qq(-7,2),F[0][4]));
 Poly jac=multiply(F[0][1],g);add(jac,scale(multiply(ell[1],h),-1));Poly expectedJac=multiply(multiply(a,C),multiply(d,d));expectedJac=multiply(qq(3200,81),expectedJac);
 assert(same(jac,expectedJac));assert(same(jac,matrix(frozen+"/coordinate-Jacobian.txt",1,1)[0][0]));
 for(int i=0;i<7;i++)for(int j=0;j<7;j++){
  Poly fn,nf;for(int k=0;k<7;k++){add(fn,multiply(F[i][k],N[k][j]));add(nf,multiply(N[i][k],F[k][j]));}
  assert(i==j?same(fn,jac):zero(fn));assert(i==j?same(nf,jac):zero(nf));
 }
 // Substitute newY=F*oldX BACK into all six saved primequadrics.
 vector<L8>back(8);for(int v=0;v<7;v++)for(int j=0;j<7;j++){E8 ex{};ex[j]=1;back[v][ex]=F[v][j];}E8 eh{};eh[7]=1;back[7][eh]=constant(1);
 Matrix Pback(6,vector<Poly>(29));for(int j=0;j<29;j++){
  L8 p;p[E8{}]=constant(1);for(int v=0;v<8;v++)for(int k=0;k<m2[j][v];k++)p=prod8(p,back[v]);
  for(auto&[ex,co]:p)for(int r=0;r<6;r++)add(Pback[r][i2.at(ex)],multiply(nums[r][j],co));
 }
 vector<Poly>ellZ(29),ellW(29);for(int v=0;v<7;v++){E8 ex{};ex[v]++;ex[5]++;add(ellZ[i2.at(ex)],ell[v]);ex={};ex[v]++;ex[6]++;add(ellW[i2.at(ex)],ell[v]);}
 for(int j=0;j<29;j++){
  assert(same(Pback[0][j],scale(multiply(dens[0][0],Qold[1][j]),8100)));
  Poly p=multiply(e,Qold[1][j]);add(p,scale(Qold[3][j],-1));assert(same(scale(multiply(a,Pback[1][j]),4),scale(multiply(dens[1][0],p),81)));
  assert(same(Pback[2][j],multiply(dens[2][0],ellZ[j]))&&same(Pback[3][j],multiply(dens[3][0],ellW[j])));
  assert(same(Pback[4][j],multiply(dens[4][0],Qold[0][j]))&&same(Pback[5][j],multiply(dens[5][0],Qold[5][j])));
  Poly z=multiply(a,ellZ[j]);add(z,scale(multiply(T,Qold[1][j]),-1));add(z,scale(multiply(em,Qold[3][j]),3));add(z,scale(multiply(a,Qold[4][j]),10));assert(zero(z));
 }
 for(auto&row:dens){assert(!zero(row[0]));Poly rem=row[0];for(auto&factor:vector<Poly>{a,d,C})while(fmpq_mpoly_total_degree_si(rem.x,qc)>0){Poly q;if(!fmpq_mpoly_divides(q.x,rem.x,factor.x,qc))break;rem=q;}assert(!zero(rem)&&fmpq_mpoly_is_fmpq(rem.x,qc));}
 vector<E8>lead(4);lead[0][0]=2;lead[1][0]=lead[1][2]=1;lead[2][1]=lead[2][5]=1;lead[3][1]=lead[3][6]=1;long correctionterms=0;bool outside=false;
 for(int r=0;r<4;r++)for(int j=0;j<29;j++){
  if(m2[j]==lead[r])assert(same(nums[r][j],dens[r][0]));else{
   if(m2[j]>lead[r])assert(zero(nums[r][j]));if(r>=2)assert(zero(nums[r][j]));
   if(r<2&&!zero(nums[r][j])){assert(m2[j][5]||m2[j][6]);correctionterms+=fmpq_mpoly_length(nums[r][j].x,qc);if(m2[j]!=lead[2]&&m2[j]!=lead[3])outside=true;}
  }
 }
 assert(correctionterms==82&&outside);
 Matrix S=matrix(frozen+"/nonzero-cubic-S-numerator.txt",1,91),Sd=matrix(frozen+"/nonzero-cubic-S-denominator.txt",1,1);vector<Poly>cross(91),exactzero(91);
 assert(!zero(Sd[0][0]));Poly sr=Sd[0][0];for(auto&factor:vector<Poly>{a,d,C})while(fmpq_mpoly_total_degree_si(sr.x,qc)>0){Poly q;if(!fmpq_mpoly_divides(q.x,sr.x,factor.x,qc))break;sr=q;}assert(!zero(sr)&&fmpq_mpoly_is_fmpq(sr.x,qc));
 for(int j=0;j<29;j++){
  auto ex=m2[j];ex[2]++;add(cross[i3.at(ex)],multiply(nums[0][j],dens[1][0]));ex=m2[j];ex[0]++;add(cross[i3.at(ex)],scale(multiply(nums[1][j],dens[0][0]),-1));
  ex=m2[j];ex[6]++;add(exactzero[i3.at(ex)],multiply(nums[2][j],dens[3][0]));ex=m2[j];ex[5]++;add(exactzero[i3.at(ex)],scale(multiply(nums[3][j],dens[2][0]),-1));
 }
 Poly common=multiply(dens[0][0],dens[1][0]);bool any=false;for(int j=0;j<91;j++){
  assert(same(multiply(S[0][j],common),multiply(cross[j],Sd[0][0])));assert(zero(exactzero[j]));if(!zero(S[0][j])){any=true;assert(m3[j][5]||m3[j][6]);}
 }assert(any&&count_terms(S)==72);
 assert(count_terms(nums)==1633);long nz=0;for(auto&r:nums)for(auto&p:r)nz+=!zero(p);assert(nz==72);
 cout<<"{\"all_passed\":true,\"all174_primequad_inverse_coordinate_identities\":true,\"ell_and_L_and_forward_frame_rebuilt_from_actualq6\":true,\"actual_Jacobian_and_both49_inverse_identities\":true,\"all_primequad_denominators_only_a_d_C_units\":true,\"all4_monic_leading_terms\":true,\"nonzero_ZW_correction_terms\":"<<correctionterms<<",\"actual_four_monomial_ideal_equality_rejected\":true,\"nonzero_S_all91_coefficients_and72terms_checked\":true,\"W_p3_minus_Z_p4_all91_coefficients_zero\":true,\"primequad_terms1633_nonzeros72\":true,\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"}"<<endl;
 }fmpq_mpoly_ctx_clear(qc);flint_cleanup();return 0;}
