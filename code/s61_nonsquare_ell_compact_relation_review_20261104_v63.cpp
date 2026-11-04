// Independently verify the saved polynomial kernel against ORIGINAL source rows.
#define main independent_Lframe_inverse_audit_main
#include "s61_nonsquare_ell_Lframe_review_20261104_v63.cpp"
#undef main
using Mon8=array<int,8>;
Mon8 sum8(Mon8 a,const Mon8&b){for(int j=0;j<8;j++)a[j]+=b[j];return a;}
void term(Poly&z,const Poly&a,const Poly&b,long c=1){Poly v=multiply(a,b);add(z,scale(v,c));}
int main(int argc,char**argv){assert(argc==5);fmpq_mpoly_ctx_init(qc,3,ORD_LEX);{
 auto started=chrono::steady_clock::now();string work=argv[1],frozen=argv[2],pre=argv[4];ifstream mf(work+"/independent-source-metadata.txt");
 int n2,n3,n4;mf>>n2>>n3>>n4;assert(n2==29&&n3==91&&n4==239);vector<Mon8>m2(29),m3(91),m4(239);for(auto*arr:{&m2,&m3,&m4})for(auto&m:*arr)for(int&v:m)mf>>v;
 map<Mon8,int>i2,i3,i4;for(int j=0;j<29;j++)i2[m2[j]]=j;for(int j=0;j<91;j++)i3[m3[j]]=j;for(int j=0;j<239;j++)i4[m4[j]]=j;
 vector<array<int,3>>defs(319);int n;mf>>n;assert(n==319);for(auto&d:defs)for(int&v:d)mf>>v;string tail;assert(!(mf>>tail));
 Matrix Q=matrix(work+"/original-six-quadrics.txt",6,29);equalmat(Q,matrix(frozen+"/six-original-quadrics.txt",6,29));
 Matrix N=matrix(frozen+"/full319-kernel.txt",319,1),qc_saved=matrix(frozen+"/quadric6-cofactors.txt",6,29);
 Matrix G(2,vector<Poly>(91));ifstream source(argv[3]);int nr,nc,header;source>>nr>>nc>>header;assert(nr==57&&nc==91&&header==57);
 for(int r=0;r<57;r++)for(int j=0;j<91;j++){Poly p;load(source,p);if(r==42||r==47)G[r==42?0:1][j]=p;}assert(!(source>>tail));
 Poly a,e,d;fmpq_mpoly_gen(a.x,0,qc);fmpq_mpoly_gen(e.x,1,qc);fmpq_mpoly_gen(d.x,2,qc);
 Poly C=scale(a,100);add(C,scale(multiply(e,e),18));add(C,scale(e,-9));add(C,constant(-2));
 Poly J=scale(multiply(a,d),5);add(J,scale(e,-6));add(J,constant(-4));
 Poly de8=multiply(d,e);add(de8,constant(8));Poly em=e;add(em,constant(-1));
 Poly T=scale(a,10);add(T,scale(multiply(e,e),3));add(T,scale(e,-3));
 // Recover ell directly from -(de+8)q1+8q2+dq3 = ell*v3.
 // This uses no transformed-q helper or elliptic-coordinate certificate.
 vector<Poly>ell(8);for(int j=0;j<29;j++){
  Poly p=scale(Q[2][j],8);term(p,d,Q[3][j]);term(p,de8,Q[1][j],-1);if(zero(p))continue;
  auto ex=m2[j];assert(ex[3]>0);ex[3]--;int v=-1;for(int k=0;k<8;k++)if(ex[k]){assert(ex[k]==1&&v<0);v=k;}assert(v>=0);add(ell[v],p);
 }
 Matrix ell_saved=matrix(frozen+"/ell-original8.txt",8,1);for(int j=0;j<8;j++)assert(same(ell[j],ell_saved[j][0]));
 // The compact E,A,B data are checked, not trusted as a division algorithm.
 ifstream compact(frozen+"/compact-polynomials.txt");getline(compact,tail);assert(tail=="E8 A29 B29 gc8 delta1 anchor1");
 vector<Poly>E(8),A(29),B(29),gc(8);for(auto*row:{&E,&A,&B,&gc})for(auto&p:*row)load(compact,p);Poly delta,anchor;load(compact,delta);load(compact,anchor);assert(!(compact>>tail));
 Poly expectDelta=scale(multiply(C,C),-64000);assert(same(delta,expectDelta));
 Matrix G40=matrix(pre+"/genuine16x40.txt",16,40),Lm=matrix(pre+"/L-coefficients.txt",5,1);vector<Poly>diff(40);for(int j=0;j<40;j++)fmpq_mpoly_sub(diff[j].x,G40[6][j].x,G40[1][j].x,qc);
 for(int j=5;j<40;j++)assert(zero(diff[j]));
 vector<Poly>E_source(8);vector<int>old{0,1,2,4,7};for(int j=0;j<5;j++)E_source[old[j]]=multiply(diff[0],Lm[j][0]);for(int j=1;j<5;j++)add(E_source[old[j==1?0:j]],diff[j]);
 for(int v=0;v<8;v++)assert(same(E_source[v],E[v]));
 // Both G row clearings are den^3/C=-64000C². All91 original cubic coefficients.
 vector<Poly>relation(91);for(int j=0;j<91;j++){Poly p;fmpq_mpoly_sub(p.x,G[1][j].x,G[0][j].x,qc);relation[j]=multiply(delta,p);}
 for(int j=0;j<29;j++)for(int v=0;v<8;v++){if(zero(E[v]))continue;auto m=m2[j];m[v]++;term(relation[i3.at(m)],Q[1][j],E[v],-8100);}
 for(int j=0;j<29;j++){
  for(auto[v,c]:vector<pair<int,int>>{{0,3},{1,-1},{6,1}}){auto m=m2[j];m[v]++;add(relation[i3.at(m)],scale(A[j],-c));}
  auto m=m2[j];m[3]++;add(relation[i3.at(m)],scale(B[j],-1));
 }
 for(auto&p:relation)assert(zero(p));
 vector<Poly>ellZ(29),ellW(29),ellE(29);
 for(int v=0;v<8;v++){if(zero(ell[v]))continue;
  for(auto[w,c]:vector<pair<int,int>>{{0,3},{1,-1},{6,1}}){Mon8 m{};m[v]++;m[w]++;add(ellZ[i2.at(m)],scale(ell[v],c));}
  Mon8 m{};m[v]++;m[3]++;add(ellW[i2.at(m)],ell[v]);
  for(int w=0;w<8;w++){if(zero(E[w]))continue;Mon8 ex{};ex[v]++;ex[w]++;term(ellE[i2.at(ex)],ell[v],E[w]);}
 }
 for(int j=0;j<29;j++){
  Poly z=multiply(a,ellZ[j]);term(z,T,Q[1][j],-1);term(z,em,Q[3][j],3);term(z,a,Q[4][j],10);assert(zero(z));
  z=ellW[j];term(z,de8,Q[1][j]);add(z,scale(Q[2][j],-8));term(z,d,Q[3][j],-1);assert(zero(z));
 }
 Matrix expected_qc(6,vector<Poly>(29));for(int j=0;j<29;j++){
  term(expected_qc[1][j],a,ellE[j],-8100);term(expected_qc[1][j],T,A[j],-1);term(expected_qc[1][j],multiply(a,de8),B[j]);
  term(expected_qc[2][j],a,B[j],-8);term(expected_qc[3][j],em,A[j],3);term(expected_qc[3][j],multiply(a,d),B[j],-1);term(expected_qc[4][j],a,A[j],10);
 }equalmat(expected_qc,qc_saved);
 for(int v=0;v<8;v++)assert(same(gc[v],multiply(multiply(a,delta),ell[v])));
 // Audit the literal saved319-vector using its ORIGINAL row definition,
 // independently of the preceding six-cofactor recipe.
 vector<Poly>identity(239);int nonzero=0,quartic_count=0;long terms=0;
 for(int r=0;r<319;r++){
  const auto[ty,id,v]=defs[r];const Poly&co=N[r][0];terms+=fmpq_mpoly_length(co.x,qc);nonzero+=!zero(co);
  if(ty==0){assert(same(co,qc_saved[id][v]));for(int j=0;j<29;j++)term(identity[i4.at(sum8(m2[v],m2[j]))],co,Q[id][j]);}
  else if(ty==1){assert(v!=5);Poly expected;if(id==43)expected=scale(gc[v],-1);else if(id==49)expected=gc[v];assert(same(co,expected));if(zero(co))continue;
   assert(id==43||id==49);for(int j=0;j<91;j++){auto m=m3[j];m[v]++;term(identity[i4.at(m)],co,G[id==43?0:1][j]);}
  }else{assert(ty==2&&zero(co));quartic_count++;}
 }
 for(auto&p:identity)assert(zero(p));assert(quartic_count==33&&nonzero==123&&terms==19694);
 for(int row:{176,177,237})assert(zero(N[row][0]));Poly separator=N[182][0];add(separator,N[186][0]);
 Poly known=multiply(multiply(a,multiply(C,multiply(C,C))),J);fmpq_t qq;fmpq_init(qq);fmpq_set_si(qq,-1024000,243);fmpq_mpoly_scalar_mul_fmpq(known.x,known.x,qq,qc);fmpq_clear(qq);
 assert(same(separator,known)&&same(anchor,known));
 // Independently verify the24 selected rational functions from acceptedR4.
 // Only known a,F,B,C factors may occur in their denominators.
 Poly F=scale(a,500);add(F,scale(multiply(e,e),90));add(F,scale(e,3));add(F,constant(-24));
 Poly Bcore=scale(multiply(a,d),50);add(Bcore,scale(multiply(multiply(e,e),d),18));add(Bcore,scale(multiply(e,d),-9));add(Bcore,scale(d,-2));add(Bcore,scale(e,60));add(Bcore,constant(-40));
 ifstream ratios(work+"/residual4-selected-ratios.txt");int rr,cc;ratios>>rr>>cc;assert(rr==4&&cc==6);Matrix nums(4,vector<Poly>(6)),dens(4,vector<Poly>(6));
 for(int r=0;r<4;r++)for(int j=0;j<6;j++){
  load(ratios,nums[r][j]);load(ratios,dens[r][j]);assert(!zero(dens[r][j]));Poly rem=dens[r][j];
  for(auto&fac:vector<Poly>{a,F,Bcore,C})while(fmpq_mpoly_total_degree_si(rem.x,qc)>0){Poly q;if(!fmpq_mpoly_divides(q.x,rem.x,fac.x,qc))break;rem=q;}
  assert(fmpq_mpoly_is_fmpq(rem.x,qc)&&!zero(rem));if(j<4)assert(r==j?same(nums[r][j],dens[r][j]):zero(nums[r][j]));
 }
 Poly last=multiply(nums[3][4],dens[3][5]);add(last,multiply(nums[3][5],dens[3][4]));assert(zero(last));assert(!(ratios>>tail));
 cout<<"{\"all_passed\":true,\"all35_Lframe_difference_coefficients_zero\":true,\"E_pulled_back_directly\":true,\"original_ell_recovered_from_q6\":true,\"all91_compact_original_cubic_coefficients\":true,\"all58_six_quad_source_change_coefficients\":true,\"all174_quadric_cofactor_entries\":true,\"literal319_vector_dot_original_source_all239_zero\":true,\"all33_quartic_cofactors_zero\":true,\"kernel_nonzero_entries\":"<<nonzero<<",\"kernel_terms\":"<<terms<<",\"known_FULLU_unit_separator_checked\":true,\"residual4_selected24_ratios_and_known_unit_denominators_checked\":true,\"outside37_gate_checked\":true,\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"}"<<endl;
 }fmpq_mpoly_ctx_clear(qc);flint_cleanup();return 0;}
