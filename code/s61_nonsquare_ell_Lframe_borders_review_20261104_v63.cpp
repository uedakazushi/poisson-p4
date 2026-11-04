// Literal replay of saved polynomial row identities. No pivot choice/gcd search.
#define main independent_Lframe_inverse_audit_main
#include "s61_nonsquare_ell_Lframe_review_20261104_v63.cpp"
#undef main
void compare(const Matrix&a,const Matrix&b){equalmat(a,b);}
int main(int argc,char**argv){assert(argc==4);fmpq_mpoly_ctx_init(qc,3,ORD_LEX);{
 auto start=chrono::steady_clock::now();string input=argv[1],saved=argv[2];ifstream meta(argv[3]);
 vector<pair<int,int>>piv(16);for(auto&[r,c]:piv)meta>>r>>c;assert(meta);
 vector<int>available21,available16;for(int j=0;j<31;j++){
  bool q=true,g=true;for(int k=0;k<10;k++)if(piv[k].second==j)q=false;for(int k=0;k<15;k++)if(piv[k].second==j)g=false;
  if(q)available21.push_back(j);if(g)available16.push_back(j);}
 assert(available21.size()==21&&available16.size()==16);
 vector<int>remaining11{10,13,14,15,16,17,19,22,23,24,25};
 Matrix Q=matrix(input+"/quad10x26.txt",10,26),G=matrix(input+"/genuine16x31.txt",16,31),A(26,vector<Poly>(31));
 for(int r=0;r<2;r++)for(int j=0;j<26;j++)A[r][j]=Q[8+r][j];for(int r=0;r<8;r++)for(int j=0;j<26;j++)A[r+2][j]=Q[r][j];for(int r=0;r<16;r++)A[10+r]=G[r];
 Matrix frozenR=matrix(saved+"/completed-R11x16.txt",11,16),frozenAfterQ=matrix(saved+"/completed-afterquad16x21.txt",16,21);
 bool afterquad=false,afterheavy=false;long normalizations=0,eliminations=0,nonconstant_normalizations=0;ifstream transcript(saved+"/operations.txt");string op;
 while(transcript>>op){
  if(op=="NORMALIZE"){
   int row;transcript>>row;assert(row>=0&&row<26);Poly g;load(transcript,g);assert(!zero(g));
   string s;transcript>>s;fmpq_t c;fmpq_init(c);assert(!fmpq_set_str(c,s.c_str(),10));fmpq_canonicalise(c);assert(!fmpq_is_zero(c));
   if(fmpq_mpoly_total_degree_si(g.x,qc)>0)nonconstant_normalizations++;
   Poly factor;fmpq_mpoly_scalar_mul_fmpq(factor.x,g.x,c,qc);assert(!zero(factor));
   for(int j=0;j<31;j++){Poly q,back;assert(fmpq_mpoly_divides(q.x,A[row][j].x,factor.x,qc));fmpq_mpoly_mul(back.x,q.x,factor.x,qc);assert(same(back,A[row][j]));A[row][j]=q;}
   normalizations++;fmpq_clear(c);
  }else{
   assert(op=="ELIMINATE");int t,p,c;transcript>>t>>p>>c;assert(t>=0&&t<26&&p>=0&&p<26&&p!=t&&p!=16&&c>=0&&c<31);
   if(!afterquad&&p==piv[10].first&&c==piv[10].second){
    Matrix snap(16,vector<Poly>(21));for(int r=0;r<16;r++)for(int j=0;j<21;j++)snap[r][j]=A[r+10][available21[j]];
    compare(snap,frozenAfterQ);afterquad=true;
   }
   if(!afterheavy&&p==piv[15].first&&c==piv[15].second){
    Matrix snap(11,vector<Poly>(16));for(int r=0;r<11;r++)for(int j=0;j<16;j++)snap[r][j]=A[remaining11[r]][available16[j]];
    compare(snap,frozenR);for(auto&entry:A[16])assert(zero(entry));afterheavy=true;
   }
   Poly u,v;load(transcript,u);load(transcript,v);assert(!zero(u)&&!zero(A[p][c]));
   Poly cancel=multiply(u,A[t][c]);Poly vv=multiply(v,A[p][c]);assert(same(cancel,vv));
   // Compute EVERY column, including the supposedly cancelled column.
   for(int j=0;j<31;j++){Poly x=multiply(u,A[t][j]),y=multiply(v,A[p][j]);fmpq_mpoly_sub(A[t][j].x,x.x,y.x,qc);}
   assert(zero(A[t][c]));eliminations++;
  }
 }
 assert(afterquad&&afterheavy&&normalizations>0&&eliminations>0);
 Matrix QE;for(int i=0;i<10;i++)QE.push_back(A[piv[i].first]);compare(QE,matrix(saved+"/completed-quad10-echelon.txt",10,31));
 compare(A,matrix(saved+"/current-matrix.txt",26,31));
 for(auto[r,c]:piv)assert(!zero(A[r][c]));
 // Saved R11 relative row4 is G6 and is zero BEFORE the final cubic pivot.
 for(auto&p:frozenR[4])assert(zero(p));for(auto&p:A[16])assert(zero(p));
 ifstream remainders(saved+"/six-remainders.txt");string status;int count;remainders>>status>>count;assert(status=="all_six_zero_exact"&&count==15);
 for(int j=0;j<count;j++){int col;remainders>>col;Poly p;load(remainders,p);assert(col>=0&&col<31&&zero(p)&&zero(A[16][col]));}string tail;assert(!(remainders>>tail));
 cout<<"{\"all_passed\":true,\"all_saved_polynomial_row_operations_replayed\":true,\"eliminations\":"<<eliminations<<",\"normalizations\":"<<normalizations<<",\"nonconstant_normalization_factors\":"<<nonconstant_normalizations<<",\"all15_generic_pivots_nonzero\":true,\"all16x21_and11x16_checkpoints_equal\":true,\"whole26x31_final_checkpoint_equal\":true,\"G6_zero_before_final_cubic_pivot\":true,\"all15_saved_remainders_zero\":true,\"no_new_pivot_choice_gcd_rank_or_point\":true,\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}"<<endl;
 }fmpq_mpoly_ctx_clear(qc);flint_cleanup();return 0;}
