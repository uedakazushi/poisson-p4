// Exact symbolic minors only; no point, prime, or modular rank.
#define main unused_frozen_qpoly_main
#include "q_new_family_v_zero_degree26_restriction_20261104_v63.cpp"
#undef main
struct Step{int k,i,j;P pivot,oldij,oldik,oldkj,divisor,quotient;};
struct Attempt{int deleted;vector<int>rows,swaps;vector<Step>steps;P determinant;};
P readminorpoly(istream&in){P p;int n;in>>n;fmpq_t c;fmpq_init(c);ulong e[4];for(int i=0;i<n;i++){string q;in>>q;for(auto&x:e)in>>x;assert(fmpq_set_str(c,q.c_str(),10)==0);fmpq_canonicalise(c);fmpq_mpoly_set_coeff_fmpq_ui(p.v,c,e,ctx);}fmpq_clear(c);return p;}
P oneminor(){P p;fmpq_mpoly_one(p.v,ctx);return p;}
Attempt determinant(const vector<vector<P>>&K,const vector<int>&rows,int deleted){
 Attempt r;r.deleted=deleted;vector<vector<P>>A;for(size_t i=0;i<K.size();i++)if(rows[i]!=deleted){A.push_back(K[i]);r.rows.push_back(rows[i]);}assert(A.size()==10);P previous=oneminor();int sign=1;
 for(int k=0;k<9;k++){
  int pivot=k;while(pivot<10&&zero(A[pivot][k]))pivot++;
  if(pivot==10)return r;
  if(pivot!=k){swap(A[pivot],A[k]);r.swaps.push_back(k);r.swaps.push_back(pivot);sign=-sign;}
  P pp=A[k][k];assert(!zero(pp)&&!zero(previous));
  for(int i=k+1;i<10;i++){
   P left=A[i][k];for(int j=k+1;j<10;j++){
    P num=add(mul(pp,A[i][j]),scale(mul(left,A[k][j]),-1)),q;
    int ok=fmpq_mpoly_divides(q.v,num.v,previous.v,ctx);assert(ok==1&&equal(mul(q,previous),num));check(q);
    r.steps.push_back(Step{k,i,j,pp,A[i][j],left,A[k][j],previous,q});A[i][j]=move(q);
   }A[i][k]=P();
  }previous=move(pp);
 }r.determinant=scale(A[9][9],sign);return r;
}
int main(int argc,char**argv){assert(argc==3);fmpq_mpoly_ctx_init(ctx,4,ORD_DEGLEX);int status=0;auto start=chrono::steady_clock::now();{
 try{
  ifstream in(argv[1]);int n,m;in>>n>>m;assert(n==11&&m==10);vector<int>rows(n);for(auto&x:rows)in>>x;vector<vector<P>>K(n,vector<P>(m));for(auto&r:K)for(auto&p:r)p=readminorpoly(in);string extra;assert(in&&!(in>>extra));
  vector<int>order={7};for(int j:rows)if(j!=7)order.push_back(j);vector<Attempt>attempts;int found=-1;
  for(int j:order){attempts.push_back(determinant(K,rows,j));cerr<<"EXACT_MINOR_ROW "<<j<<" TERMS "<<fmpq_mpoly_length(attempts.back().determinant.v,ctx)<<endl;if(!zero(attempts.back().determinant)){found=j;break;}}
  ofstream out(argv[2]);out<<"{\"all_passed\":true,\"nonzero_deleted_row\":"<<found<<",\"attempts\":[";
  for(size_t a=0;a<attempts.size();a++){auto&r=attempts[a];if(a)out<<',';out<<"{\"deleted_original_row\":"<<r.deleted<<",\"selected_original_rows\":";ints(out,r.rows);out<<",\"row_swap_pairs_flat\":";ints(out,r.swaps);out<<",\"determinant\":";writepoly(out,r.determinant);out<<",\"exact_division_steps\":[";
   for(size_t k=0;k<r.steps.size();k++){auto&t=r.steps[k];if(k)out<<',';out<<"{\"stage\":"<<t.k<<",\"row\":"<<t.i<<",\"column\":"<<t.j<<",\"pivot\":";writepoly(out,t.pivot);out<<",\"old_ij\":";writepoly(out,t.oldij);out<<",\"old_ik\":";writepoly(out,t.oldik);out<<",\"old_kj\":";writepoly(out,t.oldkj);out<<",\"divisor\":";writepoly(out,t.divisor);out<<",\"quotient\":";writepoly(out,t.quotient);out<<'}';}out<<"]}";
  }out<<"],\"all_Bareiss_divisions_exact_and_remultiplied\":true,\"no_point_prime_modular_rank\":true,\"max_terms\":"<<maxterms<<",\"products\":"<<products<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";out.close();if(filesystem::file_size(argv[2])>32*1024*1024)throw runtime_error("32MiB certificate cap");
 }catch(const exception&e){cerr<<"BOUNDED_SYMBOLIC_MINOR_STOP "<<e.what()<<endl;status=2;}
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();return status;}
