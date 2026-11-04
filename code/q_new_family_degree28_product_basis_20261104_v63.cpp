// Same saved formal plane only: no source restriction and no point evaluation.
#define main unused_immutable_universal24_main
#include "q_new_family_full_degree24_universal_basis_20261104_v63.cpp"
#undef main
P mul28(const P&a,const P&b){P c;nmod_mpoly_mul(c.v,a.v,b.v,ctx);++multiplications;long n=nmod_mpoly_length(c.v,ctx);maxterms=max(maxterms,n);if(n>435)throw runtime_error("formal degree28 support435 cap");return c;}
P read28(istream&in,int expected){int d,n;in>>d>>n;assert(d==expected&&n==(d+1)*(d+2)/2);P p;ulong e[2];for(int i=0;i<=d;i++)for(int j=0;j<=d-i;j++){ulong c;in>>c;assert(c<PRIME);e[0]=i;e[1]=j;if(c)nmod_mpoly_set_coeff_ui_ui(p.v,c,e,ctx);}return p;}
int main(int argc,char**argv){assert(argc==3);nmod_mpoly_ctx_init(ctx,2,ORD_DEGLEX,PRIME);int status=0;auto start=chrono::steady_clock::now();{
 try{
  ifstream in(argv[1]);int count;ulong prime;in>>prime>>count;assert(prime==PRIME&&count==80);
  P F=read28(in,4);vector<P>p,q12,q14,v24;for(int i=0;i<4;i++)p.push_back(read28(in,8));for(int i=0;i<4;i++)q12.push_back(read28(in,12));P q16=read28(in,16);for(int i=0;i<3;i++)q14.push_back(read28(in,14));for(int i=0;i<67;i++)v24.push_back(read28(in,24));string extra;assert(in&&!(in>>extra));
  vector<P>columns;for(auto&x:v24)columns.push_back(mul28(F,x));for(auto&x:q12)for(int i=0;i<4;i++)for(int j=i;j<4;j++)columns.push_back(mul28(x,mul28(p[i],p[j])));for(auto&x:q12)columns.push_back(mul28(q16,x));for(int i=0;i<3;i++)for(int j=i;j<3;j++)columns.push_back(mul28(q14[i],q14[j]));assert(columns.size()==117);
  Rank rank;for(int i=0;i<117;i++)rank.push(coeffs(columns[i],28),i);assert(rank.rank()<=108);
  ofstream out(argv[2]);out<<"{\"all_passed\":true,\"status\":\""<<(rank.rank()==108?"complete_QQ_Inv28_basis_certified":"product_span_rank_gate_failed")<<"\",\"prime\":32003,\"character_dimension28\":108,\"same_fixed_formal_plane\":true,\"candidate_count\":117,\"rank\":"<<rank.rank()<<",\"selected_columns\":";ints(out,rank.selected);out<<",\"selected_coefficient_rows\":";ints(out,rank.pivots);out<<",\"columns117\":";polys(out,columns,28);out<<",\"max_terms\":"<<maxterms<<",\"products\":"<<multiplications<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";out.close();cerr<<"UNIVERSAL28_RANK "<<rank.rank()<<endl;
 }catch(const exception&e){cerr<<"BOUNDED_UNIVERSAL28_STOP "<<e.what()<<endl;status=2;}
 }nmod_mpoly_ctx_clear(ctx);flint_cleanup();return status;
}
