// Exact107x107 QQ solve and ALL-coefficient identity, no source evaluation.
#define main unused_frozen_source26_main
#include "q_new_family_v_zero_degree26_restriction_20261104_v63.cpp"
#undef main
#include <flint/fmpq_mat.h>
P readprimitive(istream&in){P p;string content;long n;in>>content>>n;assert(fmpq_set_str(p.v->content,content.c_str(),10)==0);fmpq_canonicalise(p.v->content);fmpz_t c;fmpz_init(c);ulong e[4];for(long k=0;k<n;k++){string z;in>>z;for(auto&x:e)in>>x;assert(fmpz_set_str(c,z.c_str(),10)==0);fmpz_mpoly_push_term_fmpz_ui(p.v->zpoly,c,e,ctx->zctx);}fmpz_mpoly_sort_terms(p.v->zpoly,ctx->zctx);fmpz_mpoly_combine_like_terms(p.v->zpoly,ctx->zctx);fmpq_mpoly_reduce(p.v,ctx);fmpz_clear(c);return p;}
int main(int argc,char**argv){assert(argc==3);fmpq_mpoly_ctx_init(ctx,4,ORD_DEGLEX);int result=0;auto start=chrono::steady_clock::now();{
 try{
  ifstream in(argv[1]);int n,rank;in>>n>>rank;assert(n==108&&rank==107);vector<P>cols;for(int i=0;i<n;i++)cols.push_back(readprimitive(in));vector<array<ulong,4>>rows(rank);for(auto&e:rows)for(auto&x:e)in>>x;vector<int>piv(rank);for(auto&x:piv)in>>x;int free;in>>free;string extra;assert(in&&!(in>>extra));assert(free>=0&&free<n);
  fmpq_mat_t A,b,x,check;fmpq_mat_init(A,rank,rank);fmpq_mat_init(b,rank,1);fmpq_mat_init(x,rank,1);fmpq_mat_init(check,rank,1);
  for(int i=0;i<rank;i++){for(int j=0;j<rank;j++)fmpq_mpoly_get_coeff_fmpq_ui(fmpq_mat_entry(A,i,j),cols[piv[j]].v,rows[i].data(),ctx);fmpq_mpoly_get_coeff_fmpq_ui(fmpq_mat_entry(b,i,0),cols[free].v,rows[i].data(),ctx);fmpq_neg(fmpq_mat_entry(b,i,0),fmpq_mat_entry(b,i,0));}
  assert(fmpq_mat_solve(x,A,b)==1);fmpq_mat_mul(check,A,x);assert(fmpq_mat_equal(check,b));
  fmpz_t den,g,t;fmpz_init(den);fmpz_init(g);fmpz_init(t);fmpz_one(den);for(int i=0;i<rank;i++)fmpz_lcm(den,den,fmpq_denref(fmpq_mat_entry(x,i,0)));
  vector<string>kernel(n,"0");kernel[free]=zstr(den);fmpz_set(g,den);
  for(int i=0;i<rank;i++){fmpz_divexact(t,den,fmpq_denref(fmpq_mat_entry(x,i,0)));fmpz_mul(t,t,fmpq_numref(fmpq_mat_entry(x,i,0)));kernel[piv[i]]=zstr(t);fmpz_gcd(g,g,t);}
  bool negative=false;for(auto&q:kernel){fmpz_set_str(t,q.c_str(),10);if(!fmpz_is_zero(t)){negative=fmpz_sgn(t)<0;break;}}
  for(auto&q:kernel){fmpz_set_str(t,q.c_str(),10);fmpz_divexact(t,t,g);if(negative)fmpz_neg(t,t);q=zstr(t);}
  P total;long terms=0;for(int i=0;i<n;i++){terms+=fmpq_mpoly_length(cols[i].v,ctx);if(kernel[i]!="0")total=add(total,scaleq(cols[i],kernel[i]));}bool allzero=zero(total);
  ofstream out(argv[2]);out<<"{\"all_passed\":"<<(allzero?"true":"false")<<",\"exact_QQ_solve107\":true,\"all107_selected_equations_verified\":true,\"all_source_coefficients_zero\":"<<(allzero?"true":"false")<<",\"source_terms_checked\":"<<terms<<",\"primitive_integer_kernel108\":[";for(int i=0;i<n;i++){if(i)out<<',';out<<'\"'<<kernel[i]<<'\"';}out<<"],\"selected_matrix_QQ\":[";for(int i=0;i<rank;i++){if(i)out<<',';out<<'[';for(int j=0;j<rank;j++){if(j)out<<',';out<<'\"'<<qstr(fmpq_mat_entry(A,i,j))<<'\"';}out<<']';}out<<"],\"rhs_QQ\":[";for(int i=0;i<rank;i++){if(i)out<<',';out<<'\"'<<qstr(fmpq_mat_entry(b,i,0))<<'\"';}out<<"],\"solution_QQ\":[";for(int i=0;i<rank;i++){if(i)out<<',';out<<'\"'<<qstr(fmpq_mat_entry(x,i,0))<<'\"';}out<<"],\"remaining_polynomial\":";writepoly(out,total);out<<",\"max_terms\":"<<maxterms<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";out.close();
  fmpq_mat_clear(A);fmpq_mat_clear(b);fmpq_mat_clear(x);fmpq_mat_clear(check);fmpz_clear(den);fmpz_clear(g);fmpz_clear(t);if(!allzero)result=2;cerr<<"EXACT_QQ28_IDENTITY "<<allzero<<endl;
 }catch(const exception&e){cerr<<"EXACT_QQ28_STOP "<<e.what()<<endl;result=2;}
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();return result;}
