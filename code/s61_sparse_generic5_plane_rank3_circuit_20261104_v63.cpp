// ONE rank3 check. Primitive rows are ONLY for rank; adapter uses ORIGINAL W.
#define main unused_frozen_qpoly_main
#include "q_new_family_v_zero_degree26_restriction_20261104_v63.cpp"
#undef main
#include <flint/fmpz_mpoly_q.h>
#include <tuple>
using PM=vector<vector<P>>;
long ops=0,mn=0,md=0;
P rd(istream&in){P p;int n;in>>n;fmpq_t q;fmpq_init(q);ulong e[5];for(int i=0;i<n;i++){string c;in>>c;for(auto&z:e)in>>z;assert(fmpq_set_str(q,c.c_str(),10)==0);fmpq_canonicalise(q);fmpq_mpoly_set_coeff_fmpq_ui(p.v,q,e,ctx);}fmpq_clear(q);return p;}
void pp(ostream&f,const P&p){f<<"{\"content\":\""<<qstr(p.v->content)<<"\",\"terms\":[";fmpz_t c;fmpz_init(c);ulong e[5];for(long i=0;i<fmpq_mpoly_length(p.v,ctx);i++){if(i)f<<',';fmpz_mpoly_get_term_coeff_fmpz(c,p.v->zpoly,i,ctx->zctx);fmpq_mpoly_get_term_exp_ui(e,p.v,i,ctx);f<<"[\""<<zstr(c)<<"\"";for(auto z:e)f<<','<<z;f<<']';}f<<"]}";fmpz_clear(c);}
struct F{fmpz_mpoly_q_t v;F(){fmpz_mpoly_q_init(v,ctx->zctx);}F(const F&a):F(){fmpz_mpoly_q_set(v,a.v,ctx->zctx);}F(F&&a):F(){fmpz_mpoly_q_swap(v,a.v,ctx->zctx);}F&operator=(F a){fmpz_mpoly_q_swap(v,a.v,ctx->zctx);return *this;}~F(){fmpz_mpoly_q_clear(v,ctx->zctx);}};
void guardf(const F&a){mn=max(mn,fmpz_mpoly_length(fmpz_mpoly_q_numref(a.v),ctx->zctx));md=max(md,fmpz_mpoly_length(fmpz_mpoly_q_denref(a.v),ctx->zctx));if(mn>150000||md>150000)throw runtime_error("150k numerator/denominator cap");}
F cv(const P&p){F a;fmpz_mpoly_scalar_mul_fmpz(fmpz_mpoly_q_numref(a.v),p.v->zpoly,fmpq_numref(p.v->content),ctx->zctx);fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(a.v),fmpq_denref(p.v->content),ctx->zctx);fmpz_mpoly_q_canonicalise(a.v,ctx->zctx);guardf(a);return a;}
bool zz(const F&a){return fmpz_mpoly_q_is_zero(a.v,ctx->zctx);}
F fm(const F&a,const F&b){F c;fmpz_mpoly_q_mul(c.v,a.v,b.v,ctx->zctx);++ops;guardf(c);return c;}
F fs(const F&a,const F&b){F c;fmpz_mpoly_q_sub(c.v,a.v,b.v,ctx->zctx);++ops;guardf(c);return c;}
F fd(const F&a,const F&b){assert(!zz(b));F c;fmpz_mpoly_q_div(c.v,a.v,b.v,ctx->zctx);++ops;guardf(c);assert(fmpz_mpoly_q_equal(fm(c,b).v,a.v,ctx->zctx));return c;}
long sz(const F&a){return fmpz_mpoly_length(fmpz_mpoly_q_numref(a.v),ctx->zctx)+fmpz_mpoly_length(fmpz_mpoly_q_denref(a.v),ctx->zctx);}
void ip(ostream&f,const fmpz_mpoly_t a){f<<'[';fmpz_t c;fmpz_init(c);ulong e[5];for(long i=0;i<fmpz_mpoly_length(a,ctx->zctx);i++){if(i)f<<',';fmpz_mpoly_get_term_coeff_fmpz(c,a,i,ctx->zctx);fmpz_mpoly_get_term_exp_ui(e,a,i,ctx->zctx);f<<"[\""<<zstr(c)<<"\"";for(auto z:e)f<<','<<z;f<<']';}f<<']';fmpz_clear(c);}
void pf(ostream&f,const F&a){f<<"{\"numerator\":";ip(f,fmpz_mpoly_q_numref(a.v));f<<",\"denominator\":";ip(f,fmpz_mpoly_q_denref(a.v));f<<'}';}
void cap(const string&d){size_t n=0;for(auto&p:filesystem::directory_iterator(d))if(p.is_regular_file())n+=p.file_size();if(n>32*1024*1024)throw runtime_error("32MiB artifact cap");}
int main(int argc,char**argv){assert(argc==3);fmpq_mpoly_ctx_init(ctx,5,ORD_DEGLEX);auto start=chrono::steady_clock::now();int status=0;{
 try{
  string dir=argv[2];ifstream in(argv[1]);int n,m;in>>n>>m;assert(n==45&&m==3);PM original(45,vector<P>(3));for(auto&r:original)for(auto&p:r)p=rd(in);string extra;assert(in&&!(in>>extra));PM primitive=original;vector<P>factors(45);long divs=0;
  for(int i=0;i<45;i++){P g;bool found=false;for(auto&p:original[i])if(!zero(p)){if(!found){g=p;found=true;}else{P h;assert(fmpq_mpoly_gcd(h.v,g.v,p.v,ctx)==1);g=move(h);}if(fmpq_mpoly_is_fmpq(g.v,ctx))break;}
   if(found&&!fmpq_mpoly_is_fmpq(g.v,ctx)){factors[i]=g;for(int j=0;j<3;j++){P q;assert(fmpq_mpoly_divides(q.v,original[i][j].v,g.v,ctx)==1&&equal(mul(q,g),original[i][j]));primitive[i][j]=move(q);divs++;}}
   else fmpq_mpoly_one(factors[i].v,ctx);
  }
  {ofstream f(dir+"/row-primitive45.json");f<<"{\"all_passed\":true,\"factors45\":[";for(int i=0;i<45;i++){if(i)f<<',';pp(f,factors[i]);}f<<"],\"quotients45x3\":[";for(int i=0;i<45;i++){if(i)f<<',';f<<'[';for(int j=0;j<3;j++){if(j)f<<',';pp(f,primitive[i][j]);}f<<']';}f<<"],\"exact_divisions_remultiplied\":"<<divs<<",\"row_scaling_for_rank_only\":true}\n";f.close();cap(dir);}
  vector<vector<F>>A(45,vector<F>(3));for(int i=0;i<45;i++)for(int j=0;j<3;j++)A[i][j]=cv(primitive[i][j]);vector<int>rr(45),cc={0,1,2};for(int i=0;i<45;i++)rr[i]=i;vector<F>piv;int rank=0;
  for(int k=0;k<3;k++){int pi=-1,pj=-1;bool first=true;tuple<int,long,int,int>best;for(int i=k;i<45;i++)for(int j=k;j<3;j++)if(!zz(A[i][j])){bool mono=fmpz_mpoly_length(fmpz_mpoly_q_numref(A[i][j].v),ctx->zctx)==1&&fmpz_mpoly_length(fmpz_mpoly_q_denref(A[i][j].v),ctx->zctx)==1;auto sc=make_tuple(mono?0:1,sz(A[i][j]),i,j);if(first||sc<best){first=false;best=sc;pi=i;pj=j;}}if(pi<0)break;
   if(pi!=k){swap(A[pi],A[k]);swap(rr[pi],rr[k]);}if(pj!=k){for(auto&r:A)swap(r[pj],r[k]);swap(cc[pj],cc[k]);}F p=A[k][k];piv.push_back(p);for(int j=k;j<3;j++)A[k][j]=fd(A[k][j],p);vector<F>multipliers(45);for(int i=k+1;i<45;i++)if(!zz(A[i][k])){F c=A[i][k];multipliers[i]=c;for(int j=k+1;j<3;j++)A[i][j]=fs(A[i][j],fm(c,A[k][j]));A[i][k]=F();}rank++;
   {ofstream f(dir+"/field-op-"+to_string(k)+".json");f<<"{\"stage\":"<<k<<",\"swap_row_with\":"<<pi<<",\"swap_column_with\":"<<pj<<",\"pivot\":";pf(f,p);f<<",\"normalized_row3\":[";for(int j=0;j<3;j++){if(j)f<<',';pf(f,A[k][j]);}f<<"],\"multipliers45\":[";for(int i=0;i<45;i++){if(i)f<<',';pf(f,multipliers[i]);}f<<"]}\n";f.close();cap(dir);}
   cerr<<"PLANE_FIELD_PIVOT "<<rank<<" NUMDEN "<<mn<<' '<<md<<endl;
  }
  vector<int>rows(rr.begin(),rr.begin()+rank);int sign=1;for(int i=0;i<3;i++)for(int j=i+1;j<3;j++)if(cc[i]>cc[j])sign=-sign;
  {ofstream f(dir+"/rank-certificate.json");f<<"{\"all_passed\":true,\"rank\":"<<rank<<",\"selected_original_q_rows\":";ints(f,rows);f<<",\"column_order\":";ints(f,cc);f<<",\"column_permutation_sign\":"<<sign<<",\"nonzero_pivot_factors\":[";for(size_t i=0;i<piv.size();i++){if(i)f<<',';pf(f,piv[i]);}f<<"],\"selected_row_factors\":[";for(size_t i=0;i<rows.size();i++){if(i)f<<',';pp(f,factors[rows[i]]);}f<<"],\"actual_W_minor_is_sign_times_rowfactors_times_pivots\":true,\"original_W_required_for_adapter\":true,\"max_num_terms\":"<<mn<<",\"max_den_terms\":"<<md<<",\"field_operations\":"<<ops<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";f.close();cap(dir);}
 }catch(const exception&e){cerr<<"BOUNDED_PLANE_RANK3_STOP "<<e.what()<<endl;status=2;}
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();return status;}
