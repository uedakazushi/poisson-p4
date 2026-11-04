// ONE canonical generic-field solve from the frozen11-unit prefix.
#define main unused_frozen_qpoly_main
#include "q_new_family_v_zero_degree26_restriction_20261104_v63.cpp"
#undef main
#include <flint/fmpz_mpoly_q.h>
#include <tuple>
using PM=vector<vector<P>>;
long fieldops=0,maxnum=0,maxden=0;
P readp(istream&in){P p;long n;in>>n;fmpq_t c;fmpq_init(c);ulong e[5];for(long i=0;i<n;i++){string z;in>>z;for(auto&x:e)in>>x;assert(fmpq_set_str(c,z.c_str(),10)==0);fmpq_canonicalise(c);fmpq_mpoly_set_coeff_fmpq_ui(p.v,c,e,ctx);}fmpq_clear(c);return p;}
void pp(ostream&f,const P&p){f<<"{\"content\":\""<<qstr(p.v->content)<<"\",\"terms\":[";fmpz_t c;fmpz_init(c);ulong e[5];for(long i=0;i<fmpq_mpoly_length(p.v,ctx);i++){if(i)f<<',';fmpz_mpoly_get_term_coeff_fmpz(c,p.v->zpoly,i,ctx->zctx);fmpq_mpoly_get_term_exp_ui(e,p.v,i,ctx);f<<"[\""<<zstr(c)<<"\"";for(auto x:e)f<<','<<x;f<<']';}f<<"]}";fmpz_clear(c);}
struct F{fmpz_mpoly_q_t v;F(){fmpz_mpoly_q_init(v,ctx->zctx);}F(const F&a):F(){fmpz_mpoly_q_set(v,a.v,ctx->zctx);}F(F&&a):F(){fmpz_mpoly_q_swap(v,a.v,ctx->zctx);}F&operator=(F a){fmpz_mpoly_q_swap(v,a.v,ctx->zctx);return *this;}~F(){fmpz_mpoly_q_clear(v,ctx->zctx);}};
void fc(const F&a){maxnum=max(maxnum,fmpz_mpoly_length(fmpz_mpoly_q_numref(a.v),ctx->zctx));maxden=max(maxden,fmpz_mpoly_length(fmpz_mpoly_q_denref(a.v),ctx->zctx));if(maxnum>150000||maxden>150000)throw runtime_error("field numerator/denominator150k cap");}
F fp(const P&p){F a;fmpz_mpoly_scalar_mul_fmpz(fmpz_mpoly_q_numref(a.v),p.v->zpoly,fmpq_numref(p.v->content),ctx->zctx);fmpz_mpoly_set_fmpz(fmpz_mpoly_q_denref(a.v),fmpq_denref(p.v->content),ctx->zctx);fmpz_mpoly_q_canonicalise(a.v,ctx->zctx);fc(a);return a;}
bool fzero(const F&a){return fmpz_mpoly_q_is_zero(a.v,ctx->zctx);}
F fa(const F&a,const F&b){F c;fmpz_mpoly_q_add(c.v,a.v,b.v,ctx->zctx);++fieldops;fc(c);return c;}
F fs(const F&a,const F&b){F c;fmpz_mpoly_q_sub(c.v,a.v,b.v,ctx->zctx);++fieldops;fc(c);return c;}
F fm(const F&a,const F&b){F c;fmpz_mpoly_q_mul(c.v,a.v,b.v,ctx->zctx);++fieldops;fc(c);return c;}
F fd(const F&a,const F&b){assert(!fzero(b));F c;fmpz_mpoly_q_div(c.v,a.v,b.v,ctx->zctx);++fieldops;fc(c);assert(fmpz_mpoly_q_equal(fm(c,b).v,a.v,ctx->zctx));return c;}
F neg(const F&a){F b;fmpz_mpoly_q_neg(b.v,a.v,ctx->zctx);return b;}
F onef(){F a;fmpz_mpoly_q_one(a.v,ctx->zctx);return a;}
using FM=vector<vector<F>>;
void zp(ostream&f,const fmpz_mpoly_t p){f<<'[';fmpz_t c;fmpz_init(c);ulong e[5];for(long i=0;i<fmpz_mpoly_length(p,ctx->zctx);i++){if(i)f<<',';fmpz_mpoly_get_term_coeff_fmpz(c,p,i,ctx->zctx);fmpz_mpoly_get_term_exp_ui(e,p,i,ctx->zctx);f<<"[\""<<zstr(c)<<"\"";for(auto x:e)f<<','<<x;f<<']';}f<<']';fmpz_clear(c);}
void pf(ostream&f,const F&a){f<<"{\"numerator\":";zp(f,fmpz_mpoly_q_numref(a.v));f<<",\"denominator\":";zp(f,fmpz_mpoly_q_denref(a.v));f<<'}';}
void mf(ostream&f,const FM&A){f<<'[';for(size_t i=0;i<A.size();i++){if(i)f<<',';f<<'[';for(size_t j=0;j<A[i].size();j++){if(j)f<<',';pf(f,A[i][j]);}f<<']';}f<<']';}
void capdir(const string&d){size_t b=0;for(auto&p:filesystem::directory_iterator(d))if(p.is_regular_file())b+=p.file_size();if(b>32*1024*1024)throw runtime_error("artifact32MiB cap");}
long sizef(const F&a){return fmpz_mpoly_length(fmpz_mpoly_q_numref(a.v),ctx->zctx)+fmpz_mpoly_length(fmpz_mpoly_q_denref(a.v),ctx->zctx);}
int main(int argc,char**argv){assert(argc==3);fmpq_mpoly_ctx_init(ctx,5,ORD_DEGLEX);int status=0;auto start=chrono::steady_clock::now();{
 try{
  string dir=argv[2];ifstream in(argv[1]);int vars;in>>vars;assert(vars==5);PM prefix(51,vector<P>(20)),L(51,vector<P>(20)),N(45,vector<P>(20));for(auto*A:{&prefix,&L,&N})for(auto&r:*A)for(auto&p:r)p=readp(in);vector<P>center,d;for(auto*v:{&center,&d})for(int i=0;i<45;i++)v->push_back(readp(in));vector<int>co(20);for(auto&i:co)in>>i;P oldD=readp(in);PM oldT(20,vector<P>(2));for(auto&r:oldT)for(auto&p:r)p=readp(in);vector<int>oldfree(2);for(auto&i:oldfree)in>>i;string extra;assert(in&&!(in>>extra));
  PM primitive(40,vector<P>(9));vector<P>gcds(40);long gcddiv=0;
  for(int i=0;i<40;i++){
   P g;bool initialized=false;for(int j=0;j<9;j++){primitive[i][j]=prefix[11+i][11+j];if(zero(primitive[i][j]))continue;if(!initialized){g=primitive[i][j];initialized=true;}else{P h;assert(fmpq_mpoly_gcd(h.v,g.v,primitive[i][j].v,ctx)==1);g=move(h);}if(fmpq_mpoly_is_fmpq(g.v,ctx))break;}
   if(initialized&&!fmpq_mpoly_is_fmpq(g.v,ctx)){gcds[i]=g;for(int j=0;j<9;j++){P z=prefix[11+i][11+j],v;assert(fmpq_mpoly_divides(v.v,z.v,g.v,ctx)==1&&equal(mul(v,g),z));primitive[i][j]=move(v);gcddiv++;}}
   else for(int j=0;j<9;j++)primitive[i][j]=prefix[11+i][11+j];
  }
  {ofstream f(dir+"/row-primitive.json");f<<"{\"all_passed\":true,\"row_gcds40\":[";for(int i=0;i<40;i++){if(i)f<<',';pp(f,gcds[i]);}f<<"],\"primitive40x9\":[";for(int i=0;i<40;i++){if(i)f<<',';f<<'[';for(int j=0;j<9;j++){if(j)f<<',';pp(f,primitive[i][j]);}f<<']';}f<<"],\"literal_exact_division_and_remultiplication_count\":"<<gcddiv<<"}\n";f.close();capdir(dir);}
  FM A(40,vector<F>(9));for(int i=0;i<40;i++)for(int j=0;j<9;j++)A[i][j]=fp(primitive[i][j]);vector<int>tc(9),rr(40);for(int i=0;i<9;i++)tc[i]=i;for(int i=0;i<40;i++)rr[i]=i;vector<F>pivots;int rank=0;
  for(int k=0;k<7;k++){
   bool first=true;tuple<int,long,long,int,int>best;int pi=-1,pj=-1;
   for(int i=k;i<40;i++)for(int j=k;j<9;j++)if(!fzero(A[i][j])){long nn=fmpz_mpoly_length(fmpz_mpoly_q_numref(A[i][j].v),ctx->zctx),dd=fmpz_mpoly_length(fmpz_mpoly_q_denref(A[i][j].v),ctx->zctx),deg=fmpz_mpoly_total_degree_si(fmpz_mpoly_q_numref(A[i][j].v),ctx->zctx)+fmpz_mpoly_total_degree_si(fmpz_mpoly_q_denref(A[i][j].v),ctx->zctx);auto score=make_tuple(nn==1&&dd==1?0:1,nn+dd,deg,i,j);if(first||score<best){first=false;best=score;pi=i;pj=j;}}
   if(pi<0)break;if(pi!=k){swap(A[pi],A[k]);swap(rr[pi],rr[k]);}if(pj!=k){for(auto&r:A)swap(r[pj],r[k]);swap(tc[pj],tc[k]);}F pivot=A[k][k];pivots.push_back(pivot);for(int j=k;j<9;j++)A[k][j]=fd(A[k][j],pivot);assert(fmpz_mpoly_q_is_one(A[k][k].v,ctx->zctx));
   vector<F>normalized=A[k],multipliers(40);for(int i=0;i<40;i++)if(i!=k&&!fzero(A[i][k])){F c=A[i][k];multipliers[i]=c;for(int j=k+1;j<9;j++)A[i][j]=fs(A[i][j],fm(c,A[k][j]));A[i][k]=F();}rank++;
   {ofstream f(dir+"/field-op-"+to_string(k)+".json");f<<"{\"stage\":"<<k<<",\"swap_row_with\":"<<pi<<",\"swap_column_with\":"<<pj<<",\"pivot\":";pf(f,pivot);f<<",\"normalized_row9\":[";for(int j=0;j<9;j++){if(j)f<<',';pf(f,normalized[j]);}f<<"],\"elimination_multipliers40\":[";for(int i=0;i<40;i++){if(i)f<<',';pf(f,multipliers[i]);}f<<"]}\n";f.close();capdir(dir);}
   {ofstream f(dir+"/field-prefix-checkpoint.json");f<<"{\"completed_pivots\":"<<rank<<",\"row_order\":";ints(f,rr);f<<",\"column_order\":";ints(f,tc);f<<",\"rref40x9\":";mf(f,A);f<<",\"pivot_values\":[";for(size_t i=0;i<pivots.size();i++){if(i)f<<',';pf(f,pivots[i]);}f<<"]}\n";f.close();capdir(dir);}cerr<<"GENERIC_FIELD_PIVOT "<<rank<<" NUMDEN "<<maxnum<<' '<<maxden<<endl;
  }
  assert(rank==7);for(int i=7;i<40;i++)for(int j=7;j<9;j++)assert(fzero(A[i][j]));FM tail(9,vector<F>(2));for(int i=0;i<7;i++)for(int j=0;j<2;j++)tail[tc[i]][j]=neg(A[i][7+j]);for(int j=0;j<2;j++)tail[tc[7+j]][j]=onef();FM T(20,vector<F>(2));for(int i=0;i<9;i++)T[co[11+i]]=tail[i];for(int i=0;i<11;i++)for(int j=0;j<2;j++)for(int h=0;h<9;h++)if(!zero(prefix[i][11+h])&&!fzero(tail[h][j]))T[co[i]][j]=fs(T[co[i]][j],fm(fp(prefix[i][11+h]),tail[h][j]));vector<int>free={co[11+tc[7]],co[11+tc[8]]};
  for(int i=0;i<51;i++)for(int j=0;j<2;j++){F z;for(int h=0;h<20;h++)if(!zero(L[i][h])&&!fzero(T[h][j]))z=fa(z,fm(fp(L[i][h]),T[h][j]));assert(fzero(z));}for(int i=0;i<2;i++)for(int j=0;j<2;j++)assert(fmpz_mpoly_q_equal(T[free[i]][j].v,(i==j?onef():F()).v,ctx->zctx));
  {ofstream f(dir+"/generic-field-kernel.json");f<<"{\"all_passed\":true,\"tail_rank\":7,\"full_L_rank\":18,\"free_original_columns\":";ints(f,free);f<<",\"T20x2_rational\":";mf(f,T);f<<",\"tail9x2_rational\":";mf(f,tail);f<<",\"all51_kernel_and_freeI2_verified\":true,\"coefficient_field\":\"QQ(H,L,s0,q0,t0)\",\"not_FULLU\":true}\n";f.close();capdir(dir);}cerr<<"NEW_KERNEL_ALL51_PASS"<<endl;
  for(int i=0;i<20;i++)for(int j=0;j<2;j++){F x,y;for(int h=0;h<2;h++){x=fa(x,fm(T[i][h],fp(oldT[free[h]][j])));y=fa(y,fm(fp(oldT[i][h]),T[oldfree[h]][j]));}assert(fmpz_mpoly_q_equal(x.v,fp(oldT[i][j]).v,ctx->zctx));assert(fmpz_mpoly_q_equal(y.v,fm(fp(oldD),T[i][j]).v,ctx->zctx));}
  {ofstream f(dir+"/old-kernel-span.json");f<<"{\"all_passed\":true,\"both_directions_all20x2_exact_crossproducts\":true,\"same_source_field_kernel\":true}\n";f.close();capdir(dir);}cerr<<"OLD_NEW_SPAN_PASS"<<endl;
  FM V(45,vector<F>(3));long vm=0;for(int i=0;i<45;i++){for(int j=0;j<2;j++)for(int h=0;h<20;h++)if(!zero(N[i][h])&&!fzero(T[h][j]))V[i][j]=fa(V[i][j],fm(fp(N[i][h]),T[h][j]));V[i][2]=fp(d[i]);for(auto&p:V[i])vm=max(vm,sizef(p));}
  {ofstream f(dir+"/plane-rational45x3.json");f<<"{\"all_passed\":true,\"V45x3\":";mf(f,V);f<<",\"center45\":[";for(int i=0;i<45;i++){if(i)f<<',';pp(f,center[i]);}f<<"],\"columns\":[\"N*T0\",\"N*T1\",\"Q52-center\"],\"max_numden_support_sum\":"<<vm<<"}\n";f.close();capdir(dir);}
  bool planedone=false;
  if(vm<=300){FM W=V;vector<int>pr(45),pc={0,1,2};for(int i=0;i<45;i++)pr[i]=i;vector<F>p3;for(int k=0;k<3;k++){bool first=true;long best=0;int pi=-1,pj=-1;for(int i=k;i<45;i++)for(int j=k;j<3;j++)if(!fzero(W[i][j])&&(first||sizef(W[i][j])<best)){first=false;best=sizef(W[i][j]);pi=i;pj=j;}assert(pi>=0);if(pi!=k){swap(W[pi],W[k]);swap(pr[pi],pr[k]);}if(pj!=k){for(auto&r:W)swap(r[pj],r[k]);swap(pc[pj],pc[k]);}F p=W[k][k];p3.push_back(p);for(int j=k;j<3;j++)W[k][j]=fd(W[k][j],p);for(int i=k+1;i<45;i++)if(!fzero(W[i][k])){F c=W[i][k];for(int j=k+1;j<3;j++)W[i][j]=fs(W[i][j],fm(c,W[k][j]));W[i][k]=F();}}
   pr.resize(3);{ofstream f(dir+"/plane-rank3-factors.json");f<<"{\"all_passed\":true,\"rank\":3,\"selected_original_q_rows\":";ints(f,pr);f<<",\"pivot_column_order\":";ints(f,pc);f<<",\"nonzero_determinant_pivot_factors\":[";for(int i=0;i<3;i++){if(i)f<<',';pf(f,p3[i]);}f<<"],\"determinant_is_column_permutation_sign_times_product\":true,\"no_expanded3adjugate\":true}\n";f.close();capdir(dir);}planedone=true;
  }
  {ofstream f(dir+"/native-certificate.json");f<<"{\"all_passed\":true,\"tail_rank7\":true,\"all51_kernel\":true,\"old_span_verified\":true,\"plane_rank3_factored\":"<<(planedone?"true":"false")<<",\"field_operations\":"<<fieldops<<",\"max_numerator_terms\":"<<maxnum<<",\"max_denominator_terms\":"<<maxden<<",\"row_gcd_exact_divisions\":"<<gcddiv<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";f.close();capdir(dir);}
 }catch(const exception&e){cerr<<"BOUNDED_SMALL_TAIL_FIELD_STOP "<<e.what()<<endl;status=2;}
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();return status;}
