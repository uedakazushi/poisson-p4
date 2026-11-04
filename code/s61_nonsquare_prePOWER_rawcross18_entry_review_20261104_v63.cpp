// Independent DEGREVLEX62 source/11 saved values plus reversed QQ reproducts.
#define main independent_raw24_main_NOT_CALLED
#include "s61_nonsquare_prePOWER_raw8_zero_review_20261104_v63.cpp"
#undef main
void rebase14(L&v,long target){assert(target<=v.a);ulong e[3];long k=v.a-target;if(!k)return;for(long i=0;i<fmpq_mpoly_length(v.p,C);i++){fmpq_mpoly_get_term_exp_ui(e,v.p,i,C);e[0]+=k;fmpq_mpoly_set_term_exp_ui(v.p,i,e,C);}v.a=target;assert(fmpq_mpoly_is_canonical(v.p,C));}
void add14(L&out,L&x,L&y,bool sub){long k=min(x.a,y.a);rebase14(x,k);rebase14(y,k);out.a=k;if(sub)fmpq_mpoly_sub(out.p,x.p,y.p,C);else fmpq_mpoly_add(out.p,x.p,y.p,C);guard();}
bool equal14(const L&x,const L&y){assert(!y.a);if(fmpq_mpoly_length(x.p,C)!=fmpq_mpoly_length(y.p,C))return false;ulong e[3];fmpq_t a,b;fmpq_init(a);fmpq_init(b);bool ok=true;for(long i=0;i<fmpq_mpoly_length(x.p,C);i++){fmpq_mpoly_get_term_exp_ui(e,x.p,i,C);long k=long(e[0])+x.a;assert(k>=0);e[0]=k;fmpq_mpoly_get_term_coeff_fmpq(a,x.p,i,C);fmpq_mpoly_get_coeff_fmpq_ui(b,y.p,e,C);if(!fmpq_equal(a,b)){ok=false;break;}}fmpq_clear(a);fmpq_clear(b);return ok;}
int main(int argc,char**argv){assert(argc==3);fmpq_mpoly_ctx_init(C,3,ORD_DEGREVLEX);{
 ifstream f(argv[1]);int n,nr;f>>n>>nr;assert(n==76&&nr==11);vector<int>roots(nr);vector<unique_ptr<L>>saved;for(int i=0;i<nr;i++){f>>roots[i];auto v=make_unique<L>();polyread(f,*v);saved.push_back(move(v));}vector<long>uses(n);for(auto&u:uses)f>>u;vector<unique_ptr<L>>values(n);long ops=0,compared=0;
 for(int i=0;i<n;i++){int id;char op;f>>id>>op;assert(id==i);if(op=='x')continue;auto v=make_unique<L>();vector<int>args;if(op=='p'){long shift;string unit;f>>shift>>unit;polyread(f,*v);v->a=shift;fmpq_t q;fmpq_init(q);assert(!fmpq_set_str(q,unit.c_str(),10));fmpq_canonicalise(q);fmpq_mpoly_scalar_div_fmpq(v->p,v->p,q,C);fmpq_clear(q);}else{int a,b;f>>a>>b;assert(a<i&&b<i&&values[a]&&values[b]);args={a,b};if(op=='m')multiply(*v,*values[a],*values[b]);else{assert(op=='a'||op=='s');add14(*v,*values[a],*values[b],op=='s');}}
  for(int j=0;j<nr;j++)if(roots[j]==i){assert(equal14(*v,*saved[j]));compared++;}values[i]=move(v);ops++;guard();for(int j:args){assert(uses[j]);if(!--uses[j])values[j].reset();}if(!uses[i])values[i].reset();}
 assert(ops==62&&compared==11);for(auto u:uses)assert(!u);values.clear();string extra;assert(!(f>>extra));for(int i=0;i<3;i++)assert(fmpq_mpoly_is_zero(saved[i]->p,C));
 for(int i=1;i<4;i++){L a,b,d;multiply(a,*saved[7+i],*saved[3]);multiply(b,*saved[7],*saved[3+i]);add14(d,b,a,true);assert(fmpq_mpoly_is_zero(d.p,C));}
 ofstream o(argv[2]);o<<"{\"all_passed\":true,\"source_ops\":62,\"saved11_QQ_polynomials_match\":true,\"no_kappa18_evaluation\":true,\"all3_raw_cross_QQ_products_zero\":true,\"order\":\"DEGREVLEX\",\"peak_terms_with15000reserve\":"<<peak<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";
 }fmpq_mpoly_ctx_clear(C);flint_cleanup();return 0;}
