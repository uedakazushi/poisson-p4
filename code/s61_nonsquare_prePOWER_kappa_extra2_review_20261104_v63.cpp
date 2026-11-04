// Independent DEGREVLEX fixed238/240 source replay +four-product QQ identity.
// Numeric utility is our previous independent checker, never a producer header.
#define main independent_raw24_main_NOT_CALLED
#include "s61_nonsquare_prePOWER_raw8_zero_review_20261104_v63.cpp"
#undef main
void rebase_inplace(L&x,long target){assert(target<=x.a);if(target==x.a)return;ulong e[3];long k=x.a-target;for(long j=0;j<fmpq_mpoly_length(x.p,C);j++){fmpq_mpoly_get_term_exp_ui(e,x.p,j,C);e[0]+=k;fmpq_mpoly_set_term_exp_ui(x.p,j,e,C);}x.a=target;assert(fmpq_mpoly_is_canonical(x.p,C));}
void sum_inplace_basis(L&z,L&x,L&y,bool subtract){long a=min(x.a,y.a);rebase_inplace(x,a);rebase_inplace(y,a);z.a=a;if(subtract)fmpq_mpoly_sub(z.p,x.p,y.p,C);else fmpq_mpoly_add(z.p,x.p,y.p,C);guard();}
bool equal_saved(const L&x,const L&y){assert(!y.a);if(fmpq_mpoly_length(x.p,C)!=fmpq_mpoly_length(y.p,C))return false;ulong e[3];fmpq_t a,b;fmpq_init(a);fmpq_init(b);bool ok=true;for(long j=0;j<fmpq_mpoly_length(x.p,C);j++){fmpq_mpoly_get_term_exp_ui(e,x.p,j,C);long aa=long(e[0])+x.a;assert(aa>=0);e[0]=aa;fmpq_mpoly_get_term_coeff_fmpq(a,x.p,j,C);fmpq_mpoly_get_coeff_fmpq_ui(b,y.p,e,C);if(!fmpq_equal(a,b)){ok=false;break;}}fmpq_clear(a);fmpq_clear(b);return ok;}
int main(int argc,char**argv){assert(argc==3);fmpq_mpoly_ctx_init(C,3,ORD_DEGREVLEX);{
 ifstream f(argv[1]);int n,nr;f>>n>>nr;assert(n==240&&nr==9);vector<long>ids(nr);vector<unique_ptr<L>>saved;for(int i=0;i<nr;i++){f>>ids[i];auto x=make_unique<L>();polyread(f,*x);saved.push_back(move(x));}
 vector<long>uses(n);for(auto&v:uses)f>>v;vector<unique_ptr<L>>values(n);long checked=0;
 for(int i=0;i<n;i++){int id;char op;f>>id>>op;assert(id==i);if(op=='x')continue;auto x=make_unique<L>();vector<int>args;if(op=='p'){long shift;string unit;f>>shift>>unit;polyread(f,*x);x->a=shift;fmpq_t q;fmpq_init(q);assert(!fmpq_set_str(q,unit.c_str(),10));fmpq_canonicalise(q);fmpq_mpoly_scalar_div_fmpq(x->p,x->p,q,C);fmpq_clear(q);}else{int a,b;f>>a>>b;assert(a<i&&b<i&&values[a]&&values[b]);args={a,b};if(op=='m')multiply(*x,*values[a],*values[b]);else{assert(op=='a'||op=='s');sum_inplace_basis(*x,*values[a],*values[b],op=='s');}}
  for(int j=0;j<nr;j++)if(ids[j]==i){assert(equal_saved(*x,*saved[j]));checked++;}values[i]=move(x);guard();for(int j:args){assert(uses[j]>0);if(!--uses[j])values[j].reset();}if(!uses[i])values[i].reset();
 }
 assert(checked==9);for(auto u:uses)assert(!u);string extra;assert(!(f>>extra));values.clear();L total;vector<int>order={2,0,3,1};for(int k:order){L product,result;multiply(product,*saved[1+k],*saved[5+k]);sum_inplace_basis(result,total,product,false);fmpq_mpoly_set(total.p,result.p,C);total.a=result.a;guard();}assert(equal_saved(total,*saved[0]));assert(!fmpq_mpoly_is_zero(saved[0]->p,C));
 ofstream out(argv[2]);out<<"{\"all_passed\":true,\"fixed_source_operations\":238,\"all9_source_values_QQ_match\":true,\"four_products_order\":[2,0,3,1],\"literal_four_product_sum_equals_T\":true,\"T_nonzero\":true,\"T_terms\":"<<fmpq_mpoly_length(saved[0]->p,C)<<",\"different_monomial_order\":\"DEGREVLEX\",\"peak_terms_with15000reserve\":"<<peak<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";
 }fmpq_mpoly_ctx_clear(C);flint_cleanup();return 0;}
