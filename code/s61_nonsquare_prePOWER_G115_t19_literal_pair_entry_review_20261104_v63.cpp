// Independent ordinary QQ products, DEGREVLEX; no integer-lattice worker reuse.
#define main independent_raw24_main_NOT_CALLED
#include "s61_nonsquare_prePOWER_raw8_zero_review_20261104_v63.cpp"
#undef main
#include <map>
void rb(L&v,long target){assert(target<=v.a);if(target==v.a)return;ulong e[3];long k=v.a-target;for(long i=0;i<fmpq_mpoly_length(v.p,C);i++){fmpq_mpoly_get_term_exp_ui(e,v.p,i,C);e[0]+=k;fmpq_mpoly_set_term_exp_ui(v.p,i,e,C);}v.a=target;assert(fmpq_mpoly_is_canonical(v.p,C));}
void sum_into(L&out,L&x,L&y,bool sub){long k=min(x.a,y.a);rb(x,k);rb(y,k);out.a=k;if(sub)fmpq_mpoly_sub(out.p,x.p,y.p,C);else fmpq_mpoly_add(out.p,x.p,y.p,C);guard();}
bool eqsaved(const L&x,const L&y){assert(!y.a);if(fmpq_mpoly_length(x.p,C)!=fmpq_mpoly_length(y.p,C))return false;fmpq_t a,b;fmpq_init(a);fmpq_init(b);ulong e[3];bool ok=true;for(long i=0;i<fmpq_mpoly_length(x.p,C);i++){fmpq_mpoly_get_term_exp_ui(e,x.p,i,C);long k=long(e[0])+x.a;assert(k>=0);e[0]=k;fmpq_mpoly_get_term_coeff_fmpq(a,x.p,i,C);fmpq_mpoly_get_coeff_fmpq_ui(b,y.p,e,C);if(!fmpq_equal(a,b)){ok=false;break;}}fmpq_clear(a);fmpq_clear(b);return ok;}
struct Op{char kind='x';long a=-1,b=-1,shift=0,upper=0;string den;};
struct Event{int id,a,b,sign;};
int main(int argc,char**argv){assert(argc==3);fmpq_mpoly_ctx_init(C,3,ORD_DEGREVLEX);{
 ifstream in(argv[1]);int n;in>>n;assert(n==1357);vector<long>uses(n);for(auto&x:uses)in>>x;vector<Op>ops(n);vector<unique_ptr<L>>v(n);
 for(int i=0;i<n;i++){int id;in>>id>>ops[i].kind;assert(id==i);auto&o=ops[i];if(o.kind=='x')continue;in>>o.upper>>o.den;if(o.kind=='p'){string unit;in>>o.shift>>unit;v[i]=make_unique<L>();polyread(in,*v[i]);v[i]->a=o.shift;fmpq_t q;fmpq_init(q);assert(!fmpq_set_str(q,unit.c_str(),10));fmpq_canonicalise(q);fmpq_mpoly_scalar_div_fmpq(v[i]->p,v[i]->p,q,C);fmpq_clear(q);}else if(o.kind=='k'){in>>o.a>>o.shift;assert(o.a<i);}else{in>>o.a>>o.b;assert(o.a<i&&o.b<i);}}
 int ne;in>>ne;assert(ne==15);map<int,vector<Event>>events;for(int i=0;i<ne;i++){int at;Event e;in>>at>>e.id>>e.a>>e.b>>e.sign;assert(at==max(e.a,e.b));events[at].push_back(e);}L total;long checked=0,done=0;fmpz_t D;fmpz_init(D);
 for(int i=0;i<n;i++){auto&o=ops[i];if(o.kind=='x')continue;
  if(o.kind!='p'){assert(v[o.a]);v[i]=make_unique<L>();if(o.kind=='k'){fmpq_mpoly_set(v[i]->p,v[o.a]->p,C);v[i]->a=v[o.a]->a+o.shift;}else{assert(v[o.b]);if(o.kind=='m')multiply(*v[i],*v[o.a],*v[o.b]);else{assert(o.kind=='a'||o.kind=='s');sum_into(*v[i],*v[o.a],*v[o.b],o.kind=='s');}}}
  assert(fmpq_mpoly_length(v[i]->p,C)<=o.upper);assert(!fmpz_set_str(D,o.den.c_str(),10));assert(fmpz_divisible(D,fmpq_denref(v[i]->p->content)));checked++;guard();
  if(o.kind!='p'){assert(uses[o.a]);if(!--uses[o.a])v[o.a].reset();if(o.kind!='k'){assert(uses[o.b]);if(!--uses[o.b])v[o.b].reset();}}
  for(auto&e:events[i]){{assert(v[e.a]&&v[e.b]);L product;multiply(product,*v[e.a],*v[e.b]);if(e.sign<0)fmpq_mpoly_neg(product.p,product.p,C);sum_into(total,total,product,false);done++;}if(!--uses[e.a])v[e.a].reset();if(!--uses[e.b])v[e.b].reset();guard();}
  if(!uses[i])v[i].reset();
 }
 assert(checked==1325&&done==15);for(auto u:uses)assert(!u);v.clear();L expected;polyread(in,expected);string extra;assert(!(in>>extra));assert(eqsaved(total,expected)&&!fmpq_mpoly_is_zero(expected.p,C));assert(fmpq_mpoly_length(expected.p,C)==9047);
 ofstream out(argv[2]);out<<"{\"all_passed\":true,\"all1325_source_operations_and_denominator_bounds\":true,\"ordinary_QQ_laplace_products15\":15,\"saved_normalized_D6_all9047_QQ_coefficients_match\":true,\"normalized_D6_NONZERO\":true,\"order\":\"DEGREVLEX\",\"integer_lattice_producer_not_used\":true,\"peak_terms_with15000reserve\":"<<peak<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";fmpz_clear(D);
 }fmpq_mpoly_ctx_clear(C);flint_cleanup();return 0;}
