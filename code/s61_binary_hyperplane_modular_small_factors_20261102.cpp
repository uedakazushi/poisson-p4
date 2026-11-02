#include <flint/nmod_mpoly_factor.h>
#define main old_main
#include "s61_binary_hyperplane_modular_determinants_20261102.cpp"
#undef main
int main(int argc,char**argv){assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);assert(in&&out);prime=32003;nmod_mpoly_ctx_init(ctx,3,ORD_LEX,prime);{string s;int i,j;while(in>>s&&s=="HIT"){in>>i>>j;Poly p;read(in,p);nmod_mpoly_factor_t f;nmod_mpoly_factor_init(f,ctx);assert(nmod_mpoly_factor(f,p.p,ctx));assert(nmod_mpoly_factor_matches(p.p,f,ctx));out<<"HIT "<<i<<' '<<j<<" FACTORS "<<f->num<<'\n';for(int k=0;k<f->num;k++){Poly q;nmod_mpoly_factor_get_base(q.p,f,k,ctx);out<<"EXP "<<nmod_mpoly_factor_get_exp_si(f,k,ctx)<<' ';write(out,q);}out.flush();cerr<<"FACTORED "<<i<<' '<<j<<" COUNT "<<f->num<<'\n';nmod_mpoly_factor_clear(f,ctx);}assert(out.good());}nmod_mpoly_ctx_clear(ctx);flint_cleanup();}
