#include <flint/nmod_poly.h>
#include <algorithm>
#include <tuple>
vector<Rat> unitFactors;
nmod_poly_t unitAt2,unitAt3;
void evaluateCore(nmod_poly_t out,const nmod_mpoly_t in,ulong e){Rat temp;nmod_mpoly_evaluate_one_ui(temp.n,in,0,e,ctx);assert(nmod_mpoly_get_nmod_poly(out,temp.n,1,ctx));}
bool univariateUnit(const nmod_mpoly_t in,ulong e,const nmod_poly_t allowed){nmod_poly_t p,g,q,r;nmod_poly_init(p,prime);nmod_poly_init(g,prime);nmod_poly_init(q,prime);nmod_poly_init(r,prime);evaluateCore(p,in,e);bool ok=true;if(!nmod_poly_is_zero(p)){while(nmod_poly_degree(p)>0){nmod_poly_gcd(g,p,allowed);if(nmod_poly_degree(g)<=0){ok=false;break;}nmod_poly_divrem(q,r,p,g);assert(nmod_poly_is_zero(r));nmod_poly_swap(p,q);}}nmod_poly_clear(p);nmod_poly_clear(g);nmod_poly_clear(q);nmod_poly_clear(r);return ok;}
bool localizedUnit(const nmod_mpoly_t in){if(nmod_mpoly_degree_si(in,1,ctx)<=0)return true;if(!univariateUnit(in,2,unitAt2)||!univariateUnit(in,3,unitAt3))return false;Rat rem,g,tmp;nmod_mpoly_set(rem.n,in,ctx);bool changed=true;while(changed){changed=false;for(auto&f:unitFactors){assert(nmod_mpoly_gcd(g.n,rem.n,f.n,ctx));if(nmod_mpoly_total_degree_si(g.n,ctx)>0){assert(nmod_mpoly_divides(tmp.n,rem.n,g.n,ctx));nmod_mpoly_swap(tmp.n,rem.n,ctx);changed=true;}}}return nmod_mpoly_degree_si(rem.n,1,ctx)<=0;}
void loadUnitFactors(const char*path){ifstream in(path);int n;in>>n;unitFactors.resize(n);nmod_poly_init(unitAt2,prime);nmod_poly_init(unitAt3,prime);nmod_poly_one(unitAt2);nmod_poly_one(unitAt3);nmod_poly_t p,t;nmod_poly_init(p,prime);nmod_poly_init(t,prime);for(auto&f:unitFactors){readpoly(in,f.n);for(int e=2;e<=3;e++){evaluateCore(p,f.n,e);auto target=e==2?unitAt2:unitAt3;nmod_poly_mul(t,target,p);nmod_poly_swap(t,target);}}nmod_poly_clear(p);nmod_poly_clear(t);}
