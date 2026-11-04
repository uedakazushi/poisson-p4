// Independent symbolic Fp[ a,e,d ] source replay. Ordinary nmod polynomial
// products/sums for top5 (no producer integer lattice), fused outer accumulation.
#include <flint/nmod_mpoly.h>
#include <flint/fmpz.h>
#include <flint/ulong_extras.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;constexpr ulong PRIME=32003; nmod_mpoly_ctx_t ctx;
struct P;set<const P*>live;long peak=0,dense_slots=0;auto began=chrono::steady_clock::now();
struct P{nmod_mpoly_t p;long a=0;P(){nmod_mpoly_init(p,ctx);live.insert(this);}~P(){live.erase(this);nmod_mpoly_clear(p,ctx);}P(const P&)=delete;P&operator=(const P&)=delete;};
void guard(long vectors=0){long n=1000+dense_slots+vectors;for(auto*x:live)n+=nmod_mpoly_length(x->p,ctx);peak=max(peak,n);if(n>150000)throw runtime_error("150k owned coefficient guard");if(chrono::duration<double>(chrono::steady_clock::now()-began).count()>49)throw runtime_error("49s independent guard");}
void copy_shift(P&out,const P&x,long target){assert(target<=x.a);nmod_mpoly_set(out.p,x.p,ctx);out.a=target;if(target!=x.a){ulong ex[3];for(long j=0;j<nmod_mpoly_length(out.p,ctx);j++){nmod_mpoly_get_term_exp_ui(ex,out.p,j,ctx);ex[0]+=x.a-target;nmod_mpoly_set_term_exp_ui(out.p,j,ex,ctx);}assert(nmod_mpoly_is_canonical(out.p,ctx));}}
void sum(P&out,const P&x,const P&y,bool subtract){long a=min(x.a,y.a);P yy;copy_shift(yy,y,a);if(&out!=&x)copy_shift(out,x,a);else if(out.a!=a){P xx;copy_shift(xx,x,a);nmod_mpoly_swap(out.p,xx.p,ctx);out.a=a;}if(subtract)nmod_mpoly_sub(out.p,out.p,yy.p,ctx);else nmod_mpoly_add(out.p,out.p,yy.p,ctx);out.a=a;guard();}
void mul(P&out,const P&x,const P&y){nmod_mpoly_mul(out.p,x.p,y.p,ctx);out.a=x.a+y.a;guard();}
struct Op{char k='x';int x=-1,y=-1;long shift=0,upper=0;};
struct Inner{int x,y,sign;vector<int>need;};struct Event{int w,sign;vector<int>need;vector<Inner>inner;};
struct Term{long a,e,d;ulong c;};
vector<Term>terms(const P&p){vector<Term>v;ulong ex[3];for(long j=0;j<nmod_mpoly_length(p.p,ctx);j++){nmod_mpoly_get_term_exp_ui(ex,p.p,j,ctx);v.push_back({long(ex[0])+p.a,long(ex[1]),long(ex[2]),nmod_mpoly_get_term_coeff_ui(p.p,j,ctx)});}return v;}
ulong modstr(const string&s){fmpz_t z;fmpz_init(z);assert(!fmpz_set_str(z,s.c_str(),10));ulong v=fmpz_fdiv_ui(z,PRIME);fmpz_clear(z);return v;}
int main(int argc,char**argv){assert(argc==4);nmod_mpoly_ctx_init(ctx,3,ORD_DEGREVLEX,PRIME);ofstream progress(argv[3]);try{{
 ifstream in(argv[1]);int n;in>>n;assert(n==11523);vector<Op>ops(n);vector<unique_ptr<P>>atoms(n);long atomcount=0;
 for(int i=0;i<n;i++){auto&o=ops[i];int id;in>>id>>o.k;assert(id==i);if(o.k=='x')continue;in>>o.upper;if(o.k=='p'){long cnt;in>>o.shift>>cnt;atoms[i]=make_unique<P>();atoms[i]->a=o.shift;for(long j=0;j<cnt;j++){ulong ex[3],c;in>>ex[0]>>ex[1]>>ex[2]>>c;assert(c<PRIME);nmod_mpoly_set_coeff_ui_ui(atoms[i]->p,c,ex,ctx);}assert(nmod_mpoly_length(atoms[i]->p,ctx)<=o.upper);atomcount++;}else if(o.k=='k'){in>>o.x>>o.shift;assert(o.x<i);}else{in>>o.x>>o.y;assert(o.x<i&&o.y<i);}}
 long low[3],high[3],dglo,dghi;for(int i=0;i<3;i++)in>>low[i]>>high[i];in>>dglo>>dghi;long ne;in>>ne;assert(ne==21);vector<Event>ev(ne);
 for(auto&e:ev){long nn;in>>e.w>>e.sign>>nn;e.need.resize(nn);for(auto&i:e.need)in>>i;long ni;in>>ni;assert(ni==10);e.inner.resize(ni);for(auto&t:e.inner){in>>t.x>>t.y>>t.sign>>nn;t.need.resize(nn);for(auto&i:t.need)in>>i;}}
 long E=high[1]-low[1]+1,D=high[2]-low[2]+1,A=high[0]-low[0]+1;vector<int>lookup(A*E*D,-1);vector<ulong>acc;for(long a=low[0];a<=high[0];a++)for(long e=low[1];e<=high[1];e++)for(long d=low[2];d<=high[2];d++)if(a+e+d>=dglo&&a+e+d<=dghi){lookup[((a-low[0])*E+e-low[1])*D+d-low[2]]=acc.size();acc.push_back(0);}dense_slots=acc.size();assert(dense_slots==78338);guard();
 auto index=[&](long a,long e,long d){assert(a>=low[0]&&a<=high[0]&&e>=low[1]&&e<=high[1]&&d>=low[2]&&d<=high[2]);int j=lookup[((a-low[0])*E+e-low[1])*D+d-low[2]];assert(j>=0);return j;};long evalcount=0;
 auto evaluate=[&](const vector<int>&need,const vector<int>&roots){vector<int>uses(n);for(int i:need){auto&o=ops[i];if(o.k!='p'){uses[o.x]++;if(o.k!='k')uses[o.y]++;}}for(int i:roots)uses[i]++;vector<unique_ptr<P>>v(n);
  auto ptr=[&](int i)->P*{return atoms[i]?atoms[i].get():v[i].get();};
  for(int i:need){auto&o=ops[i];assert(o.k!='x');if(o.k!='p'){assert(ptr(o.x));v[i]=make_unique<P>();if(o.k=='k'){nmod_mpoly_set(v[i]->p,ptr(o.x)->p,ctx);v[i]->a=ptr(o.x)->a+o.shift;}else{assert(ptr(o.y));if(o.k=='m')mul(*v[i],*ptr(o.x),*ptr(o.y));else sum(*v[i],*ptr(o.x),*ptr(o.y),o.k=='s');}assert(nmod_mpoly_length(v[i]->p,ctx)<=o.upper);guard();assert(uses[o.x]);if(!--uses[o.x])v[o.x].reset();if(o.k!='k'){assert(uses[o.y]);if(!--uses[o.y])v[o.y].reset();}}evalcount++;}
  vector<unique_ptr<P>>out;for(int i:roots){assert(ptr(i));if(atoms[i]){auto z=make_unique<P>();nmod_mpoly_set(z->p,atoms[i]->p,ctx);z->a=atoms[i]->a;out.push_back(move(z));}else out.push_back(move(v[i]));assert(!--uses[i]);}for(int u:uses)assert(!u);guard();return out;};
 long updates=0,innerdone=0,outerdone=0;
 for(auto&e:ev){auto w=evaluate(e.need,{e.w});P top;
  for(auto&t:e.inner){auto r=evaluate(t.need,{t.x,t.y});P product;mul(product,*r[0],*r[1]);if(t.sign<0)nmod_mpoly_neg(product.p,product.p,ctx);sum(top,top,product,false);innerdone++;guard();}
  auto tt=terms(top),ww=terms(*w[0]);guard(tt.size()+ww.size());
  // No ordinary final product coexists with the dense accumulated output.
  for(auto&t:tt){for(auto&u:ww){long j=index(t.a+u.a,t.e+u.e,t.d+u.d);ulong c=(t.c*u.c)%PRIME;acc[j]=(acc[j]+(e.sign>0?c:PRIME-c))%PRIME;updates++;}if((updates&262143)==0)guard(tt.size()+ww.size());}
  outerdone++;progress<<"outer "<<outerdone<<" inner "<<innerdone<<" owned_peak "<<peak<<"\n"<<flush;guard();}
 assert(innerdone==210&&outerdone==21);atoms.clear();ops.clear();ev.clear();guard();
 string ds;long count;in>>ds>>count;ulong gd=modstr(ds);assert(gd);ulong inv=n_invmod(gd,PRIME);long nonzero=0,compared=0;for(ulong x:acc)if(x)nonzero++;for(long i=0;i<count;i++){long a,e,d;string cs;in>>a>>e>>d>>cs;ulong c=(modstr(cs)*inv)%PRIME;assert(acc[index(a,e,d)]==c);if(c)compared++;}assert(nonzero==compared&&count==14986&&nonzero>0);string extra;assert(!(in>>extra));
 long wi=-1;for(long i=0;i<long(acc.size());i++)if(acc[i]){wi=i;break;}long wa=0,we=0,wd=0;for(long a=low[0];a<=high[0];a++)for(long e=low[1];e<=high[1];e++)for(long d=low[2];d<=high[2];d++)if(a+e+d>=dglo&&a+e+d<=dghi&&index(a,e,d)==wi){wa=a;we=e;wd=d;}
 ofstream out(argv[2]);out<<"{\"all_passed\":true,\"prime\":32003,\"source_atom_count\":"<<atomcount<<",\"source_evaluations\":"<<evalcount<<",\"ordinary_nmod_inner_products\":"<<innerdone<<",\"outer_fused_products\":"<<outerdone<<",\"all14986_stored_QQ_coefficients_reduce_to_same_polynomial\":true,\"nonzero_reduced_terms\":"<<nonzero<<",\"nonzero_coefficient\":[["<<wa<<","<<we<<","<<wd<<"],"<<acc[wi]<<"],\"output_global_denominator_modp\":"<<gd<<",\"actual_scalar_updates\":"<<updates<<",\"peak_visible_slots_with1000reserve\":"<<peak<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<",\"full_QQ_coefficient_identity_NOT_claimed\":true}\n";
 }}catch(const exception&e){progress<<"STOP "<<e.what()<<"\n";cerr<<e.what()<<"\n";return 2;}nmod_mpoly_ctx_clear(ctx);flint_cleanup();}
