#define main frozen_fixed3_main_NOT_CALLED
#include "s61_nonsquare_prePOWER_h22_pure3_fixed_minor_20261104_v63.cpp"
#undef main
#include <filesystem>
#include <memory>
#include <stdexcept>
long precise_peak=0;int lastnode=-1;string folder;
void bounded(long scratch=0){long s=reserve_terms+15000+scratch;for(auto p:alive)s+=fmpq_mpoly_length(p->p,ctx);precise_peak=max(precise_peak,s);if(s>150000)throw runtime_error("150000_total_coefficient_slots");if(chrono::duration<double>(chrono::steady_clock::now()-began).count()>50)throw runtime_error("50s_native_guard");}
void rebase(L&f,long target){assert(target<=f.a);long k=f.a-target;if(!k||zero(f)){f.a=zero(f)?0:target;return;}ulong e[3];long old=f.a;for(long i=0;i<fmpq_mpoly_length(f.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(e,f.p,i,ctx);long actual=long(e[0])+old;e[0]+=k;fmpq_mpoly_set_term_exp_ui(f.p,i,e,ctx);assert(long(e[0])+target==actual);}f.a=target;assert(fmpq_mpoly_is_canonical(f.p,ctx));}
void normal(L&f){if(zero(f)){f.a=0;return;}ulong e[3];long m=LONG_MAX;for(long i=0;i<fmpq_mpoly_length(f.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(e,f.p,i,ctx);m=min(m,long(e[0]));}if(!m)return;for(long i=0;i<fmpq_mpoly_length(f.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(e,f.p,i,ctx);e[0]-=m;fmpq_mpoly_set_term_exp_ui(f.p,i,e,ctx);}f.a+=m;assert(fmpq_mpoly_is_canonical(f.p,ctx));}
L operate(L&x,L&y,char op){L z;if(op=='m'){fmpq_mpoly_mul(z.p,x.p,y.p,ctx);z.a=x.a+y.a;}else if(zero(x)){z=y;if(op=='s')fmpq_mpoly_neg(z.p,z.p,ctx);}else if(zero(y))z=x;else{long s=min(x.a,y.a);rebase(x,s);rebase(y,s);if(op=='a')fmpq_mpoly_add(z.p,x.p,y.p,ctx);else fmpq_mpoly_sub(z.p,x.p,y.p,ctx);z.a=s;}normal(z);bounded(fmpq_mpoly_length(z.p,ctx));return z;}
void files_guard(){long long s=0;for(auto&p:filesystem::directory_iterator(folder))if(p.is_regular_file()&&p.path().filename()!="native-worker")s+=p.file_size();if(s>31*1024*1024)throw runtime_error("31MiB_artifact_guard");}
void savepoly(int id,const L&p){ofstream j(folder+"/root"+to_string(id)+".json");j<<"{\"root\":"<<id<<",\"polynomial\":";jp(j,p);j<<"}\n";j.close();ofstream f(folder+"/root"+to_string(id)+".poly");f<<p.a<<'\n'<<fmpq_mpoly_length(p.p,ctx)<<'\n';fmpq_t c;fmpq_init(c);ulong e[3];for(long i=0;i<fmpq_mpoly_length(p.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(e,p.p,i,ctx);fmpq_mpoly_get_term_coeff_fmpq(c,p.p,i,ctx);char*s=fmpq_get_str(nullptr,10,c);f<<e[0]<<' '<<e[1]<<' '<<e[2]<<' '<<s<<'\n';flint_free(s);}fmpq_clear(c);f.close();files_guard();}
unique_ptr<L>loadpoly(int id){ifstream f(folder+"/root"+to_string(id)+".poly");long s;f>>s;auto p=make_unique<L>(readpoly(f));p->a+=s;string extra;assert(!(f>>extra));return p;}

#include "s61_nonsquare_prePOWER_rawcross18_entry_main_20261104_v63.hpp"
