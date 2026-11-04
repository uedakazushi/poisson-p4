// SAME fixed K120/t14 three-term numerator; nested integer lattices.
#define main frozen_fixed3_main_NOT_CALLED
#include "s61_nonsquare_prePOWER_h22_pure3_fixed_minor_20261104_v63.cpp"
#undef main
#include <flint/fmpz_vec.h>
#include <memory>
#include <stdexcept>
#include <array>
#include <sstream>
#include <filesystem>
struct Z{fmpz_t z;Z(){fmpz_init(z);}~Z(){fmpz_clear(z);}Z(const Z&)=delete;};
struct Dense;
set<const Dense*>dense_alive;long actual_peak=0;unsigned long long scalar_ops=0;int currentnode=-1,outer_done=0,inner_done=0;long leftpos=-1,rightpos=-1;string currentkind="parse",folder;
struct Dense{
 int lo[3],hi[3],dl,dh,ne,nd;vector<int>index;vector<array<long,3>>exps;fmpz*coeff=nullptr;long count=0;Z den;
 Dense(istream&f){for(int j=0;j<3;j++)f>>lo[j]>>hi[j];f>>dl>>dh;string s;f>>s;assert(!fmpz_set_str(den.z,s.c_str(),10));ne=hi[1]-lo[1]+1;nd=hi[2]-lo[2]+1;index.assign((hi[0]-lo[0]+1)*ne*nd,-1);for(long a=lo[0];a<=hi[0];a++)for(long e=lo[1];e<=hi[1];e++)for(long d=lo[2];d<=hi[2];d++)if(a+e+d>=dl&&a+e+d<=dh){index[((a-lo[0])*ne+e-lo[1])*nd+d-lo[2]]=count++;exps.push_back({a,e,d});}coeff=_fmpz_vec_init(count);dense_alive.insert(this);}
 ~Dense(){dense_alive.erase(this);if(coeff)_fmpz_vec_clear(coeff,count);}
 long slot(long a,long e,long d)const{assert(a>=lo[0]&&a<=hi[0]&&e>=lo[1]&&e<=hi[1]&&d>=lo[2]&&d<=hi[2]);long i=index[((a-lo[0])*ne+e-lo[1])*nd+d-lo[2]];assert(i>=0&&i<count);return i;}
};
unique_ptr<Dense>output_dense,bottom_dense;
void bounded(long scratch=0){long live=15000+scratch;for(auto p:alive)live+=fmpq_mpoly_length(p->p,ctx);for(auto p:dense_alive)live+=p->count;actual_peak=max(actual_peak,live);if(live>150000)throw runtime_error("150000_total_coefficient_slots");if(chrono::duration<double>(chrono::steady_clock::now()-began).count()>50)throw runtime_error("50s_native_guard");}
void rebase(L&f,long target){assert(target<=f.a);long k=f.a-target;if(!k||zero(f)){f.a=zero(f)?0:target;return;}ulong e[3];long old=f.a;for(long i=0;i<fmpq_mpoly_length(f.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(e,f.p,i,ctx);long before=long(e[0])+old;e[0]+=k;fmpq_mpoly_set_term_exp_ui(f.p,i,e,ctx);assert(long(e[0])+target==before);}f.a=target;assert(fmpq_mpoly_is_canonical(f.p,ctx));}
void normalize(L&f){if(zero(f)){f.a=0;return;}ulong e[3];long m=LONG_MAX;for(long i=0;i<fmpq_mpoly_length(f.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(e,f.p,i,ctx);m=min(m,long(e[0]));}if(!m)return;for(long i=0;i<fmpq_mpoly_length(f.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(e,f.p,i,ctx);e[0]-=m;fmpq_mpoly_set_term_exp_ui(f.p,i,e,ctx);}f.a+=m;assert(fmpq_mpoly_is_canonical(f.p,ctx));}
L opvalue(L&x,L&y,char op){L z;if(op=='m'){fmpq_mpoly_mul(z.p,x.p,y.p,ctx);z.a=x.a+y.a;}else if(zero(x)){z=y;if(op=='s')fmpq_mpoly_neg(z.p,z.p,ctx);}else if(zero(y))z=x;else{long a=min(x.a,y.a);rebase(x,a);rebase(y,a);if(op=='a')fmpq_mpoly_add(z.p,x.p,y.p,ctx);else fmpq_mpoly_sub(z.p,x.p,y.p,ctx);z.a=a;}normalize(z);bounded(fmpq_mpoly_length(z.p,ctx));return z;}
void progress(){ofstream f(folder+"/progress.json");f<<"{\"node\":"<<currentnode<<",\"kind\":\""<<currentkind<<"\",\"completed_outer_terms\":"<<outer_done<<",\"completed_inner_terms\":"<<inner_done<<",\"scalar_updates\":"<<scalar_ops<<",\"left_index\":"<<leftpos<<",\"right_index\":"<<rightpos<<",\"peak_slots\":"<<actual_peak<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";}
void jp_dense(ostream&o,const Dense&D){o<<'[';bool first=true;fmpq_t q;fmpq_init(q);for(long i=0;i<D.count;i++)if(!fmpz_is_zero(D.coeff+i)){if(!first)o<<',';first=false;fmpz_set(fmpq_numref(q),D.coeff+i);fmpz_set(fmpq_denref(q),D.den.z);fmpq_canonicalise(q);char*s=fmpq_get_str(nullptr,10,q);auto e=D.exps[i];o<<"[["<<e[0]<<','<<e[1]<<','<<e[2]<<"],\""<<s<<"\"]";flint_free(s);if(o.tellp()>29*1024*1024)throw runtime_error("28MiB_single_dense_output_guard");}fmpq_clear(q);o<<']';}
vector<array<long,3>>exponents(const L&p){vector<array<long,3>>v;ulong e[3];for(long i=0;i<fmpq_mpoly_length(p.p,ctx);i++){fmpq_mpoly_get_term_exp_ui(e,p.p,i,ctx);v.push_back({long(e[0])+p.a,long(e[1]),long(e[2])});}return v;}
void check_tick(){scalar_ops++;if(!(scalar_ops&65535))bounded();}
void actual_scale_record(ostream&o,int k,const fmpz_t scale){char*s=fmpz_get_str(nullptr,10,scale);o<<"{\"event\":"<<k<<",\"actual_integer_content_scale\":\""<<s<<"\"}\n";flint_free(s);}
void fuse_sparse(Dense&D,const L&a,const L&b,int sign,int event,ostream&scales){if(zero(a)||zero(b))return;fmpq_t q;fmpq_init(q);fmpq_mul(q,a.p->content,b.p->content);fmpq_mul_fmpz(q,q,D.den.z);assert(fmpz_is_one(fmpq_denref(q)));Z scale,left;fmpz_set(scale.z,fmpq_numref(q));if(sign<0)fmpz_neg(scale.z,scale.z);actual_scale_record(scales,event,scale.z);auto ae=exponents(a),be=exponents(b);
 for(long i=0;i<long(ae.size());i++){leftpos=i;fmpz_mul(left.z,a.p->zpoly->coeffs+i,scale.z);for(long j=0;j<long(be.size());j++){rightpos=j;long k=D.slot(ae[i][0]+be[j][0],ae[i][1]+be[j][1],ae[i][2]+be[j][2]);fmpz_addmul(D.coeff+k,left.z,b.p->zpoly->coeffs+j);check_tick();}}fmpq_clear(q);bounded();}
void fuse_dense(Dense&D,const L&a,const Dense&b,int sign,int event,ostream&scales){if(zero(a))return;Z ratio,scale,left;assert(fmpz_divisible(D.den.z,b.den.z));fmpz_divexact(ratio.z,D.den.z,b.den.z);fmpq_t q;fmpq_init(q);fmpq_mul_fmpz(q,a.p->content,ratio.z);assert(fmpz_is_one(fmpq_denref(q)));fmpz_set(scale.z,fmpq_numref(q));if(sign<0)fmpz_neg(scale.z,scale.z);actual_scale_record(scales,event,scale.z);auto ae=exponents(a);
 for(long i=0;i<long(ae.size());i++){leftpos=i;fmpz_mul(left.z,a.p->zpoly->coeffs+i,scale.z);for(long j=0;j<b.count;j++)if(!fmpz_is_zero(b.coeff+j)){rightpos=j;auto e=b.exps[j];long k=D.slot(ae[i][0]+e[0],ae[i][1]+e[1],ae[i][2]+e[2]);fmpz_addmul(D.coeff+k,left.z,b.coeff+j);check_tick();}}fmpq_clear(q);bounded();}

#include "s61_nonsquare_prePOWER_K120_t14_nested_entry_main_20261104_v63.hpp"
