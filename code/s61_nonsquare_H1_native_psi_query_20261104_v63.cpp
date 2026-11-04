// ONE factored-unit coefficient query: psi row113 ONLY. No Gaussian or generic GCD.
#include <flint/fmpq.h>
#include <flint/fmpq_mpoly.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>
using namespace std;fmpq_mpoly_ctx_t ctx;using E=array<long,6>;
struct P{fmpq_mpoly_t p;P(){fmpq_mpoly_init(p,ctx);}P(const P&v):P(){fmpq_mpoly_set(p,v.p,ctx);}P(P&&v):P(){fmpq_mpoly_swap(p,v.p,ctx);}P&operator=(P v){fmpq_mpoly_swap(p,v.p,ctx);return*this;}~P(){fmpq_mpoly_clear(p,ctx);}};
struct V{P n;E den{};};struct Atom{P n;E den{};string c;};struct Node{long id,x=-1,y=-1;char op;};
vector<P>known;map<pair<int,long>,P>powers;chrono::steady_clock::time_point start;
long maxinter=0,maxoutput=0,divisions=0,normchecks=0,nodechecks=0,evalnodes=0,maxliveterms=0,maxlivenodes=0,oldcross=0;
long pendingid=-1;string pendingstage="initial";ofstream trace;
double seconds(){return chrono::duration<double>(chrono::steady_clock::now()-start).count();}
void timer(){if(seconds()>50)throw runtime_error("50s_native_arithmetic_boundary");}
long terms(const P&p){return fmpq_mpoly_length(p.p,ctx);}
void intermediate(const P&p){long n=terms(p);maxinter=max(maxinter,n);if(n>150000)throw runtime_error("150000_intermediate_terms");timer();}
void output(const P&p){long n=terms(p);maxoutput=max(maxoutput,n);if(n>100000)throw runtime_error("100000_output_terms");timer();}
long double totalbound(long d){if(d<0)return 0;return(long double)(d+1)*(d+2)*(d+3)/6;}
P mul(const P&a,const P&b){
 P z;if(!terms(a)||!terms(b))return z;
 long d=fmpq_mpoly_total_degree_si(a.p,ctx)+fmpq_mpoly_total_degree_si(b.p,ctx);
 long double bound=min((long double)terms(a)*terms(b),totalbound(d));if(bound>150000)throw runtime_error("150000_preproduct_support_bound");
 fmpq_mpoly_mul(z.p,a.p,b.p,ctx);intermediate(z);return z;
}
P addpoly(const P&a,const P&b,int sign){
 if(terms(a)+terms(b)>150000)throw runtime_error("150000_preaddition_support_bound");
 P z;if(sign>0)fmpq_mpoly_add(z.p,a.p,b.p,ctx);else fmpq_mpoly_sub(z.p,a.p,b.p,ctx);intermediate(z);return z;
}
const P&power(int k,long n){
 assert(n>=0);auto key=make_pair(k,n);auto it=powers.find(key);if(it!=powers.end())return it->second;
 long d=fmpq_mpoly_total_degree_si(known[k].p,ctx);if(totalbound(d*n)>150000)throw runtime_error("150000_named_power_support_bound");
 P z;assert(fmpq_mpoly_pow_ui(z.p,known[k].p,n,ctx));intermediate(z);auto ans=powers.emplace(key,move(z));return ans.first->second;
}
P align(P n,const E&v){
 for(int k=0;k<6;k++){assert(v[k]>=0);if(v[k])n=mul(n,power(k,v[k]));}return n;
}
void normalize(V&v){
 if(!terms(v.n)){v.den.fill(0);return;}
 P before=v.n;E old=v.den;pendingstage="named_factor_strip";
 for(int k=0;k<6;k++){
  if(fmpq_mpoly_is_fmpq(v.n.p,ctx))break;
  P q;
  while(fmpq_mpoly_divides(q.p,v.n.p,known[k].p,ctx)){
   P check=mul(q,known[k]);assert(fmpq_mpoly_equal(check.p,v.n.p,ctx));v.n=q;v.den[k]--;divisions++;
   if(abs(v.den[k])>1000)throw runtime_error("1000_signed_unit_exponent_guard");timer();
  }
 }
 E stripped{};for(int k=0;k<6;k++){stripped[k]=old[k]-v.den[k];assert(stripped[k]>=0);}
 P rebuilt=align(v.n,stripped);assert(fmpq_mpoly_equal(rebuilt.p,before.p,ctx));normchecks++;output(v.n);
}
V product(const V&a,const V&b){
 V z;if(!terms(a.n)||!terms(b.n))return z;
 for(int k=0;k<6;k++){z.den[k]=a.den[k]+b.den[k];if(abs(z.den[k])>1000)throw runtime_error("1000_signed_unit_exponent_guard");}
 pendingstage="polynomial_product";z.n=mul(a.n,b.n);normalize(z);return z;
}
V sum(const V&a,const V&b,int sign){
 if(!terms(b.n))return a;
 if(!terms(a.n)){V z=b;if(sign<0)fmpq_mpoly_neg(z.n.p,z.n.p,ctx);return z;}
 V z;E da{},db{};for(int k=0;k<6;k++){z.den[k]=max(a.den[k],b.den[k]);da[k]=z.den[k]-a.den[k];db[k]=z.den[k]-b.den[k];}
 pendingstage="lcm_unit_alignment";P x=align(a.n,da),y=align(b.n,db);z.n=addpoly(x,y,sign);normalize(z);return z;
}
void readp(istream&f,P&q){
 long n;f>>n;assert(f&&n>=0);fmpq_t c;fmpq_init(c);ulong ex[3];
 for(long i=0;i<n;i++){string z;f>>z>>ex[0]>>ex[1]>>ex[2];assert(f&&!fmpq_set_str(c,z.c_str(),10));fmpq_canonicalise(c);assert(!fmpq_is_zero(c));fmpq_mpoly_push_term_fmpq_ui(q.p,c,ex,ctx);}
 fmpq_mpoly_sort_terms(q.p,ctx);fmpq_mpoly_combine_like_terms(q.p,ctx);assert(terms(q)==n);fmpq_clear(c);
}
V sourcefraction(const P&n,const P&d){
 P r=d,q;V z;z.n=n;fmpq_t c;fmpq_init(c);
 for(int k=0;k<6;k++)while(fmpq_mpoly_divides(q.p,r.p,known[k].p,ctx)){r=q;z.den[k]++;}
 assert(fmpq_mpoly_is_fmpq(r.p,ctx));fmpq_mpoly_get_fmpq(c,r.p,ctx);assert(!fmpq_is_zero(c));fmpq_mpoly_scalar_div_fmpq(z.n.p,z.n.p,c,ctx);normalize(z);
 // Literal cross-product with the stored fraction, allowing signed unit exponents.
 E plus{},minus{};for(int k=0;k<6;k++){plus[k]=max(z.den[k],0L);minus[k]=max(-z.den[k],0L);}
 P left=align(mul(z.n,d),minus),right=align(n,plus);assert(fmpq_mpoly_equal(left.p,right.p,ctx));oldcross++;fmpq_clear(c);return z;
}
long byteswritten=0;
void writev(ostream&f,const V&v){
 for(long d:v.den)f<<d<<' ';f<<terms(v.n);fmpq_t c;fmpq_init(c);ulong ex[3];
 for(long i=0;i<terms(v.n);i++){
  fmpq_mpoly_get_term_coeff_fmpq(c,v.n.p,i,ctx);fmpq_mpoly_get_term_exp_ui(ex,v.n.p,i,ctx);char*z=fmpq_get_str(nullptr,10,c);byteswritten+=strlen(z)+100;if(byteswritten>30000000){flint_free(z);throw runtime_error("30MB_output_artifact_guard");}
  f<<' '<<z;flint_free(z);for(ulong x:ex)f<<' '<<x;if(i%100==0)timer();
 }f<<'\n';fmpq_clear(c);
}
int main(int argc,char**argv){
 assert(argc==3);start=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);string dir=argv[2];trace.open(dir+"/native-progress.jsonl");bool passed=false;string stop="";long expected=0,rawnonzero=0,reducednonzero=0;
 try{
  ifstream f(argv[1]);long nk,na,nn,nr,no,nmax;f>>nk>>na>>nn>>nr>>no>>nmax;assert(nk==6&&na==631&&nn==31860&&nr==23&&no==7);expected=nn;
  for(int k=0;k<6;k++){known.emplace_back();readp(f,known.back());}
  map<long,Atom>atoms;
  for(int i=0;i<na;i++){long id;f>>id;Atom z;f>>z.c;for(long&e:z.den)f>>e;readp(f,z.n);atoms.emplace(id,move(z));}
  vector<Node>nodes;vector<long>refs(nmax,0);for(int i=0;i<nn;i++){Node n;f>>n.id>>n.op;assert(0<=n.id&&n.id<nmax);if(n.op=='c')f>>n.x;else if(n.op!='z'){f>>n.x>>n.y;assert(0<=n.x&&n.x<n.id&&0<=n.y&&n.y<n.id);refs[n.x]++;refs[n.y]++;}nodes.push_back(n);}
  vector<long>roots(nr);for(long&r:roots){f>>r;refs[r]++;}
  vector<long>pivots(no);vector<vector<V>>old(no);
  for(int i=0;i<no;i++){f>>pivots[i];for(int j=0;j<nr;j++){P n,d;readp(f,n);readp(f,d);old[i].push_back(sourcefraction(n,d));}}
  string extra;assert(!(f>>extra));vector<unique_ptr<V>>values(nmax);long livet=0,liven=0;
  for(const Node&n:nodes){
   pendingid=n.id;pendingstage="node";timer();auto v=make_unique<V>();
   if(n.op=='z'){}
   else if(n.op=='c'){
    const Atom&t=atoms.at(n.x);v->n=t.n;v->den=t.den;fmpq_t c;fmpq_init(c);assert(!fmpq_set_str(c,t.c.c_str(),10));fmpq_canonicalise(c);fmpq_mpoly_scalar_div_fmpq(v->n.p,v->n.p,c,ctx);fmpq_clear(c);normalize(*v);
   }else{
    assert(values[n.x]&&values[n.y]);if(n.op=='m')*v=product(*values[n.x],*values[n.y]);else{assert(n.op=='a'||n.op=='s');*v=sum(*values[n.x],*values[n.y],n.op=='a'?1:-1);}nodechecks++;
   }
   livet+=terms(v->n);liven++;values[n.id]=move(v);evalnodes++;maxliveterms=max(maxliveterms,livet);maxlivenodes=max(maxlivenodes,liven);
   if(n.op!='z'&&n.op!='c')for(long dep:{n.x,n.y}){refs[dep]--;if(!refs[dep]){livet-=terms(values[dep]->n);liven--;values[dep].reset();}}
   if(evalnodes%1000==0)trace<<"{\"node\":"<<n.id<<",\"evaluated\":"<<evalnodes<<",\"live_terms\":"<<livet<<",\"live_nodes\":"<<liven<<",\"seconds\":"<<seconds()<<"}\n"<<flush;
  }
  vector<V>raw,row,cofs;for(long root:roots){assert(values[root]);raw.push_back(*values[root]);}row=raw;
  for(int i=0;i<no;i++){
   V c=row[pivots[i]];cofs.push_back(c);
   for(int j=0;j<nr;j++)row[j]=sum(row[j],product(c,old[i][j]),-1);
   assert(!terms(row[pivots[i]].n));
  }
  // All23 final full-source Schur identities are additionally checked by reconstruction.
  for(int j=0;j<nr;j++){
   V z=row[j];for(int i=0;i<no;i++)z=sum(z,product(cofs[i],old[i][j]),1);
   V difference=sum(z,raw[j],-1);assert(!terms(difference.n));
  }
  ofstream rf(dir+"/raw23-factored.txt"),sf(dir+"/schur23-factored.txt"),cf(dir+"/schur7-coefficients-factored.txt");rf<<23<<'\n';sf<<23<<'\n';cf<<7<<'\n';
  for(const V&v:raw){writev(rf,v);rawnonzero+=terms(v.n)>0;}for(const V&v:row){writev(sf,v);reducednonzero+=terms(v.n)>0;}for(const V&v:cofs)writev(cf,v);
  passed=true;known.clear();powers.clear();
 }catch(const exception&e){stop=e.what();}
 struct rusage ru;getrusage(RUSAGE_SELF,&ru);ofstream a(dir+"/native-audit.json");
 a<<"{\"all_passed\":"<<(passed?"true":"false")<<",\"stop_reason\":\""<<stop<<"\",\"expected_DAG_nodes\":"<<expected<<",\"evaluated_DAG_nodes\":"<<evalnodes<<",\"arithmetic_node_identities\":"<<nodechecks<<",\"normalization_crossproduct_checks\":"<<normchecks<<",\"exact_named_factor_divisions\":"<<divisions<<",\"old7_input_fraction_crossproducts\":"<<oldcross<<",\"raw23_nonzero\":"<<rawnonzero<<",\"schur23_nonzero\":"<<reducednonzero<<",\"all23_final_schur_identities_checked\":"<<(passed?"true":"false")<<",\"max_intermediate_polynomial_terms\":"<<maxinter<<",\"max_output_polynomial_terms\":"<<maxoutput<<",\"max_live_scalar_polynomial_terms\":"<<maxliveterms<<",\"max_live_scalar_nodes\":"<<maxlivenodes<<",\"peak_RSS_bytes\":"<<ru.ru_maxrss*1024L<<",\"pending_node\":"<<pendingid<<",\"pending_stage\":\""<<pendingstage<<"\",\"native_seconds\":"<<seconds()<<"}\n";
 known.clear();powers.clear();fmpq_mpoly_ctx_clear(ctx);flint_cleanup();return 0;
}
