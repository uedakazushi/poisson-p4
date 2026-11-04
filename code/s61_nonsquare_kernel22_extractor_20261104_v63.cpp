// Exact scalar/linear circuit evaluation on all22 frozen K basis vectors.
#include "s61_binary_hyperplane_cubic_qq3_20261102.hpp"
#include <chrono>
#include <memory>
struct SNode{char t;int i=0,j=0;};
struct LNode{char t;int row=0;vector<pair<int,int>> terms;};
int main(int argc,char**argv){assert(argc==3);auto start=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,3,ORD_LEX);{
 ifstream f(argv[1]);ofstream log(argv[2]);int na,ns,nl,nc,no;f>>na>>ns>>nl>>nc>>no;
 vector<Rat>atoms(na);const char*names[]={"a","e","d"};vector<Rat>units;
 for(const char*t:{"a","d","500*a+90*e^2+3*e-24","50*a*d+18*d*e^2-9*d*e-2*d+60*e-40","100*a+18*e^2-9*e-2","5*a*d-6*e-4"}){Rat p;assert(!fmpq_mpoly_set_str_pretty(p.n,t,names,ctx));units.push_back(p);}
 fmpq_mpoly_t den,tmp;fmpq_mpoly_init(den,ctx);fmpq_mpoly_init(tmp,ctx);
 for(auto&a:atoms){readpoly(f,a.n);readpoly(f,a.d);assert(!fmpq_mpoly_is_zero(a.d,ctx));fmpq_mpoly_set(den,a.d,ctx);for(auto&u:units)while(fmpq_mpoly_divides(tmp,den,u.n,ctx))fmpq_mpoly_swap(tmp,den,ctx);assert(fmpq_mpoly_is_fmpq(den,ctx));}
 fmpq_mpoly_clear(den,ctx);fmpq_mpoly_clear(tmp,ctx);
 vector<SNode>sn(ns);for(auto&n:sn){f>>n.t;if(n.t=='c')f>>n.i;else if(n.t!='z')f>>n.i>>n.j;}
 vector<LNode>ln(nl);vector<int>uses(nl);for(int i=0;i<nl;i++){auto&n=ln[i];f>>n.t;if(n.t=='v')f>>n.row;else if(n.t=='l'){int count;f>>count;while(count--){int s,k;f>>s>>k;assert(s<ns&&k<i);n.terms.push_back({s,k});uses[k]++;}}else assert(n.t=='z');}
 vector<pair<int,int>>checks(nc),outputs(no);for(auto&v:checks)f>>v.first>>v.second;for(auto&v:outputs){f>>v.first>>v.second;uses[v.first]++;}
 int rows,cols;f>>rows>>cols;assert(rows==174&&cols==22);vector<vector<Rat>>K(rows,vector<Rat>(cols));for(auto&r:K)for(auto&v:r){readpoly(f,v.n);readpoly(f,v.d);}string extra;assert(!(f>>extra));
 vector<Rat>sv(ns);for(int i=0;i<ns;i++){auto&n=sn[i];if(n.t=='z')continue;if(n.t=='c')sv[i]=atoms[n.i];else{assert(n.i<i&&n.j<i);if(n.t=='m')sv[i]=mul(sv[n.i],sv[n.j]);else if(n.t=='s')sv[i]=sub(sv[n.i],sv[n.j]);else{assert(n.t=='a');sv[i]=sub(sv[n.i],sub(Rat(),sv[n.j]));}}}
 for(int i=0;i<nc;i++){auto [root,val]=checks[i];Rat expected;fmpq_mpoly_set_si(expected.n,val,ctx);assert(sub(sv[root],expected).zero());log<<"SCALAR_CHECK "<<i<<" EXACT "<<val<<'\n';}log.flush();
 vector<unique_ptr<vector<Rat>>>lv(nl);
 for(int i=0;i<nl;i++){
  auto&n=ln[i];lv[i]=make_unique<vector<Rat>>(22);
  if(n.t=='v')*lv[i]=K[n.row];
  if(n.t=='l')for(auto [s,k]:n.terms){if(!sv[s].zero())for(int j=0;j<22;j++)if(!(*lv[k])[j].zero())(*lv[i])[j]=sub((*lv[i])[j],sub(Rat(),mul(sv[s],(*lv[k])[j])));if(!--uses[k])lv[k].reset();}
  if(!uses[i])lv[i].reset();
 }
 for(int i=0;i<no;i++){auto [root,want]=outputs[i];assert(lv[root]);for(int j=0;j<22;j++){Rat expected;if(j==want)fmpq_mpoly_one(expected.n,ctx);assert(sub((*lv[root])[j],expected).zero());}log<<"LINEAR_OUTPUT "<<i<<" EXACT_UNIT_ROW "<<want<<" ALL22\n";if(!--uses[root])lv[root].reset();}
 log<<"ALL_FULLU_ATOM_DENOMINATORS_KNOWN_UNITS\nEXACT_SCALAR_CHECKS "<<nc<<" LINEAR_OUTPUTS "<<no<<" K_COLUMNS22 PASS SEC "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';log.flush();
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();}
