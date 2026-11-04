// Independent full67 coefficient extraction, 10 actual generators, reversed variables.
#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <flint/nmod_mpoly.h>
using namespace std;
nmod_mpoly_ctx_t ctx;
long maxterms=0;
struct P {
 nmod_mpoly_t z;
 P(){nmod_mpoly_init(z,ctx);}
 P(const P&a):P(){nmod_mpoly_set(z,a.z,ctx);}
 P(P&&a):P(){nmod_mpoly_swap(z,a.z,ctx);}
 P&operator=(P a){nmod_mpoly_swap(z,a.z,ctx);return *this;}
 ~P(){nmod_mpoly_clear(z,ctx);}
};
void guard(const P&a){long n=nmod_mpoly_length(a.z,ctx);maxterms=max(maxterms,n);if(n>100000)throw runtime_error("100000 term cap");}
P product(const P&a,const P&b){P z;nmod_mpoly_mul(z.z,a.z,b.z,ctx);guard(z);return z;}
P unit(){P z;nmod_mpoly_one(z.z,ctx);return z;}
int main(int argc,char**argv){
 assert(argc==3);ifstream f(argv[1]);ulong prime;int nv,ng,nc,nr;
 f>>prime>>nv>>ng>>nc>>nr;assert(f&&prime==32003&&nv==4&&ng==10&&nc==67&&nr==67);
 // Author order was (a,Y,W,B), DEGLEX. This engine uses (B,W,Y,a), LEX,
 // with modular coefficients throughout, and never constructs integer columns.
 nmod_mpoly_ctx_init(ctx,4,ORD_LEX,prime);
 {
  vector<P> generators(10);vector<int> bpower(10);
  for(int j=0;j<10;j++){
   long n;f>>bpower[j]>>n;assert(n>=0);ulong mon[4];
   for(long k=0;k<n;k++){ulong c;f>>c;for(auto&e:mon)f>>e;assert(f&&c<prime);nmod_mpoly_push_term_ui_ui(generators[j].z,c,mon,ctx);}
   nmod_mpoly_sort_terms(generators[j].z,ctx);nmod_mpoly_combine_like_terms(generators[j].z,ctx);guard(generators[j]);
  }
  vector<array<int,10>> columns(67);for(auto&e:columns)for(int&k:e)f>>k;
  vector<array<ulong,4>> rows(67);for(auto&e:rows)for(auto&k:e)f>>k;
  assert(f);string extra;assert(!(f>>extra));
  vector<vector<P>> powers(10);
  for(int j=0;j<10;j++){
   int m=0;for(auto&e:columns)m=max(m,e[j]);powers[j].push_back(unit());
   for(int k=1;k<=m;k++)powers[j].push_back(product(powers[j].back(),generators[j]));
  }
  vector<vector<ulong>> matrix(67,vector<ulong>(67));
  for(int j=0;j<67;j++){
   vector<P> factors;int bs=0;
   for(int i=0;i<10;i++){assert(columns[j][i]>=0);bs+=bpower[i]*columns[j][i];if(columns[j][i])factors.push_back(powers[i][columns[j][i]]);}
   assert(bs<=12);if(bs<12){P B,powB;nmod_mpoly_gen(B.z,0,ctx);assert(nmod_mpoly_pow_ui(powB.z,B.z,12-bs,ctx));factors.push_back(powB);}
   // Multiplication starts with the sparsest factors and folds from the right.
   sort(factors.begin(),factors.end(),[](const P&a,const P&b){return nmod_mpoly_length(a.z,ctx)<nmod_mpoly_length(b.z,ctx);});
   P result=unit();for(auto it=factors.rbegin();it!=factors.rend();++it)result=product(*it,result);
   for(int i=0;i<67;i++)matrix[i][j]=nmod_mpoly_get_coeff_ui_ui(result.z,rows[i].data(),ctx);
  }
  ofstream out(argv[2]);assert(out);out<<"{\"all_passed\":true,\"prime\":32003,\"variable_order\":[\"B\",\"W\",\"Y\",\"a\"],\"term_order\":\"LEX\",\"max_terms\":"<<maxterms<<",\"matrix67\":[";
  for(int i=0;i<67;i++){if(i)out<<',';out<<'[';for(int j=0;j<67;j++){if(j)out<<',';out<<matrix[i][j];}out<<']';}out<<"]}\n";
 }
 nmod_mpoly_ctx_clear(ctx);flint_cleanup();
}
