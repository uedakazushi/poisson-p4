// Bounded direct summand lift test on a pinned INTEGER polynomial matrix.
// A degree<=1 cofactor q has basis(1,a,e), and a direct row-target lift would
// satisfy Mq=unit*e_r.  Evaluation at fixed integer pairs gives NECESSARY
// integer linear equations. Full column rank modulo p of all rows except r
// proves q=0 over QQ; lower modular ranks are only diagnostics.
// CLI input output-json. Matrix format has the two standard metadata lines.
#include <flint/fmpq.h>
#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
int main(int argc,char**argv){assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);int n,m,k;in>>n>>m>>k;int discard;for(int i=0;i<k+n;i++)in>>discard;const ulong p=32003;nmod_t mod;nmod_init(&mod,p);
 vector<ulong>as={1,2,3,5,7,11,13,17,19,23},es={2,3,5,7,11,13,17,19,23,29};int samples=as.size(),cols=3*m;nmod_mat_t N;nmod_mat_init(N,samples*n,cols,p);fmpq_t q;fmpq_init(q);
 for(int i=0;i<n;i++)for(int j=0;j<m;j++){
  int nt;in>>nt;assert(in&&nt>=0);vector<ulong>v(samples);
  for(int t=0;t<nt;t++){string c;ulong a,e,d;in>>c>>a>>e>>d;assert(in&&d==0&&!fmpq_set_str(q,c.c_str(),10));fmpq_canonicalise(q);assert(fmpz_is_one(fmpq_denref(q)));ulong u=fmpz_fdiv_ui(fmpq_numref(q),p);
   for(int s=0;s<samples;s++){ulong term=nmod_mul(u,nmod_mul(nmod_pow_ui(as[s],a,mod),nmod_pow_ui(es[s],e,mod),mod),mod);v[s]=nmod_add(v[s],term,mod);}}
  for(int s=0;s<samples;s++){int r=s*n+i;nmod_mat_entry(N,r,j)=v[s];nmod_mat_entry(N,r,m+j)=nmod_mul(as[s],v[s],mod);nmod_mat_entry(N,r,2*m+j)=nmod_mul(es[s],v[s],mod);}
 }
 nmod_mat_t E;nmod_mat_init(E,samples*n,cols,p);nmod_mat_set(E,N);int rank=nmod_mat_rref(E);cerr<<"EVALUATION_FULL_RANK "<<rank<<" COLUMNS "<<cols<<endl;nmod_mat_clear(E);vector<int>ranks;
 for(int target=0;target<n;target++){nmod_mat_init(E,samples*(n-1),cols,p);int at=0;for(int s=0;s<samples;s++)for(int i=0;i<n;i++)if(i!=target){for(int j=0;j<cols;j++)nmod_mat_entry(E,at,j)=nmod_mat_entry(N,s*n+i,j);at++;}int r=nmod_mat_rref(E);ranks.push_back(r);nmod_mat_clear(E);cerr<<"EXCLUDED_TARGET_ROW "<<target<<" RANK "<<r<<endl;}
 out<<"{\"prime\":"<<p<<",\"source_shape\":["<<n<<','<<m<<"],\"cofactor_basis\":[\"1\",\"a\",\"e\"],\"unknowns\":"<<cols<<",\"evaluation_pairs\":[";
 for(int s=0;s<samples;s++){if(s)out<<',';out<<'['<<as[s]<<','<<es[s]<<']';}out<<"],\"full_rank\":"<<rank<<",\"excluded_row_ranks\":[";for(int i=0;i<n;i++){if(i)out<<',';out<<ranks[i];}out<<"],\"evaluation_matrix\":[";
 for(int i=0;i<samples*n;i++){if(i)out<<',';out<<'[';for(int j=0;j<cols;j++){if(j)out<<',';out<<nmod_mat_entry(N,i,j);}out<<']';}out<<"]}\n";
 fmpq_clear(q);nmod_mat_clear(N);flint_cleanup();}
