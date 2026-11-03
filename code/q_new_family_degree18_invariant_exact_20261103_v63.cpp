// Exact QQ theta/H traces of all60 length9 necklaces and a selected9minor.
#include <flint/fmpq_mat.h>
#include <flint/fmpq.h>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
using namespace std;
struct M{
 fmpq_mat_t a;M(long n=16,long m=16){fmpq_mat_init(a,n,m);}
 M(const M& b):M(fmpq_mat_nrows(b.a),fmpq_mat_ncols(b.a)){fmpq_mat_set(a,b.a);}
 M(M&& b):M(){fmpq_mat_swap(a,b.a);}
 M& operator=(M b){fmpq_mat_swap(a,b.a);return *this;}
 ~M(){fmpq_mat_clear(a);}
};
M mul(const M& a,const M& b){M c(fmpq_mat_nrows(a.a),fmpq_mat_ncols(b.a));fmpq_mat_mul(c.a,a.a,b.a);return c;}
M trans(const M& a){M c(fmpq_mat_ncols(a.a),fmpq_mat_nrows(a.a));fmpq_mat_transpose(c.a,a.a);return c;}
M read(ifstream& f){long n,m;f>>n>>m;assert(f&&n>0&&m>0);M a(n,m);string c;for(long i=0;i<n;i++)for(long j=0;j<m;j++){f>>c;assert(f&&fmpq_set_str(fmpq_mat_entry(a.a,i,j),c.c_str(),10)==0);fmpq_canonicalise(fmpq_mat_entry(a.a,i,j));}return a;}
string str(const fmpq_t a){char* p=fmpq_get_str(nullptr,10,a);string r(p);flint_free(p);return r;}
int main(int argc,char**argv){
 assert(argc==3);{
 ifstream f(argv[1]);ofstream out(argv[2]);int count,width;f>>count>>width;assert(count==9&&width==60);
 M O=read(f),Oi=read(f),Gi=read(f);vector<M> rep;for(int i=0;i<10;i++)rep.push_back(read(f));
 vector<string> words(width);for(auto& w:words){f>>w;assert(w.size()==9&&w.find_first_not_of("01")==string::npos);}
 vector<int> selected(count);for(auto& i:selected){f>>i;assert(0<=i&&i<width);}
 M minor(count,count);fmpq_t half,tmp,sum;fmpq_init(half);fmpq_set_si(half,1,2);fmpq_init(tmp);fmpq_init(sum);
 out<<"{\"all_passed\":true,\"evaluation\":[";
 for(int sample=0;sample<count;sample++){
  M b=read(f);assert(fmpq_mat_nrows(b.a)==16&&fmpq_mat_ncols(b.a)==1);M V(16,10);
  for(int k=0;k<10;k++){M v=mul(mul(O,rep[k]),b);for(int i=0;i<16;i++)fmpq_set(fmpq_mat_entry(V.a,i,k),fmpq_mat_entry(v.a,i,0));}
  M mu=trans(mul(trans(b),V));for(int k=0;k<10;k++)fmpq_mul(fmpq_mat_entry(mu.a,k,0),fmpq_mat_entry(mu.a,k,0),half);
  M co=mul(Gi,mu),theta;
  for(int i=0;i<16;i++)for(int j=0;j<16;j++)for(int k=0;k<10;k++){
   fmpq_mul(tmp,fmpq_mat_entry(co.a,k,0),fmpq_mat_entry(rep[k].a,i,j));fmpq_add(fmpq_mat_entry(theta.a,i,j),fmpq_mat_entry(theta.a,i,j),tmp);
  }
  M H=mul(mul(mul(Oi,V),Gi),trans(V));fmpq_mat_add(H.a,H.a,theta.a);
  map<string,M> table;M I;fmpq_mat_one(I.a);table.emplace("",I);
  for(int n=1;n<=5;n++)for(int v=0;v<(1<<n);v++){
   string w;for(int k=n-1;k>=0;k--)w+=char('0'+((v>>k)&1));
   table.emplace(w,mul(table.at(w.substr(0,n-1)),w.back()=='0'?theta:H));
  }
  M values(1,width);
  for(int k=0;k<width;k++){
   const M& a=table.at(words[k].substr(0,4));const M& b=table.at(words[k].substr(4));auto q=fmpq_mat_entry(values.a,0,k);
   for(int i=0;i<16;i++)for(int j=0;j<16;j++){
    fmpq_mul(tmp,fmpq_mat_entry(a.a,i,j),fmpq_mat_entry(b.a,j,i));fmpq_add(q,q,tmp);
   }
  }
  if(sample)out<<',';out<<'[';
  for(int k=0;k<width;k++){if(k)out<<',';out<<'"'<<str(fmpq_mat_entry(values.a,0,k))<<'"';}out<<']';
  for(int k=0;k<count;k++)fmpq_set(fmpq_mat_entry(minor.a,sample,k),fmpq_mat_entry(values.a,0,selected[k]));
 }
 string extra;assert(!(f>>extra));fmpq_mat_det(sum,minor.a);assert(!fmpq_is_zero(sum));
 out<<"],\"exact_selected9_minor\":\""<<str(sum)<<"\"}\n";
 cerr<<"EXACT_ALL540_QQ_TRACE9_VALUES_AND_9_MINOR_PASS\n";
 fmpq_clear(half);fmpq_clear(tmp);fmpq_clear(sum);
 }flint_cleanup();
}
