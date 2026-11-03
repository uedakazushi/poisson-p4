// Exact rational moment/Hessian trace invariants and a selected20x20 minor.
#include <flint/fmpq_mat.h>
#include <flint/fmpq.h>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
using namespace std;
struct M {
 fmpq_mat_t a;
 M(long n=16,long m=16){fmpq_mat_init(a,n,m);}
 M(const M& b):M(fmpq_mat_nrows(b.a),fmpq_mat_ncols(b.a)){fmpq_mat_set(a,b.a);}
 M(M&& b):M(){fmpq_mat_swap(a,b.a);}
 M& operator=(M b){fmpq_mat_swap(a,b.a);return *this;}
 ~M(){fmpq_mat_clear(a);}
};
M mul(const M& a,const M& b){M c(fmpq_mat_nrows(a.a),fmpq_mat_ncols(b.a));fmpq_mat_mul(c.a,a.a,b.a);return c;}
M transpose(const M& a){M c(fmpq_mat_ncols(a.a),fmpq_mat_nrows(a.a));fmpq_mat_transpose(c.a,a.a);return c;}
M read(ifstream& f){long n,m;f>>n>>m;M a(n,m);string word;for(long i=0;i<n;i++)for(long j=0;j<m;j++){f>>word;assert(!fmpq_set_str(fmpq_mat_entry(a.a,i,j),word.c_str(),10));fmpq_canonicalise(fmpq_mat_entry(a.a,i,j));}assert(f.good());return a;}
string text(const fmpq_t x){char* p=fmpq_get_str(nullptr,10,x);string r(p);flint_free(p);return r;}
struct Def{char type;string word;int i=0,j=0;};
int main(int argc,char**argv){
 assert(argc==3);ifstream f(argv[1]);ofstream out(argv[2]);int count;f>>count;assert(count==20);
 M O=read(f),Oi=read(f),Gi=read(f);vector<M> rep;for(int i=0;i<10;i++)rep.push_back(read(f));
 vector<int> basis8(5);for(int& x:basis8)f>>x;
 vector<string>w4;for(int i=0;i<16;i++){string w;for(int j=3;j>=0;j--)w+=char('0'+((i>>j)&1));bool least=true;for(int r=1;r<4;r++)if(w.substr(r)+w.substr(0,r)<w)least=false;if(least)w4.push_back(w);}assert(w4.size()==6);
 int width;f>>width;vector<Def> defs(width);for(auto& d:defs){f>>d.type;if(d.type=='P')f>>d.i>>d.j;else f>>d.word;}
 vector<int> selected(count);for(int& x:selected)f>>x;
 M minor(count,count);fmpq_t half,term,sum,four,tmp;fmpq_init(half);fmpq_set_si(half,1,2);fmpq_init(term);fmpq_init(sum);fmpq_init(four);fmpq_init(tmp);
 out<<"{\"evaluation\":[";
 for(int sample=0;sample<count;sample++){
  M b=read(f);assert(fmpq_mat_nrows(b.a)==16&&fmpq_mat_ncols(b.a)==1);M V(16,10);
  for(int k=0;k<10;k++){M v=mul(mul(O,rep[k]),b);for(int i=0;i<16;i++)fmpq_set(fmpq_mat_entry(V.a,i,k),fmpq_mat_entry(v.a,i,0));}
  M mu=transpose(mul(transpose(b),V));for(int k=0;k<10;k++)fmpq_mul(fmpq_mat_entry(mu.a,k,0),fmpq_mat_entry(mu.a,k,0),half);
  M co=mul(Gi,mu),theta;
  for(int i=0;i<16;i++)for(int j=0;j<16;j++)for(int k=0;k<10;k++){fmpq_mul(term,fmpq_mat_entry(co.a,k,0),fmpq_mat_entry(rep[k].a,i,j));fmpq_add(fmpq_mat_entry(theta.a,i,j),fmpq_mat_entry(theta.a,i,j),term);}
  M H=mul(mul(mul(Oi,V),Gi),transpose(V));fmpq_mat_add(H.a,H.a,theta.a);
  map<string,M> table;M identity;fmpq_mat_one(identity.a);table.emplace("",identity);
  for(int n=1;n<=4;n++)for(int v=0;v<(1<<n);v++){string w;for(int j=n-1;j>=0;j--)w+=char('0'+((v>>j)&1));table.emplace(w,mul(table.at(w.substr(0,n-1)),w.back()=='0'?theta:H));}
  auto tr=[&](const string& w,fmpq_t target){const M& left=table.at(w.substr(0,w.size()/2));const M& right=table.at(w.substr(w.size()/2));fmpq_zero(target);for(int i=0;i<16;i++)for(int j=0;j<16;j++){fmpq_mul(term,fmpq_mat_entry(left.a,i,j),fmpq_mat_entry(right.a,j,i));fmpq_add(target,target,term);}};
  fmpq_zero(four);for(int k=0;k<10;k++){fmpq_mul(term,fmpq_mat_entry(mu.a,k,0),fmpq_mat_entry(co.a,k,0));fmpq_add(four,four,term);}
  M e8(7,1);for(int i=0;i<6;i++)tr(w4[i],fmpq_mat_entry(e8.a,i,0));fmpq_mul(fmpq_mat_entry(e8.a,6,0),four,four);
  M values(1,width);
  for(int j=0;j<width;j++){
   const Def& d=defs[j];auto q=fmpq_mat_entry(values.a,0,j);
   if(d.type=='T')tr(d.word,q);
   else if(d.type=='F'){tr(d.word,q);fmpq_mul(q,q,four);}
   else{assert(d.type=='P');fmpq_mul(q,fmpq_mat_entry(e8.a,basis8[d.i],0),fmpq_mat_entry(e8.a,basis8[d.j],0));}
  }
  if(sample)out<<',';out<<'[';for(int j=0;j<width;j++){if(j)out<<',';out<<'"'<<text(fmpq_mat_entry(values.a,0,j))<<'"';}out<<']';
  for(int j=0;j<count;j++)fmpq_set(fmpq_mat_entry(minor.a,sample,j),fmpq_mat_entry(values.a,0,selected[j]));
 }
 fmpq_mat_det(sum,minor.a);assert(!fmpq_is_zero(sum));out<<"],\"exact_selected20_minor\":\""<<text(sum)<<"\",\"all_passed\":true}\n";
 cerr<<"EXACT_QQ_TRACE_INVARIANTS_AND_20_MINOR_PASS\n";
 fmpq_clear(half);fmpq_clear(term);fmpq_clear(sum);fmpq_clear(four);fmpq_clear(tmp);flint_cleanup();
}
