// One exact QQ21minor of authorized24 degree22 products at existing points.
// Degree14 gradients use DIRECT dtheta/dH, independently of modular contractions.
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
M add(const M& a,const M& b){M c(a);fmpq_mat_add(c.a,a.a,b.a);return c;}
M read(ifstream& f){long n,m;f>>n>>m;assert(f&&n>0&&m>0);M a(n,m);string c;for(long i=0;i<n;i++)for(long j=0;j<m;j++){f>>c;assert(f&&fmpq_set_str(fmpq_mat_entry(a.a,i,j),c.c_str(),10)==0);fmpq_canonicalise(fmpq_mat_entry(a.a,i,j));}return a;}
string str(const fmpq_t a){char* p=fmpq_get_str(nullptr,10,a);string r(p);flint_free(p);return r;}
void tracepair(fmpq_t result,const M& A,const M& B){fmpq_zero(result);fmpq_t t;fmpq_init(t);for(int i=0;i<16;i++)for(int j=0;j<16;j++){fmpq_mul(t,fmpq_mat_entry(A.a,i,j),fmpq_mat_entry(B.a,j,i));fmpq_add(result,result,t);}fmpq_clear(t);}
void outputrow(ofstream& out,const M& A,int row,int width){out<<'[';for(int i=0;i<width;i++){if(i)out<<',';out<<'"'<<str(fmpq_mat_entry(A.a,row,i))<<'"';}out<<']';}
int main(int argc,char**argv){
 assert(argc==3);{
 ifstream in(argv[1]);ofstream out(argv[2]);int count,width;in>>count>>width;assert(count==21&&width==24);
 M O=read(in),Oi=read(in),Gi=read(in);vector<M> rep,Or;for(int i=0;i<10;i++){rep.push_back(read(in));Or.push_back(mul(O,rep.back()));}
 vector<string>w4(5),w9(9);for(auto&w:w4){in>>w;assert(w.size()==4);}for(auto&w:w9){in>>w;assert(w.size()==9);}
 vector<pair<int,int>>brackets(3);for(auto&ij:brackets){in>>ij.first>>ij.second;assert(1<=ij.first&&ij.first<ij.second&&ij.second<=3);}
 vector<int>selected(count);for(auto&i:selected){in>>i;assert(0<=i&&i<width);}
 vector<M> DV;for(int k=0;k<16;k++){M d(16,10);for(int j=0;j<10;j++)for(int i=0;i<16;i++)fmpq_set(fmpq_mat_entry(d.a,i,j),fmpq_mat_entry(Or[j].a,i,k));DV.push_back(d);}
 M eval(count,width),factors(count,18),minor(count,count);fmpq_t half,tmp,two,sum;fmpq_init(half);fmpq_set_si(half,1,2);fmpq_init(tmp);fmpq_init(two);fmpq_init(sum);
 for(int sample=0;sample<count;sample++){
  M b=read(in);assert(fmpq_mat_nrows(b.a)==16&&fmpq_mat_ncols(b.a)==1);M V(16,10);
  for(int k=0;k<10;k++){M v=mul(Or[k],b);for(int i=0;i<16;i++)fmpq_set(fmpq_mat_entry(V.a,i,k),fmpq_mat_entry(v.a,i,0));}
  M VT=trans(V),mu=trans(mul(trans(b),V));for(int k=0;k<10;k++)fmpq_mul(fmpq_mat_entry(mu.a,k,0),fmpq_mat_entry(mu.a,k,0),half);
  M co=mul(Gi,mu),theta;
  for(int i=0;i<16;i++)for(int j=0;j<16;j++)for(int k=0;k<10;k++){fmpq_mul(tmp,fmpq_mat_entry(co.a,k,0),fmpq_mat_entry(rep[k].a,i,j));fmpq_add(fmpq_mat_entry(theta.a,i,j),fmpq_mat_entry(theta.a,i,j),tmp);}
  M U=mul(V,Gi),UT=trans(U),H=add(theta,mul(mul(Oi,U),VT)),dc=mul(Gi,VT);
  map<string,M>table;M I;fmpq_mat_one(I.a);table.emplace("",I);
  for(int n=1;n<=5;n++)for(int z=0;z<(1<<n);z++){string w;for(int k=n-1;k>=0;k--)w+=char('0'+((z>>k)&1));table.emplace(w,mul(table.at(w.substr(0,n-1)),w.back()=='0'?theta:H));}
  M F4=mul(trans(mu),co);fmpq_set(fmpq_mat_entry(factors.a,sample,0),fmpq_mat_entry(F4.a,0,0));
  for(int k=0;k<5;k++)tracepair(fmpq_mat_entry(factors.a,sample,1+k),table.at(w4[k].substr(0,2)),table.at(w4[k].substr(2)));
  for(int k=0;k<9;k++)tracepair(fmpq_mat_entry(factors.a,sample,9+k),table.at(w9[k].substr(0,4)),table.at(w9[k].substr(4)));
  vector<M>K0,K1;for(int z=1;z<=3;z++){M A,B;string w=w4[z];for(int r=0;r<4;r++){string cyc=w.substr(r+1)+w.substr(0,r);M& K=w[r]=='0'?A:B;fmpq_mat_add(K.a,K.a,table.at(cyc).a);}K0.push_back(A);K1.push_back(B);}
  M grads(16,3);
  for(int k=0;k<16;k++){
   M dt;for(int i=0;i<16;i++)for(int j=0;j<16;j++)for(int r=0;r<10;r++){fmpq_mul(tmp,fmpq_mat_entry(dc.a,r,k),fmpq_mat_entry(rep[r].a,i,j));fmpq_add(fmpq_mat_entry(dt.a,i,j),fmpq_mat_entry(dt.a,i,j),tmp);}
   M dh=add(dt,mul(Oi,add(mul(DV[k],UT),mul(U,trans(DV[k])))));
   for(int z=0;z<3;z++){tracepair(sum,K0[z],dt);tracepair(two,K1[z],dh);fmpq_add(fmpq_mat_entry(grads.a,k,z),sum,two);}
  }
  M bracket=mul(mul(trans(grads),Oi),grads);
  for(int k=0;k<3;k++)fmpq_set(fmpq_mat_entry(factors.a,sample,6+k),fmpq_mat_entry(bracket.a,brackets[k].first-1,brackets[k].second-1));
  for(int k=0;k<9;k++)fmpq_mul(fmpq_mat_entry(eval.a,sample,k),fmpq_mat_entry(factors.a,sample,0),fmpq_mat_entry(factors.a,sample,9+k));
  for(int i=0;i<5;i++)for(int j=0;j<3;j++)fmpq_mul(fmpq_mat_entry(eval.a,sample,9+3*i+j),fmpq_mat_entry(factors.a,sample,1+i),fmpq_mat_entry(factors.a,sample,6+j));
  for(int k=0;k<count;k++)fmpq_set(fmpq_mat_entry(minor.a,sample,k),fmpq_mat_entry(eval.a,sample,selected[k]));
 }
 string extra;assert(!(in>>extra));fmpq_mat_det(sum,minor.a);assert(!fmpq_is_zero(sum));
 out<<"{\"all_passed\":true,\"evaluation\":[";for(int i=0;i<count;i++){if(i)out<<',';outputrow(out,eval,i,width);}out<<"],\"factors_F4_Inv8five_Inv14three_Inv18nine\":[";
 for(int i=0;i<count;i++){if(i)out<<',';outputrow(out,factors,i,18);}out<<"],\"exact_selected21_minor\":\""<<str(sum)<<"\"}\n";
 cerr<<"EXACT_ALL504_QQ_DEGREE22_PRODUCTS_DIRECT_GRADIENTS_AND_21_MINOR_PASS\n";
 fmpq_clear(half);fmpq_clear(tmp);fmpq_clear(two);fmpq_clear(sum);
 }flint_cleanup();
}
