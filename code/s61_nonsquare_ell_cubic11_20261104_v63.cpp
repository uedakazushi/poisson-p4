#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
using namespace std;
struct Mat {fmpq_mat_t a; Mat(slong n,slong m){fmpq_mat_init(a,n,m);} Mat(const Mat&b):Mat(fmpq_mat_nrows(b.a),fmpq_mat_ncols(b.a)){fmpq_mat_set(a,b.a);} ~Mat(){fmpq_mat_clear(a);}};
vector<int> range(int n){vector<int> v;for(int i=0;i<n;i++)v.push_back(i);return v;}
Mat sub(const Mat&A,const vector<int>&r,const vector<int>&c){Mat B(r.size(),c.size());for(int i=0;i<(int)r.size();i++)for(int j=0;j<(int)c.size();j++)fmpq_set(fmpq_mat_entry(B.a,i,j),fmpq_mat_entry(A.a,r[i],c[j]));return B;}
vector<int> pivcols(const Mat&A){Mat B(A);slong rank=fmpq_mat_rref(B.a,A.a);vector<int> p;for(int i=0;i<rank;i++){int j=0;while(j<fmpq_mat_ncols(B.a)&&fmpq_is_zero(fmpq_mat_entry(B.a,i,j)))j++;assert(j<fmpq_mat_ncols(B.a));p.push_back(j);}return p;}
vector<int> pivrows(const Mat&A){Mat T(fmpq_mat_ncols(A.a),fmpq_mat_nrows(A.a));fmpq_mat_transpose(T.a,A.a);return pivcols(T);}
string determinant(const Mat&A){fmpq_t d;fmpq_init(d);fmpq_mat_det(d,A.a);char*z=fmpq_get_str(nullptr,10,d);string out(z);flint_free(z);fmpq_clear(d);return out;}
void ints(ostream&o,const vector<int>&v){o<<'[';for(int i=0;i<(int)v.size();i++){if(i)o<<',';o<<v[i];}o<<']';}
void matrixjson(ostream&o,const Mat&A){o<<'[';for(int i=0;i<fmpq_mat_nrows(A.a);i++){if(i)o<<',';o<<'[';for(int j=0;j<fmpq_mat_ncols(A.a);j++){if(j)o<<',';char*z=fmpq_get_str(nullptr,10,fmpq_mat_entry(A.a,i,j));o<<'"'<<z<<'"';flint_free(z);}o<<']';}o<<']';}
void readmat(istream&i,Mat&A){int n,m;i>>n>>m;assert(n==fmpq_mat_nrows(A.a)&&m==fmpq_mat_ncols(A.a));for(int r=0;r<n;r++)for(int c=0;c<m;c++){string s;i>>s;assert(!fmpq_set_str(fmpq_mat_entry(A.a,r,c),s.c_str(),10));fmpq_canonicalise(fmpq_mat_entry(A.a,r,c));}}
int main(int argc,char**argv){assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);string tag;in>>tag;assert(tag=="ELLCOLON11_FROZEN_QQ_POINT");Mat Q(20,35),G(16,40),H7(16,7);readmat(in,Q);readmat(in,G);readmat(in,H7);string extra;assert(!(in>>extra));
 auto qr=pivrows(Q);assert(qr.size()==19);Mat QB=sub(Q,qr,range(35));auto qc=pivcols(QB);assert(qc.size()==19);vector<int> remc;for(int j=0;j<35;j++)if(find(qc.begin(),qc.end(),j)==qc.end())remc.push_back(j);assert(remc.size()==16);
 Mat QP=sub(QB,range(19),qc),Qi(19,19);assert(fmpq_mat_inv(Qi.a,QP.a));vector<int> hc={35,36,37,38,39};Mat H5=sub(G,range(16),hc);auto gr=pivrows(H5);assert(gr.size()==5);vector<int> rg;for(int i=0;i<16;i++)if(find(gr.begin(),gr.end(),i)==gr.end())rg.push_back(i);assert(rg.size()==11);
 Mat HP=sub(H5,gr,range(5)),Hi(5,5);assert(fmpq_mat_inv(Hi.a,HP.a));Mat GP=sub(G,gr,range(40)),F(5,40);fmpq_mat_mul(F.a,Hi.a,GP.a);Mat L=sub(H5,rg,range(5)),LF(11,40);fmpq_mat_mul(LF.a,L.a,F.a);Mat R=sub(G,rg,range(40));fmpq_mat_sub(R.a,R.a,LF.a);for(int i=0;i<11;i++)for(int j:hc)assert(fmpq_is_zero(fmpq_mat_entry(R.a,i,j)));
 Mat RP=sub(R,range(11),qc),RPi(11,19);fmpq_mat_mul(RPi.a,RP.a,Qi.a);Mat Qrem=sub(QB,range(19),remc),cor(11,16);fmpq_mat_mul(cor.a,RPi.a,Qrem.a);Mat S=sub(R,range(11),remc);fmpq_mat_sub(S.a,S.a,cor.a);auto sc=pivcols(S);int rank=sc.size();
 Mat full(36,40);for(int i=0;i<20;i++)for(int j=0;j<35;j++)fmpq_set(fmpq_mat_entry(full.a,i,j),fmpq_mat_entry(Q.a,i,j));for(int i=0;i<16;i++)for(int j=0;j<40;j++)fmpq_set(fmpq_mat_entry(full.a,20+i,j),fmpq_mat_entry(G.a,i,j));assert((int)pivcols(full).size()==24+rank);
 auto h7r=pivrows(H7);assert(h7r.size()==7);Mat h7sq=sub(H7,h7r,range(7));out<<"{\"quad_product_rank\":19,\"restricted_heavy_rank\":5,\"original_heavy_rank\":7,\"cubic11_rank\":"<<rank<<",\"full_restricted_rank\":"<<24+rank<<",\"quad_pivot_rows\":";ints(out,qr);out<<",\"quad_pivot_columns\":";ints(out,qc);out<<",\"quotient16_columns\":";ints(out,remc);out<<",\"quad19_minor\":\""<<determinant(QP)<<"\",\"restricted_heavy5_pivot_rows\":";ints(out,gr);out<<",\"remaining11_genuine_rows\":";ints(out,rg);out<<",\"heavy5_minor\":\""<<determinant(HP)<<"\",\"heavy7_pivot_rows\":";ints(out,h7r);out<<",\"heavy7_minor\":\""<<determinant(h7sq)<<"\",\"actual11x16\":";matrixjson(out,S);out<<",\"actual11_pivot_columns\":";ints(out,sc);if(rank==11){Mat SM=sub(S,range(11),sc);out<<",\"nonzero11_minor\":\""<<determinant(SM)<<"\"";vector<int> fr=qr,fc=qc;for(int i:gr)fr.push_back(20+i);for(int i:rg)fr.push_back(20+i);for(int j:hc)fc.push_back(j);for(int j:sc)fc.push_back(remc[j]);Mat full35=sub(full,fr,fc);out<<",\"full35_minor_rows\":";ints(out,fr);out<<",\"full35_minor_columns\":";ints(out,fc);out<<",\"full35_minor_determinant\":\""<<determinant(full35)<<"\"";}out<<"}\n";cerr<<"QUAD19_HEAVY5_CUBIC11_RANK "<<rank<<" FULL_RESTRICTED "<<24+rank<<'\n';return 0;}
