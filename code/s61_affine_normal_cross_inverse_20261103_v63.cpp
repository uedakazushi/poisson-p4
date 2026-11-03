#include <flint/fmpq_mat.h>
#include <flint/fmpq.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
using namespace std;
static const long p=32003;
long power(long x,long e){long y=1;for(;e;e>>=1,x=x*x%p)if(e&1)y=y*x%p;return y;}
void getq(istream&in,fmpq_t q){string x;if(!(in>>x))throw runtime_error("input EOF");if(fmpq_set_str(q,x.c_str(),10))throw runtime_error("bad rational");fmpq_canonicalise(q);}
void putq(ostream&out,const fmpq_t q){char*c=fmpq_get_str(nullptr,10,q);out<<c;flint_free(c);}
long residue(const fmpq_t q){long n=fmpz_fdiv_ui(fmpq_numref(q),p),d=fmpz_fdiv_ui(fmpq_denref(q),p);if(!d)throw runtime_error("nonunit input denominator");return n*power(d,p-2)%p;}
int main(int argc,char**argv){if(argc!=3)return 2;ifstream in(argv[1]);ofstream out(argv[2]);long nr,nc,nd;in>>nr>>nc>>nd;if(nr!=129||nc!=126||nd!=945)throw runtime_error("shape");fmpq_mat_t C,D,S,T,V,W,P,Q;fmpq_mat_init(C,nr,nc);fmpq_mat_init(D,nr,nd);fmpq_mat_init(S,nc,nc);fmpq_mat_init(T,nc,nd);fmpq_mat_init(V,nc,nc);fmpq_mat_init(W,nc,nc);fmpq_mat_init(P,nc,nd);fmpq_mat_init(Q,nc,nd);
 vector<vector<long>> a(nr,vector<long>(nc));for(long i=0;i<nr;i++)for(long j=0;j<nc;j++){getq(in,fmpq_mat_entry(C,i,j));a[i][j]=residue(fmpq_mat_entry(C,i,j));}for(long i=0;i<nr;i++)for(long j=0;j<nd;j++)getq(in,fmpq_mat_entry(D,i,j));string extra;if(in>>extra)throw runtime_error("extra");
 vector<long> ids(nr);for(long i=0;i<nr;i++)ids[i]=i;long rank=0;for(long j=0;j<nc;j++){long row=rank;while(row<nr&&!a[row][j])row++;if(row==nr)continue;swap(a[row],a[rank]);swap(ids[row],ids[rank]);long inv=power(a[rank][j],p-2);for(long k=j;k<nc;k++)a[rank][k]=a[rank][k]*inv%p;for(long i=rank+1;i<nr;i++){long c=a[i][j];if(c)for(long k=j;k<nc;k++)a[i][k]=(a[i][k]-c*a[rank][k]%p+p)%p;}rank++;}if(rank!=nc)throw runtime_error("rank not126");
 for(long i=0;i<nc;i++){for(long j=0;j<nc;j++)fmpq_set(fmpq_mat_entry(S,i,j),fmpq_mat_entry(C,ids[i],j));for(long j=0;j<nd;j++)fmpq_set(fmpq_mat_entry(T,i,j),fmpq_mat_entry(D,ids[i],j));}
 if(!fmpq_mat_inv(V,S))throw runtime_error("QQ inverse failed");fmpq_mat_mul(W,S,V);if(!fmpq_mat_is_one(W))throw runtime_error("S V notI");fmpq_mat_mul(W,V,S);if(!fmpq_mat_is_one(W))throw runtime_error("V S notI");for(long i=0;i<nc;i++)for(long j=0;j<nc;j++)residue(fmpq_mat_entry(V,i,j));
 fmpq_mat_mul(P,V,T);fmpq_mat_neg(P,P);fmpq_mat_mul(Q,S,P);fmpq_mat_add(Q,Q,T);if(!fmpq_mat_is_zero(Q))throw runtime_error("selected mixed reconstruction");
 out<<nc<<" "<<nd<<"\n";for(long i=0;i<nc;i++)out<<ids[i]<<" ";out<<"\n";for(long i=0;i<nc;i++){for(long j=0;j<nc;j++){putq(out,fmpq_mat_entry(V,i,j));out<<" ";}out<<"\n";}for(long i=0;i<nc;i++){for(long j=0;j<nd;j++){putq(out,fmpq_mat_entry(P,i,j));out<<" ";}out<<"\n";}out.close();cout<<"QQ rank126, both126x126 inverse identities, p-unitinverse, all126x945 mixed-coefficient identities PASS\n";
 fmpq_mat_clear(C);fmpq_mat_clear(D);fmpq_mat_clear(S);fmpq_mat_clear(T);fmpq_mat_clear(V);fmpq_mat_clear(W);fmpq_mat_clear(P);fmpq_mat_clear(Q);return 0;}
