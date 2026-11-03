#include <flint/nmod_mat.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdexcept>
using namespace std;
int main(int argc,char**argv){if(argc!=2)return 2;ifstream in(argv[1]);long nr,nq,nl,nn;in>>nr>>nq>>nl>>nn;if(nr!=670||nl!=225)throw runtime_error("shape");nmod_mat_t Q,M,R,S;nmod_mat_init(Q,nr,nq,32003);nmod_mat_init(M,nr,nq+nl,32003);nmod_mat_init(R,nr,nq,32003);for(long k=0;k<nn;k++){long i,j,c;in>>i>>j>>c;if(i<0||i>=nr||j<0||j>=nq+nl||c<=0||c>=32003)throw runtime_error("entry");nmod_mat_entry(M,i,j)=c;if(j<nq)nmod_mat_entry(Q,i,j)=c;}nmod_mat_set(R,Q);long rq=nmod_mat_rref(R),rm=rq==nr?nr:nmod_mat_rank(M);if(rq!=670)throw runtime_error("notfullrank");vector<long> columns;for(long i=0;i<rq;i++){long j=0;while(j<nq&&!nmod_mat_entry(R,i,j))j++;if(j==nq)throw runtime_error("missingpivot");columns.push_back(j);}nmod_mat_init(S,rq,rq,32003);for(long i=0;i<rq;i++)for(long j=0;j<rq;j++)nmod_mat_entry(S,i,j)=nmod_mat_entry(Q,i,columns[j]);ulong d=nmod_mat_det(S);if(!d)throw runtime_error("minorzero");cout<<"RANK_QUADRATIC "<<rq<<" RANK_TOTAL "<<rm<<" MINOR_DETERMINANT "<<d<<"\nCOLUMNS";for(long j:columns)cout<<" "<<j;cout<<"\n";nmod_mat_clear(Q);nmod_mat_clear(M);nmod_mat_clear(R);nmod_mat_clear(S);return 0;}
