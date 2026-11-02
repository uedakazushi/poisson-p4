#include <flint/nmod_mat.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <vector>

// Reconstruct and check an actual representation of w1^3 on the central P25.
int main(int argc,char**argv) {
    assert(argc==3); std::ifstream f(argv[1]);std::ofstream out(argv[2]);
    constexpr int n=26, rows=134, full=3276, size=3271, p=32003;
    std::vector<bool> omitted(full,false);
    for(int i=0;i<5;i++){int v;f>>v;assert(f&&v>=0&&v<full&&!omitted[v]);omitted[v]=true;}
    std::vector<int> selection(size), row_ids(rows*n,-1);
    for(int i=0;i<size;i++){f>>selection[i];assert(f&&row_ids[selection[i]]<0);row_ids[selection[i]]=i;}
    std::vector<int> triples(n*n*n), kept(full,-1);
    int k=0,c=0;
    for(int i=0;i<n;i++)for(int j=i;j<n;j++)for(int l=j;l<n;l++) {
        std::array<int,3> mon{i,j,l};
        do{triples[(mon[0]*n+mon[1])*n+mon[2]]=k;}while(std::next_permutation(mon.begin(),mon.end()));
        if(!omitted[k])kept[k]=c++;
        ++k;
    }
    assert(c==size&&k==full);
    nmod_mat_t M,A,rhs,coeff,check;
    nmod_mat_init(M,full,size,p);
    for(int r=0;r<rows;r++)for(int i=0;i<n;i++)for(int j=i;j<n;j++){
        ulong v;f>>v;assert(f&&v<p);
        for(int l=0;l<n;l++)if(row_ids[r*n+l]>=0)nmod_mat_entry(M,triples[(i*n+j)*n+l],row_ids[r*n+l])=v;
    }
    nmod_mat_init(A,size,size,p);
    for(int i=0;i<full;i++)if(kept[i]>=0)for(int j=0;j<size;j++)nmod_mat_entry(A,kept[i],j)=nmod_mat_entry(M,i,j);
    nmod_mat_init(rhs,size,1,p);
    int target=triples[(4*n+4)*n+4];assert(kept[target]>=0);
    nmod_mat_entry(rhs,kept[target],0)=1;
    nmod_mat_init(coeff,size,1,p);assert(nmod_mat_solve(coeff,A,rhs));
    nmod_mat_init(check,full,1,p);nmod_mat_mul(check,M,coeff);
    for(int i=0;i<full;i++)assert(nmod_mat_entry(check,i,0)==ulong(i==target));
    out<<target<<'\n';for(int i=0;i<size;i++)out<<nmod_mat_entry(coeff,i,0)<<' ';out<<'\n';
    out.flush();assert(out.good());
    nmod_mat_clear(M);nmod_mat_clear(A);nmod_mat_clear(rhs);nmod_mat_clear(coeff);nmod_mat_clear(check);flint_cleanup();
}
