#include <flint/nmod_mat.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <numeric>
#include <vector>

// Omit five monomials that independently evaluate on the known five points.
int main(int argc,char** argv) {
    assert(argc==3);
    std::ifstream f(argv[1]);
    std::ofstream output(argv[2]);
    constexpr slong n=26, rows=134, dim=3276, columns=3271, prime=32003;
    std::vector<bool> omitted(dim,false);
    for (int i=0;i<5;++i) {
        slong id; f >> id; assert(f && id>=0 && id<dim && !omitted[id]);
        omitted[id]=true;
    }
    std::vector<slong> index(n*n*n,-1);
    slong k=0, full=0;
    for(slong i=0;i<n;i++) for(slong j=i;j<n;j++) for(slong l=j;l<n;l++) {
        if(omitted[full++]) continue;
        std::array<slong,3> mon{i,j,l};
        do {index[(mon[0]*n+mon[1])*n+mon[2]]=k;} while(std::next_permutation(mon.begin(),mon.end()));
        ++k;
    }
    assert(k==columns && full==dim);
    nmod_mat_t M,E,S;
    nmod_mat_init(M,rows*n,columns,prime);
    for(slong r=0;r<rows;r++) for(slong i=0;i<n;i++) for(slong j=i;j<n;j++) {
        ulong value; f>>value;assert(f&&value<prime);
        for(slong l=0;l<n;l++) {
            slong col=index[(i*n+j)*n+l];
            if(col>=0) nmod_mat_entry(M,r*n+l,col)=value;
        }
    }
    std::string extra; assert(!(f>>extra));
    nmod_mat_init_set(E,M);
    std::vector<slong> perm(rows*n);
    std::iota(perm.begin(),perm.end(),0);
    assert(nmod_mat_lu(perm.data(),E,0)==columns);
    nmod_mat_clear(E);
    nmod_mat_init(S,columns,columns,prime);
    for(slong i=0;i<columns;++i) for(slong j=0;j<columns;++j)
        nmod_mat_entry(S,i,j)=nmod_mat_entry(M,perm[i],j);
    nmod_mat_clear(M);
    ulong det=nmod_mat_det(S);
    assert(det);
    output << rows*n << ' ' << dim << ' ' << columns << ' ' <<det<<'\n';
    for(slong i=0;i<columns;++i) output << perm[i] << ' ';
    output << '\n';output.flush();assert(output.good());
    nmod_mat_clear(S);flint_cleanup();
}
