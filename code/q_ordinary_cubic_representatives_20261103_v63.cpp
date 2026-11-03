#include <flint/nmod_mat.h>
#include <flint/ulong_extras.h>
#include <cassert>
#include <fstream>
#include <map>
#include <vector>
using namespace std;
static const unsigned long P=32003;

int main(int argc,char**argv){
    assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);int blocks,total=0;in>>blocks;
    for(int block=0;block<blocks;block++){
        int weight,cols,nW,nG;in>>weight>>cols>>nW>>nG;
        nmod_mat_t W;nmod_mat_init(W,nW,cols,P);vector<vector<unsigned long>>G(nG,vector<unsigned long>(cols));
        for(int row=0;row<nW+nG;row++){
            int length;in>>length;
            for(int term=0;term<length;term++){
                int col;unsigned long c;in>>col>>c;
                if(row<nW)nmod_mat_entry(W,row,col)=c;else G[row-nW][col]=c;
            }
        }
        long rank=nmod_mat_rref(W);vector<int>Wpivots(rank);
        for(int row=0;row<rank;row++){
            int pivot=0;while(pivot<cols&&!nmod_mat_entry(W,row,pivot))pivot++;
            assert(pivot<cols&&nmod_mat_entry(W,row,pivot)==1);Wpivots[row]=pivot;
        }
        map<int,vector<unsigned long>>extension;
        for(int index=0;index<nG;index++){
            auto row=G[index];
            for(int b=0;b<rank;b++){
                auto factor=row[Wpivots[b]];if(!factor)continue;
                for(int col=Wpivots[b];col<cols;col++)row[col]=(row[col]+P-factor*nmod_mat_entry(W,b,col)%P)%P;
            }
            while(true){
                int pivot=0;while(pivot<cols&&!row[pivot])pivot++;if(pivot==cols)break;
                auto it=extension.find(pivot);
                if(it==extension.end()){
                    auto coefficient=row[pivot],inverse=n_invmod(coefficient,P);
                    for(auto&c:row)c=c*inverse%P;extension.emplace(pivot,row);
                    out<<weight<<' '<<index<<' '<<pivot<<' '<<coefficient<<'\n';total++;break;
                }
                auto factor=row[pivot];for(int col=pivot;col<cols;col++)row[col]=(row[col]+P-factor*it->second[col]%P)%P;
            }
        }
        out<<"BLOCK "<<weight<<' '<<rank<<' '<<extension.size()<<'\n';nmod_mat_clear(W);
    }
    assert(total==17);out<<"TOTAL "<<total<<'\n';flint_cleanup();
}
