#include <flint/nmod_mat.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <fstream>
#include <map>
#include <vector>
using namespace std;
using U=unsigned long;static const U P=32003;
using Sparse=vector<pair<int,U>>;

int main(int argc,char**argv){
    assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);int nq,nc;in>>nq>>nc;assert(nq==57&&nc==17);
    vector<Sparse>L(35);for(int i=0;i<12;i++)for(int j=0;j<23;j++){long c;in>>c;L[i].push_back({j,(c%long(P)+P)%P});}
    for(int j=0;j<23;j++)L[12+j].push_back({j,1});
    vector<array<int,2>>pairs;int pairIndex[23][23];
    for(int i=0;i<23;i++)for(int j=i;j<23;j++){pairIndex[i][j]=pairIndex[j][i]=pairs.size();pairs.push_back({i,j});}
    vector<array<int,3>>triples;map<array<int,3>,int>tripleIndex;
    for(int i=0;i<23;i++)for(int j=i;j<23;j++)for(int k=j;k<23;k++){tripleIndex[{i,j,k}]=triples.size();triples.push_back({i,j,k});}
    vector<vector<Sparse>>pairImages(35,vector<Sparse>(35));
    for(int i=0;i<35;i++)for(int j=i;j<35;j++){
        vector<U>dense(276);for(auto[a,c]:L[i])for(auto[b,d]:L[j]){int index=pairIndex[a][b];dense[index]=(dense[index]+c*d)%P;}
        for(int k=0;k<276;k++)if(dense[k])pairImages[i][j].push_back({k,dense[k]});
    }
    map<array<int,3>,Sparse>cubicImages;
    auto cubicImage=[&](array<int,3>mon)->const Sparse&{
        auto it=cubicImages.find(mon);if(it!=cubicImages.end())return it->second;
        vector<U>dense(2300);
        for(auto[index,c]:pairImages[mon[0]][mon[1]])for(auto[k,d]:L[mon[2]]){
            array<int,3>target{pairs[index][0],pairs[index][1],k};sort(target.begin(),target.end());
            int column=tripleIndex.at(target);dense[column]=(dense[column]+c*d)%P;
        }
        Sparse result;for(int j=0;j<2300;j++)if(dense[j])result.push_back({j,dense[j]});
        return cubicImages.emplace(mon,std::move(result)).first->second;
    };
    nmod_mat_t Q,C,W,A;nmod_mat_init(Q,57,276,P);nmod_mat_init(C,17,2300,P);
    for(int row=0;row<57;row++){
        int length;in>>length;for(int term=0;term<length;term++){
            int i,j;U c;in>>i>>j>>c;for(auto[column,d]:pairImages[i][j])nmod_mat_entry(Q,row,column)=(nmod_mat_entry(Q,row,column)+c*d)%P;
        }
    }
    for(int row=0;row<17;row++){
        int length;in>>length;for(int term=0;term<length;term++){
            array<int,3>mon;U c;in>>mon[0]>>mon[1]>>mon[2]>>c;
            for(auto[column,d]:cubicImage(mon))nmod_mat_entry(C,row,column)=(nmod_mat_entry(C,row,column)+c*d)%P;
        }
    }
    out<<"POLYNOMIALS\n";
    for(int row=0;row<74;row++){
        int degree=row<57?2:3,columns=row<57?276:2300;long terms=0;
        for(int j=0;j<columns;j++)terms+=(row<57?nmod_mat_entry(Q,row,j):nmod_mat_entry(C,row-57,j))!=0;
        out<<row<<' '<<degree<<' '<<terms<<'\n';
        for(int j=0;j<columns;j++){
            U c=row<57?nmod_mat_entry(Q,row,j):nmod_mat_entry(C,row-57,j);if(!c)continue;
            out<<c;if(degree==2)for(int k:pairs[j])out<<' '<<k;else for(int k:triples[j])out<<' '<<k;out<<'\n';
        }
    }
    nmod_mat_init(W,57*23,2300,P);nmod_mat_init(A,57*23+17,2300,P);
    for(int q=0;q<57;q++)for(int variable=0;variable<23;variable++)for(int column=0;column<276;column++){
        U c=nmod_mat_entry(Q,q,column);if(!c)continue;
        array<int,3>mon{pairs[column][0],pairs[column][1],variable};sort(mon.begin(),mon.end());int target=tripleIndex.at(mon);
        nmod_mat_entry(W,23*q+variable,target)=c;nmod_mat_entry(A,23*q+variable,target)=c;
    }
    for(int row=0;row<17;row++)for(int column=0;column<2300;column++)nmod_mat_entry(A,57*23+row,column)=nmod_mat_entry(C,row,column);
    auto start=chrono::steady_clock::now();long rQ=nmod_mat_rref(Q),rW=nmod_mat_rref(W),rA=nmod_mat_rref(A);
    out<<"RANKS "<<rQ<<' '<<rW<<' '<<rA<<" SECONDS "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';
    nmod_mat_clear(Q);nmod_mat_clear(C);nmod_mat_clear(W);nmod_mat_clear(A);flint_cleanup();
}
