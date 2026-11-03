#include <flint/nmod_mat.h>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
using namespace std;

int main(int argc,char**argv){
    assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);int blocks;in>>blocks;
    long totalW=0,totalFull=0;auto start=chrono::steady_clock::now();
    for(int block=0;block<blocks;block++){
        int weight,cols,rowsW,rowsG;in>>weight>>cols>>rowsW>>rowsG;
        nmod_mat_t W,A;nmod_mat_init(W,rowsW,cols,32003);nmod_mat_init(A,rowsW+rowsG,cols,32003);
        for(int row=0;row<rowsW+rowsG;row++){
            int length;in>>length;
            for(int term=0;term<length;term++){
                int col;unsigned long coefficient;in>>col>>coefficient;assert(col<cols&&coefficient<32003);
                nmod_mat_entry(A,row,col)=coefficient;
                if(row<rowsW)nmod_mat_entry(W,row,col)=coefficient;
            }
        }
        long rankW=nmod_mat_rref(W),rankFull=nmod_mat_rref(A);assert(rankFull>=rankW);
        out<<weight<<' '<<cols<<' '<<rowsW<<' '<<rowsG<<' '<<rankW<<' '<<rankFull<<' '<<rankFull-rankW<<'\n';out.flush();
        totalW+=rankW;totalFull+=rankFull;
        nmod_mat_clear(W);nmod_mat_clear(A);
    }
    out<<"TOTAL "<<totalW<<' '<<totalFull<<' '<<totalFull-totalW<<" SECONDS "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';
    flint_cleanup();
}
