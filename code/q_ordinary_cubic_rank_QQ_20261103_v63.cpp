#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

int main(int argc,char**argv){
    assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);int blocks;in>>blocks;
    long totalW=0,totalFull=0;auto start=chrono::steady_clock::now();
    for(int block=0;block<blocks;block++){
        int weight,cols,rowsW,rowsG,modW,modFull;in>>weight>>cols>>rowsW>>rowsG>>modW>>modFull;
        fmpz_mat_t W,A;fmpz_mat_init(W,rowsW,cols);fmpz_mat_init(A,rowsW+rowsG,cols);
        for(int row=0;row<rowsW+rowsG;row++){
            int length;in>>length;
            for(int term=0;term<length;term++){
                int col;string coefficient;in>>col>>coefficient;assert(col<cols);
                assert(fmpz_set_str(fmpz_mat_entry(A,row,col),coefficient.c_str(),10)==0);
                if(row<rowsW)fmpz_set(fmpz_mat_entry(W,row,col),fmpz_mat_entry(A,row,col));
            }
        }
        // A modular full-row-rank certificate already proves the QQ rank in
        // these blocks; only rank-deficient cases require exact integer work.
        long rankW=modW==rowsW?rowsW:fmpz_mat_rank(W);
        long rankFull=rowsG?fmpz_mat_rank(A):rankW;
        assert(rankW>=modW&&rankFull>=modFull);
        out<<weight<<' '<<rankW<<' '<<rankFull<<' '<<rankFull-rankW<<'\n';out.flush();
        cerr<<"QQ BLOCK "<<weight<<" RANKS "<<rankW<<' '<<rankFull<<" SECONDS "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;
        totalW+=rankW;totalFull+=rankFull;fmpz_mat_clear(W);fmpz_mat_clear(A);
    }
    out<<"TOTAL "<<totalW<<' '<<totalFull<<' '<<totalFull-totalW<<" SECONDS "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';
    flint_cleanup();
}
