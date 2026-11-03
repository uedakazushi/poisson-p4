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
    auto started=chrono::steady_clock::now();long total=0;
    for(int block=0;block<blocks;block++){
        long weight,rows,cols,expected;in>>weight>>rows>>cols>>expected;
        fmpz_mat_t A,B,check;fmpz_mat_init(A,cols,rows);fmpz_mat_init(B,rows,rows);
        for(long i=0;i<rows;i++){
            long count;in>>count;for(long j=0;j<count;j++){
                long col;string c;in>>col>>c;assert(col<cols);
                assert(fmpz_set_str(fmpz_mat_entry(A,col,i),c.c_str(),10)==0);
            }
        }
        long nullity=fmpz_mat_nullspace(B,A);assert(nullity==rows-expected);
        assert(fmpz_mat_rank(B)==nullity);
        fmpz_mat_init(check,cols,rows);fmpz_mat_mul(check,A,B);assert(fmpz_mat_is_zero(check));
        out<<weight<<' '<<rows<<' '<<nullity<<'\n';
        for(long j=0;j<nullity;j++){
            long count=0;for(long i=0;i<rows;i++)count+=!fmpz_is_zero(fmpz_mat_entry(B,i,j));out<<count;
            for(long i=0;i<rows;i++)if(!fmpz_is_zero(fmpz_mat_entry(B,i,j))){
                char*c=fmpz_get_str(nullptr,10,fmpz_mat_entry(B,i,j));out<<' '<<i<<' '<<c;flint_free(c);
            }out<<'\n';
        }out.flush();total+=nullity;
        cerr<<"KERNEL "<<weight<<" NULLITY "<<nullity<<" SECONDS "<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<endl;
        fmpz_mat_clear(A);fmpz_mat_clear(B);fmpz_mat_clear(check);
    }out<<"TOTAL "<<total<<'\n';flint_cleanup();
}
