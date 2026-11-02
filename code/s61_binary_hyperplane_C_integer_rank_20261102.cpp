#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <cassert>
#include <fstream>
#include <string>

int main(int argc,char **argv) {
    assert(argc==3);
    std::ifstream input(argv[1]);
    std::ofstream output(argv[2]);
    assert(input && output);
    for (int weight=2;weight<=4;++weight) {
        int w,rows,columns;
        input>>w>>rows>>columns;
        assert(w==weight&&rows>0&&columns>0);
        fmpz_mat_t matrix;
        fmpz_mat_init(matrix,rows,columns);
        for(int i=0;i<rows;++i) {
            int count;
            input>>count;
            for(int k=0;k<count;++k) {
                int j;
                std::string value;
                input>>j>>value;
                assert(j>=0&&j<columns);
                assert(!fmpz_set_str(fmpz_mat_entry(matrix,i,j),value.c_str(),10));
            }
        }
        assert(input);
        output<<weight<<' '<<rows<<' '<<columns<<' '<<fmpz_mat_rank(matrix)<<'\n';
        fmpz_mat_clear(matrix);
    }
    std::string extra;
    assert(!(input>>extra));
    output.flush();assert(output.good());
    flint_cleanup();
}
