#include <flint/nmod_mat.h>
#include <cassert>
#include <fstream>
#include <string>

int main(int argc, char **argv) {
    assert(argc == 3);
    std::ifstream input(argv[1]);
    std::ofstream output(argv[2]);
    assert(input && output);
    int n;
    ulong p;
    input >> n >> p;
    assert(n == 2814 && p == 32003);
    nmod_mat_t matrix;
    nmod_mat_init(matrix, n, n, p);
    for (int i=0; i<n; ++i) {
        int count;
        input >> count;
        for (int k=0; k<count; ++k) {
            int j;
            ulong value;
            input >> j >> value;
            assert(j>=0 && j<n && value>0 && value<p);
            assert(nmod_mat_entry(matrix,i,j)==0);
            nmod_mat_entry(matrix,i,j)=value;
        }
    }
    std::string extra;
    assert(input && !(input>>extra));
    const ulong determinant=nmod_mat_det(matrix);
    assert(determinant);
    output << determinant << '\n';
    output.flush();
    assert(output.good());
    nmod_mat_clear(matrix);
    flint_cleanup();
}
