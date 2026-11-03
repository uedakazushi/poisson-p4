// Numeric specialization of the polynomial Macaulay input; find a cokernel basis.
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <fstream>
#include <iostream>
#include <cassert>

int main(int argc, char** argv) {
    assert(argc == 3);
    std::ifstream in(argv[1]);
    ulong p, unused0, unused1, value = std::stoul(argv[2]);
    int rows, columns, extra;
    in >> p >> rows >> columns >> extra >> unused0 >> unused1;
    assert(in && p == 32003 && extra == 0);
    nmod_mat_t matrix;
    nmod_mat_init(matrix, rows, columns, p);
    nmod_poly_t polynomial;
    nmod_poly_init(polynomial, p);
    for (int i = 0; i < rows; ++i) {
        int count;
        in >> count;
        for (int j = 0; j < count; ++j) {
            int column, size;
            in >> column >> size;
            assert(column >= 0 && column < columns && size > 0);
            nmod_poly_zero(polynomial);
            for (int k = 0; k < size; ++k) {
                ulong coefficient;
                in >> coefficient;
                assert(coefficient < p);
                nmod_poly_set_coeff_ui(polynomial, k, coefficient);
            }
            nmod_mat_entry(matrix, i, column) = nmod_poly_evaluate_nmod(polynomial, value);
        }
    }
    assert(in);
    int rank = nmod_mat_rref(matrix);
    std::cout << rank;
    int previous = -1;
    for (int i = 0; i < rank; ++i) {
        int pivot = 0;
        while (pivot < columns && !nmod_mat_entry(matrix, i, pivot)) ++pivot;
        assert(pivot > previous && pivot < columns);
        for (int j = previous + 1; j < pivot; ++j) std::cout << ' ' << j;
        previous = pivot;
    }
    for (int j = previous + 1; j < columns; ++j) std::cout << ' ' << j;
    std::cout << '\n';
    nmod_poly_clear(polynomial);
    nmod_mat_clear(matrix);
}
