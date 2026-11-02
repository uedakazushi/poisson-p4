// Compute the determinant of a weak Popov row basis over F_p[b].
#include <flint/nmod_poly_mat.h>
#include <flint/nmod_poly.h>
#include <flint/ulong_extras.h>
#include <cassert>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

namespace {
ulong prime;
struct Polynomial {
    nmod_poly_t value;
    Polynomial() { nmod_poly_init(value, prime); }
    Polynomial(const Polynomial& p) : Polynomial() { nmod_poly_set(value, p.value); }
    Polynomial(Polynomial&& p) noexcept : Polynomial() { nmod_poly_swap(value, p.value); }
    Polynomial& operator=(Polynomial p) { nmod_poly_swap(value, p.value); return *this; }
    ~Polynomial() { nmod_poly_clear(value); }
};
using Row = std::vector<Polynomial>;

std::pair<int,int> leading(const Row& row) {
    int degree = -1, position = -1;
    for (int j = 0; j < static_cast<int>(row.size()); ++j) {
        const int current = nmod_poly_degree(row[j].value);
        if (current >= 0 && current >= degree) { degree = current; position = j; }
    }
    return {degree, position};
}
}

int main(int argc, char** argv) {
    assert(argc == 2);
    std::ifstream input(argv[1]);
    int rows, columns;
    input >> prime >> rows >> columns;
    assert(input && prime == 32003 && columns == 386 && rows > columns);
    std::vector<Row> basis(columns);
    std::vector<int> degrees(columns, -1);
    Polynomial multiple, shifted;
    long operations = 0;
    for (int i = 0; i < rows; ++i) {
        Row row(columns);
        int count;
        input >> count;
        for (int k = 0; k < count; ++k) {
            int col, size;
            input >> col >> size;
            assert(col >= 0 && col < columns && size > 0);
            for (int j = 0; j < size; ++j) {
                ulong value;
                input >> value;
                assert(value < prime);
                nmod_poly_set_coeff_ui(row[col].value, j, value);
            }
        }
        assert(input);
        while (true) {
            auto [degree,j] = leading(row);
            if (j < 0) break;
            if (basis[j].empty()) {
                basis[j] = std::move(row); degrees[j] = degree; break;
            }
            if (degree < degrees[j]) {
                std::swap(row, basis[j]); std::swap(degree, degrees[j]);
            }
            // Only nonzero field constants are inverted. Each output row remains a
            // polynomial combination of the input rows, even after an exchange.
            const ulong c = nmod_poly_get_coeff_ui(row[j].value, degree);
            const ulong pivot = nmod_poly_get_coeff_ui(basis[j][j].value, degrees[j]);
            const ulong factor = (c * n_invmod(pivot, prime)) % prime;
            for (int k = 0; k < columns; ++k) {
                nmod_poly_scalar_mul_nmod(multiple.value, basis[j][k].value, factor);
                nmod_poly_shift_left(shifted.value, multiple.value, degree-degrees[j]);
                nmod_poly_sub(row[k].value, row[k].value, shifted.value);
            }
            ++operations;
        }
    }
    nmod_poly_mat_t matrix;
    nmod_poly_mat_init(matrix, columns, columns, prime);
    int degree_bound = 0;
    for (int i = 0; i < columns; ++i) {
        assert(!basis[i].empty());
        degree_bound += degrees[i];
        for (int j = 0; j < columns; ++j)
            nmod_poly_set(nmod_poly_mat_entry(matrix,i,j), basis[i][j].value);
    }
    Polynomial determinant;
    nmod_poly_mat_det_interpolate(determinant.value, matrix);
    assert(nmod_poly_degree(determinant.value) == degree_bound);
    nmod_poly_mat_clear(matrix);
    std::cout << rows << ' ' << columns << ' ' << operations << ' ' << degree_bound;
    for (int j = 0; j <= degree_bound; ++j)
        std::cout << ' ' << nmod_poly_get_coeff_ui(determinant.value,j);
    std::cout << '\n';
}
