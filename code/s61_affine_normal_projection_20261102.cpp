#include <flint/nmod_mat.h>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <numeric>
#include <vector>

// Linearized relations have shape 129 x (42 normals) x (4 plane coordinates).
// Their products by the normal coordinates have bidegree (1,2).
int main(int argc, char **argv) {
    assert(argc == 3);
    std::ifstream input(argv[1]);
    std::ofstream output(argv[2]);
    constexpr slong count = 129, normals = 42, planes = 4, p = 32003;
    slong pairs[normals][normals], size = 0;
    for (slong i = 0; i < normals; ++i)
        for (slong j = i; j < normals; ++j)
            pairs[i][j] = pairs[j][i] = size++;
    const slong columns = size * planes, rows = count * normals;
    nmod_mat_t matrix, echelon, square, transpose, rhs, solution, check;
    nmod_mat_init(matrix, rows, columns, p);
    for (slong i = 0; i < count; ++i)
        for (slong k = 0; k < normals; ++k)
            for (slong j = 0; j < planes; ++j) {
                ulong u; input >> u; assert(input && u < p);
                for (slong l = 0; l < normals; ++l)
                    nmod_mat_entry(matrix, l * count+i, j * size+pairs[l][k]) = u;
            }
    std::string extra;
    assert(!(input >> extra));
    nmod_mat_init_set(echelon, matrix);
    std::vector<slong> perm(rows);
    std::iota(perm.begin(), perm.end(), 0);
    assert(nmod_mat_lu(perm.data(), echelon, 0) == columns);
    nmod_mat_clear(echelon);
    nmod_mat_init(square, columns, columns, p);
    for (slong i = 0; i < columns; ++i)
        for (slong j = 0; j < columns; ++j)
            nmod_mat_entry(square,i,j) = nmod_mat_entry(matrix, perm[i], j);
    nmod_mat_clear(matrix);
    ulong det = nmod_mat_det(square);
    assert(det);
    output << rows << ' ' << columns << ' ' << det << '\n';
    for (slong i = 0; i < columns; ++i) output << perm[i] << ' ';
    output << '\n';
    nmod_mat_init(transpose, columns, columns, p);
    nmod_mat_transpose(transpose, square);
    nmod_mat_clear(square);
    nmod_mat_init(rhs, columns, planes, p);
    for (slong j = 0; j < planes; ++j)
        nmod_mat_entry(rhs, j * size+pairs[0][0],j) = 1;
    nmod_mat_init(solution, columns, planes, p);
    assert(nmod_mat_solve(solution, transpose, rhs));
    nmod_mat_init(check, columns, planes, p);
    nmod_mat_mul(check, transpose, solution);
    assert(nmod_mat_equal(check, rhs));
    for (slong i = 0; i < columns; ++i) {
        for (slong j = 0; j < planes; ++j) output << nmod_mat_entry(solution,i,j) << ' ';
        output << '\n';
    }
    output.flush(); assert(output.good());
    nmod_mat_clear(check); nmod_mat_clear(solution); nmod_mat_clear(rhs); nmod_mat_clear(transpose);
    flint_cleanup();
}
