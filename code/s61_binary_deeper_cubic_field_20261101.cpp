// Verify the quadratic-times-linear rank and remove all nuisance coefficients
// of the two tracked cubic classes. No floating point arithmetic is used.
#include "s61_binary_deeper_cubic_field_20261101.hpp"
#include <set>
using namespace cubic_field;

int main(int argc, char** argv) {
    assert(argc == 4);
    fmpq_mpoly_ctx_init(ctx, 2, ORD_LEX);
    flint_set_num_threads(2);
    {
        ifstream matrix_input(argv[1]);
        int row_count, column_count;
        matrix_input >> row_count >> column_count;
        assert(row_count == 216 && column_count == 165);
        vector<Row> active(row_count), basis;
        vector<int> pivots;
        for (auto& row : active) {
            int count;
            matrix_input >> count;
            for (int j = 0; j < count; ++j) {
                int column;
                matrix_input >> column;
                row.emplace(column, read_fraction(matrix_input));
            }
        }
        while (true) {
            int selected = -1, pivot = -1;
            long score = 0;
            for (int i = 0; i < row_count; ++i)
                for (const auto& [column, value] : active[i])
                    if (selected < 0 || value.size() < score) {
                        selected = i;
                        pivot = column;
                        score = value.size();
                    }
            if (selected < 0) break;
            Row row = std::move(active[selected]);
            Fraction leading = row[pivot];
            if (basis.size() < 114) {
                assert(fmpq_mpoly_is_fmpq(leading.numerator, ctx));
                assert(fmpq_mpoly_is_fmpq(leading.denominator, ctx));
            }
            for (auto& [column, value] : row) value = divide(value, leading);
            for (auto& other : active) eliminate(other, pivot, row);
            pivots.push_back(pivot);
            basis.push_back(std::move(row));
            if (basis.size() >= 114)
                cerr << "quadratic-linear rank " << basis.size() << '\n';
        }
        assert(basis.size() == 129);
        for (const auto& row : active) assert(row.empty());
        cerr << "exact quadratic-linear rank 129 PASS\n";

        map<vector<int>, int> cubic_index;
        int index = 0;
        for (int a = 0; a < 9; ++a)
            for (int b = a; b < 9; ++b)
                for (int c = b; c < 9; ++c) {
                    vector<int> monomial(9);
                    ++monomial[a]; ++monomial[b]; ++monomial[c];
                    cubic_index[monomial] = index++;
                }
        ifstream cubic_input(argv[2]);
        int cases;
        cubic_input >> cases;
        assert(cases == 2);
        vector<int> ground_sizes;
        for (int case_index = 0; case_index < cases; ++case_index) {
            int id;
            long count;
            cubic_input >> id >> count;
            assert(id == 92 || id == 134);
            map<int, Row> coefficients;
            fmpq_t value;
            fmpq_init(value);
            ulong core_exponents[2];
            for (long i = 0; i < count; ++i) {
                string rational;
                cubic_input >> rational >> core_exponents[0] >> core_exponents[1];
                assert(fmpq_set_str(value, rational.c_str(), 10) == 0);
                fmpq_canonicalise(value);
                int outer = 0, power = 1;
                for (int j = 0; j < 5; ++j) {
                    int exponent;
                    cubic_input >> exponent;
                    assert(exponent >= 0 && exponent <= 1);
                    outer += exponent * power;
                    power *= 3;
                }
                // The only possible nuisance monomials are q,h,k, each linear.
                assert(outer == 0 || outer == 3 || outer == 27 || outer == 81);
                vector<int> normal(9);
                for (int& exponent : normal) cubic_input >> exponent;
                assert(cubic_index.count(normal));
                fmpq_mpoly_push_term_fmpq_ui(
                    coefficients[outer][cubic_index[normal]].numerator,
                    value, core_exponents, ctx);
            }
            long denominator_terms;
            cubic_input >> denominator_terms;
            Fraction denominator;
            for (long i = 0; i < denominator_terms; ++i) {
                string rational;
                cubic_input >> rational >> core_exponents[0] >> core_exponents[1];
                assert(fmpq_set_str(value, rational.c_str(), 10) == 0);
                fmpq_canonicalise(value);
                for (int j = 0; j < 14; ++j) {
                    int exponent;
                    cubic_input >> exponent;
                    assert(exponent == 0);
                }
                fmpq_mpoly_push_term_fmpq_ui(denominator.numerator,
                                            value, core_exponents, ctx);
            }
            fmpq_mpoly_sort_terms(denominator.numerator, ctx);
            fmpq_mpoly_combine_like_terms(denominator.numerator, ctx);
            assert(!denominator.zero());
            for (auto& [outer, row] : coefficients) {
                for (auto& [column, coefficient] : row) {
                    fmpq_mpoly_sort_terms(coefficient.numerator, ctx);
                    fmpq_mpoly_combine_like_terms(coefficient.numerator, ctx);
                    fmpq_mpoly_set(coefficient.denominator, denominator.numerator, ctx);
                    coefficient.reduce();
                }
                // These are genuine multiples of the quadratic initial forms.
                // No backward reduction or conjectured normal form is needed.
                for (int i = 0; i < 114; ++i)
                    eliminate(row, pivots[i], basis[i]);
                if (outer != 0) assert(row.empty());
                else {
                    assert(!row.empty());
                    ground_sizes.push_back(static_cast<int>(row.size()));
                }
            }
            fmpq_clear(value);
            cerr << "all nuisance coefficients of cubic " << id << " vanish PASS\n";
        }
        assert(ground_sizes.size() == 2);
        ofstream report(argv[3]);
        assert(report);
        report << "{\"all_passed\":true,\"quadratic_linear_multiple_rank\":129,"
               << "\"constant_pivots_used_for_nuisance_removal\":114,"
               << "\"all_nuisance_coefficients_zero\":true,\"ground_sizes\":["
               << ground_sizes[0] << ',' << ground_sizes[1] << "]}\n";
    }
    fmpq_mpoly_ctx_clear(ctx);
    flint_cleanup();
    return 0;
}
