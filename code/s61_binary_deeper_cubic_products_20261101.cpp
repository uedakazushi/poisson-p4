// Exact polynomial certificate checks and cubic contraction over QQ.
// FLINT 3.5; input and output use sparse coefficients with sixteen exponents.
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <flint/fmpq_mpoly.h>

using namespace std;
constexpr int NVARS = 16;

void read_poly(istream& in, fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx) {
    fmpq_mpoly_zero(poly, ctx);
    long count;
    in >> count;
    assert(in && count >= 0);
    fmpq_t coefficient;
    fmpq_init(coefficient);
    ulong exponents[NVARS];
    for (long i = 0; i < count; ++i) {
        string value;
        in >> value;
        assert(fmpq_set_str(coefficient, value.c_str(), 10) == 0);
        fmpq_canonicalise(coefficient);
        for (int j = 0; j < NVARS; ++j) in >> exponents[j];
        assert(in);
        fmpq_mpoly_push_term_fmpq_ui(poly, coefficient, exponents, ctx);
    }
    fmpq_mpoly_sort_terms(poly, ctx);
    fmpq_mpoly_combine_like_terms(poly, ctx);
    fmpq_clear(coefficient);
}

void write_poly(ostream& out, const fmpq_mpoly_t poly,
                const fmpq_mpoly_ctx_t ctx) {
    const long count = fmpq_mpoly_length(poly, ctx);
    out << count << '\n';
    fmpq_t coefficient;
    fmpq_init(coefficient);
    ulong exponents[NVARS];
    for (long i = 0; i < count; ++i) {
        fmpq_mpoly_get_term_coeff_fmpq(coefficient, poly, i, ctx);
        char* value = fmpq_get_str(nullptr, 10, coefficient);
        out << value;
        flint_free(value);
        fmpq_mpoly_get_term_exp_ui(exponents, poly, i, ctx);
        for (int j = 0; j < NVARS; ++j) out << ' ' << exponents[j];
        out << '\n';
    }
    fmpq_clear(coefficient);
}

int main(int argc, char** argv) {
    assert(argc == 3);
    ifstream in(argv[1]);
    ofstream out(argv[2]);
    assert(in && out);
    flint_set_num_threads(2);
    fmpq_mpoly_ctx_t ctx;
    fmpq_mpoly_ctx_init(ctx, NVARS, ORD_LEX);
    fmpq_mpoly_t denominator, multiplier, product, residual;
    for (auto ptr : {denominator, multiplier, product, residual})
        fmpq_mpoly_init(ptr, ctx);
    int nbasis, nrows, nrelations;
    in >> nbasis >> nrows >> nrelations;
    assert(nbasis == 24 && nrows == 53 && nrelations == 2);
    read_poly(in, denominator, ctx);
    vector<fmpq_mpoly_struct> basis(nbasis), quadratic_rows(nrows);
    for (auto& row : basis) {
        fmpq_mpoly_init(&row, ctx);
        read_poly(in, &row, ctx);
    }
    map<int, int> row_index;
    for (int i = 0; i < nrows; ++i) {
        int id;
        in >> id;
        row_index[id] = i;
        fmpq_mpoly_init(&quadratic_rows[i], ctx);
        read_poly(in, &quadratic_rows[i], ctx);
        fmpq_mpoly_mul(residual, denominator, &quadratic_rows[i], ctx);
        for (int j = 0; j < nbasis; ++j) {
            read_poly(in, multiplier, ctx);
            fmpq_mpoly_mul(product, multiplier, &basis[j], ctx);
            fmpq_mpoly_sub(residual, residual, product, ctx);
        }
        assert(fmpq_mpoly_is_zero(residual, ctx));
    }
    cerr << "all 53 quadratic rowspace identities PASS\n";
    for (int i = 0; i < nrelations; ++i) {
        int count;
        in >> count;
        fmpq_mpoly_zero(residual, ctx);
        for (int j = 0; j < count; ++j) {
            int id;
            in >> id;
            assert(row_index.count(id));
            read_poly(in, multiplier, ctx);
            fmpq_mpoly_mul(product, multiplier,
                           &quadratic_rows[row_index[id]], ctx);
            fmpq_mpoly_add(residual, residual, product, ctx);
        }
        assert(fmpq_mpoly_is_zero(residual, ctx));
    }
    cerr << "both tracked quadratic-row syzygies PASS\n";
    int nlift;
    in >> nlift;
    assert(nlift == 17);
    vector<fmpq_mpoly_struct> lifts(nlift);
    for (auto& lift : lifts) {
        fmpq_mpoly_init(&lift, ctx);
        read_poly(in, &lift, ctx);
    }
    int cases;
    in >> cases;
    assert(cases == 2);
    out << cases << '\n';
    for (int i = 0; i < cases; ++i) {
        int id;
        in >> id;
        assert(id == 92 || id == 134);
        read_poly(in, denominator, ctx);
        fmpq_mpoly_zero(residual, ctx);
        for (int j = 0; j < nlift; ++j) {
            read_poly(in, multiplier, ctx);
            fmpq_mpoly_mul(product, multiplier, &lifts[j], ctx);
            fmpq_mpoly_add(residual, residual, product, ctx);
        }
        fmpq_mpoly_t gcd, numerator, reduced_denominator;
        for (auto ptr : {gcd, numerator, reduced_denominator})
            fmpq_mpoly_init(ptr, ctx);
        assert(fmpq_mpoly_gcd_cofactors(gcd, numerator, reduced_denominator,
                                       residual, denominator, ctx));
        out << id << '\n';
        write_poly(out, numerator, ctx);
        write_poly(out, reduced_denominator, ctx);
        out.flush();
        cerr << "tracked cubic " << id << " PASS\n";
        for (auto ptr : {gcd, numerator, reduced_denominator})
            fmpq_mpoly_clear(ptr, ctx);
    }
    for (auto& row : basis) fmpq_mpoly_clear(&row, ctx);
    for (auto& row : quadratic_rows) fmpq_mpoly_clear(&row, ctx);
    for (auto& row : lifts) fmpq_mpoly_clear(&row, ctx);
    for (auto ptr : {denominator, multiplier, product, residual})
        fmpq_mpoly_clear(ptr, ctx);
    fmpq_mpoly_ctx_clear(ctx);
    flint_cleanup();
    return 0;
}
