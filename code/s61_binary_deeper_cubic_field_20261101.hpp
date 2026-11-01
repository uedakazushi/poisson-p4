#pragma once
// Exact QQ(e,d) arithmetic. Every polynomial operation below uses FLINT.
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <flint/fmpq_mpoly.h>

namespace cubic_field {
using namespace std;
inline fmpq_mpoly_ctx_t ctx;

struct Fraction {
    fmpq_mpoly_struct numerator[1], denominator[1];
    Fraction() {
        fmpq_mpoly_init(numerator, ctx);
        fmpq_mpoly_init(denominator, ctx);
        fmpq_mpoly_one(denominator, ctx);
    }
    Fraction(const Fraction& other) : Fraction() {
        fmpq_mpoly_set(numerator, other.numerator, ctx);
        fmpq_mpoly_set(denominator, other.denominator, ctx);
    }
    Fraction(Fraction&& other) : Fraction() {
        fmpq_mpoly_swap(numerator, other.numerator, ctx);
        fmpq_mpoly_swap(denominator, other.denominator, ctx);
    }
    Fraction& operator=(Fraction other) {
        fmpq_mpoly_swap(numerator, other.numerator, ctx);
        fmpq_mpoly_swap(denominator, other.denominator, ctx);
        return *this;
    }
    ~Fraction() {
        fmpq_mpoly_clear(numerator, ctx);
        fmpq_mpoly_clear(denominator, ctx);
    }
    bool zero() const { return fmpq_mpoly_is_zero(numerator, ctx); }
    long size() const {
        return fmpq_mpoly_length(numerator, ctx)
             + fmpq_mpoly_length(denominator, ctx);
    }
    void make_denominator_monic() {
        if (zero()) {
            fmpq_mpoly_one(denominator, ctx);
            return;
        }
        fmpq_t leading;
        fmpq_init(leading);
        fmpq_mpoly_get_term_coeff_fmpq(leading, denominator, 0, ctx);
        fmpq_mpoly_scalar_div_fmpq(numerator, numerator, leading, ctx);
        fmpq_mpoly_scalar_div_fmpq(denominator, denominator, leading, ctx);
        fmpq_clear(leading);
    }
    void reduce() {
        if (zero()) {
            fmpq_mpoly_one(denominator, ctx);
            return;
        }
        fmpq_mpoly_t gcd, a, b;
        for (auto ptr : {gcd, a, b}) fmpq_mpoly_init(ptr, ctx);
        assert(fmpq_mpoly_gcd_cofactors(gcd, a, b, numerator, denominator, ctx));
        fmpq_mpoly_swap(numerator, a, ctx);
        fmpq_mpoly_swap(denominator, b, ctx);
        for (auto ptr : {gcd, a, b}) fmpq_mpoly_clear(ptr, ctx);
        make_denominator_monic();
    }
};

inline Fraction multiply(const Fraction& a, const Fraction& b) {
    Fraction result;
    if (a.zero() || b.zero()) return result;
    Fraction cross_ab, cross_ba;
    fmpq_mpoly_t gcd;
    fmpq_mpoly_init(gcd, ctx);
    assert(fmpq_mpoly_gcd_cofactors(gcd, cross_ab.numerator, cross_ab.denominator,
                                   a.numerator, b.denominator, ctx));
    assert(fmpq_mpoly_gcd_cofactors(gcd, cross_ba.numerator, cross_ba.denominator,
                                   b.numerator, a.denominator, ctx));
    fmpq_mpoly_mul(result.numerator, cross_ab.numerator, cross_ba.numerator, ctx);
    fmpq_mpoly_mul(result.denominator, cross_ab.denominator, cross_ba.denominator, ctx);
    fmpq_mpoly_clear(gcd, ctx);
    result.make_denominator_monic();
    return result;
}

inline Fraction divide(const Fraction& a, const Fraction& b) {
    assert(!b.zero());
    Fraction inverse;
    fmpq_mpoly_set(inverse.numerator, b.denominator, ctx);
    fmpq_mpoly_set(inverse.denominator, b.numerator, ctx);
    return multiply(a, inverse);
}

inline Fraction subtract(const Fraction& a, const Fraction& b) {
    if (b.zero()) return a;
    Fraction result;
    if (a.zero()) {
        result = b;
        fmpq_mpoly_neg(result.numerator, result.numerator, ctx);
        return result;
    }
    // If denominators are g*A and g*B with gcd(A,B)=1, form
    // (a_num*B - b_num*A)/(g*A*B). Only a factor of g can then cancel.
    Fraction cofactors;
    fmpq_mpoly_t gcd, cancellation, temporary, reduced_num, reduced_gcd;
    for (auto ptr : {gcd, cancellation, temporary, reduced_num, reduced_gcd})
        fmpq_mpoly_init(ptr, ctx);
    assert(fmpq_mpoly_gcd_cofactors(gcd, cofactors.numerator, cofactors.denominator,
                                   a.denominator, b.denominator, ctx));
    fmpq_mpoly_mul(result.numerator, a.numerator, cofactors.denominator, ctx);
    fmpq_mpoly_mul(temporary, b.numerator, cofactors.numerator, ctx);
    fmpq_mpoly_sub(result.numerator, result.numerator, temporary, ctx);
    if (!result.zero()) {
        assert(fmpq_mpoly_gcd_cofactors(cancellation, reduced_num, reduced_gcd,
                                       result.numerator, gcd, ctx));
        fmpq_mpoly_swap(result.numerator, reduced_num, ctx);
        fmpq_mpoly_mul(result.denominator, cofactors.numerator, cofactors.denominator, ctx);
        fmpq_mpoly_mul(result.denominator, result.denominator, reduced_gcd, ctx);
        result.make_denominator_monic();
    }
    for (auto ptr : {gcd, cancellation, temporary, reduced_num, reduced_gcd})
        fmpq_mpoly_clear(ptr, ctx);
    return result;
}

inline void read_poly(istream& in, fmpq_mpoly_t poly) {
    fmpq_mpoly_zero(poly, ctx);
    int count;
    in >> count;
    assert(in && count >= 0);
    fmpq_t value;
    fmpq_init(value);
    ulong exponents[2];
    for (int i = 0; i < count; ++i) {
        string coefficient;
        in >> coefficient >> exponents[0] >> exponents[1];
        assert(in && fmpq_set_str(value, coefficient.c_str(), 10) == 0);
        fmpq_canonicalise(value);
        fmpq_mpoly_push_term_fmpq_ui(poly, value, exponents, ctx);
    }
    fmpq_mpoly_sort_terms(poly, ctx);
    fmpq_mpoly_combine_like_terms(poly, ctx);
    fmpq_clear(value);
}

inline Fraction read_fraction(istream& in) {
    Fraction value;
    read_poly(in, value.numerator);
    read_poly(in, value.denominator);
    assert(!fmpq_mpoly_is_zero(value.denominator, ctx));
    value.reduce();
    return value;
}
using Row = map<int, Fraction>;

inline void eliminate(Row& row, int pivot, const Row& basis_row) {
    auto found = row.find(pivot);
    if (found == row.end()) return;
    Fraction multiplier = found->second;
    for (const auto& [column, coefficient] : basis_row) {
        Fraction result = subtract(row[column], multiply(multiplier, coefficient));
        if (result.zero()) row.erase(column);
        else row[column] = std::move(result);
    }
}
}  // namespace cubic_field
