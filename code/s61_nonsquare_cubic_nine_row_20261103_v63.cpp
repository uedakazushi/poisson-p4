#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/flint.h>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

static void putq(std::ostream& out, const fmpq_t x) {
    char* value = fmpq_get_str(nullptr, 10, x);
    out << value;
    flint_free(value);
}
static std::vector<slong> pivots(const fmpq_mat_t A) {
    fmpq_mat_t R;
    fmpq_mat_init(R, fmpq_mat_nrows(A), fmpq_mat_ncols(A));
    slong rank = fmpq_mat_rref(R, A);
    std::vector<slong> out;
    for (slong i = 0; i < rank; ++i) {
        slong j = 0;
        while (j < fmpq_mat_ncols(R) && fmpq_is_zero(fmpq_mat_entry(R, i, j))) ++j;
        assert(j < fmpq_mat_ncols(R));
        out.push_back(j);
    }
    fmpq_mat_clear(R);
    return out;
}
static void select(fmpq_mat_t B, const fmpq_mat_t A,
                   const std::vector<slong>& rows, const std::vector<slong>& cols) {
    assert(fmpq_mat_nrows(B) == (slong)rows.size());
    assert(fmpq_mat_ncols(B) == (slong)cols.size());
    for (slong i = 0; i < (slong)rows.size(); ++i)
        for (slong j = 0; j < (slong)cols.size(); ++j)
            fmpq_set(fmpq_mat_entry(B, i, j), fmpq_mat_entry(A, rows[i], cols[j]));
}
static void indices(std::ostream& out, const std::vector<slong>& v) {
    for (auto i : v) out << i << ' ';
    out << '\n';
}
int main(int argc, char** argv) {
    assert(argc == 3);
    std::ifstream in(argv[1]); std::ofstream out(argv[2]);
    slong n, m; in >> n >> m; assert(n == 57 && m == 91);
    std::vector<slong> heavy(7), light(84), prod(41), genuine(16);
    for (auto& i : heavy) in >> i;
    for (slong i=0; i<41; ++i) prod[i]=i;
    for (slong i=0; i<16; ++i) genuine[i]=41+i;
    for (slong j=0,k=0; j<91; ++j) {
        bool is_heavy=false; for (auto h : heavy) if (h==j) is_heavy=true;
        if (!is_heavy) light[k++]=j;
    }
    fmpq_mat_t A, P, Ht;
    fmpq_mat_init(A,n,m); fmpq_mat_init(P,41,84); fmpq_mat_init(Ht,7,16);
    std::string token;
    for (slong i=0; i<n; ++i) for (slong j=0; j<m; ++j) {
        in >> token; assert(!token.empty());
        assert(fmpq_set_str(fmpq_mat_entry(A,i,j), token.c_str(),10)==0);
        fmpq_canonicalise(fmpq_mat_entry(A,i,j));
    }
    assert(!(in >> token));
    select(P,A,prod,light);
    for (slong i=0; i<7; ++i) for (slong j=0; j<16; ++j)
        fmpq_set(fmpq_mat_entry(Ht,i,j),fmpq_mat_entry(A,genuine[j],heavy[i]));
    auto pc=pivots(P), gr=pivots(Ht);
    assert(pc.size()==41 && gr.size()==7);
    std::vector<slong> cols48, rows48=prod, rows9, cols43;
    for (auto j : pc) cols48.push_back(light[j]);
    for (auto j : heavy) cols48.push_back(j);
    for (auto i : gr) rows48.push_back(genuine[i]);
    for (auto i : genuine) {
        bool picked=false; for (auto j : rows48) if(i==j) picked=true;
        if(!picked) rows9.push_back(i);
    }
    for (auto j : light) {
        bool picked=false; for (auto k : cols48) if(j==k) picked=true;
        if(!picked) cols43.push_back(j);
    }
    assert(rows9.size()==9 && cols43.size()==43);
    fmpq_mat_t D, Di, E, F, G, tmp, S, HH, PP;
    fmpq_mat_init(D,48,48); fmpq_mat_init(Di,48,48);
    fmpq_mat_init(E,48,43); fmpq_mat_init(F,9,48);
    fmpq_mat_init(G,9,43); fmpq_mat_init(tmp,9,48); fmpq_mat_init(S,9,43);
    fmpq_mat_init(HH,7,7); fmpq_mat_init(PP,41,41);
    select(D,A,rows48,cols48); select(E,A,rows48,cols43);
    select(F,A,rows9,cols48); select(G,A,rows9,cols43);
    assert(fmpq_mat_inv(Di,D));
    fmpq_mat_mul(tmp,F,Di); fmpq_mat_mul(S,tmp,E); fmpq_mat_sub(S,G,S);
    std::vector<slong> p41(cols48.begin(),cols48.begin()+41), g7(rows48.begin()+41,rows48.end());
    select(PP,A,prod,p41); select(HH,A,g7,heavy);
    fmpq_t dp,dh,d48,ds,check;
    fmpq_init(dp);fmpq_init(dh);fmpq_init(d48);fmpq_init(ds);fmpq_init(check);
    fmpq_mat_det(dp,PP);fmpq_mat_det(dh,HH);fmpq_mat_det(d48,D);
    fmpq_mul(check,dp,dh);assert(fmpq_equal(check,d48));
    auto sc=pivots(S); assert(sc.size()==9);
    fmpq_mat_t S9;
    fmpq_mat_init(S9,9,9); std::vector<slong> r9;
    for(slong i=0;i<9;++i)r9.push_back(i);
    select(S9,S,r9,sc);fmpq_mat_det(ds,S9);assert(!fmpq_is_zero(ds));
    out << "EXACT_QQ_NINE_ROW_CRITERION 41 7 48 9 57\n";
    indices(out,p41);indices(out,g7);indices(out,rows9);indices(out,cols43);indices(out,sc);
    putq(out,dp);out << '\n';putq(out,dh);out << '\n';putq(out,d48);out << '\n';putq(out,ds);out << '\n';
    for(slong i=0;i<9;++i)for(slong j=0;j<43;++j){putq(out,fmpq_mat_entry(S,i,j));out << '\n';}
    std::cout << "EXACT_QQ_PRODUCT41_HEAVY7_SCHUR9_ALL387_ENTRIES_PASS\n";
}
