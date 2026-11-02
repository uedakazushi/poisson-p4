#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/nmod_mat.h>
#include <flint/ulong_extras.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

struct Entry { int column; long num,den; };
struct Row { std::string label; std::vector<Entry> entries; long scale=1; };

int main(int argc,char**argv) {
    assert(argc==3);
    const auto start=std::chrono::steady_clock::now();
    const auto seconds=[&](){return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
    const ulong p=32003;
    std::ifstream input(argv[1]);std::ofstream output(argv[2]);
    assert(input&&output);
    int nr,nt,nc;input>>nr>>nt>>nc;
    std::vector<Row> rows(nr+nt);
    for(auto& row:rows){
        int count;input>>row.label>>count;
        for(int k=0;k<count;++k){
            Entry term;input>>term.column>>term.num>>term.den;
            assert(term.column>=0&&term.column<nc&&term.den>0&&term.den%p);
            row.scale=std::lcm(row.scale,term.den);assert(row.scale>0&&row.scale<100000000L);
            row.entries.push_back(term);
        }
        for(auto& term:row.entries){term.num*=row.scale/term.den;term.den=1;}
    }
    std::string extra;assert(input&&!(input>>extra));
    std::vector<int> order(nr);std::iota(order.begin(),order.end(),0);
    // Relations L have no pure-beta correction. Select those first to shorten the lifts.
    const auto priority=[&](int i) {
        const char type=rows[i].label.front();
        const int group=type=='L' ? 0 : (type=='E' ? 1 : 2);
        return std::make_pair(group,rows[i].entries.size());
    };
    std::stable_sort(order.begin(),order.end(),[&](int a,int b){return priority(a)<priority(b);});
    nmod_mat_t mod;nmod_mat_init(mod,nr,nc,p);
    for(int i=0;i<nr;++i)for(auto term:rows[order[i]].entries)nmod_mat_entry(mod,i,term.column)=(term.num%(long)p+(long)p)%p;
    std::vector<slong> pr(nr);std::iota(pr.begin(),pr.end(),0);
    int rank=nmod_mat_lu(pr.data(),mod,0);nmod_mat_clear(mod);
    std::vector<int> generators(rank);for(int i=0;i<rank;++i)generators[i]=order[pr[i]];
    nmod_mat_t transpose;nmod_mat_init(transpose,nc,rank,p);
    for(int i=0;i<rank;++i)for(auto term:rows[generators[i]].entries)nmod_mat_entry(transpose,term.column,i)=(term.num%(long)p+(long)p)%p;
    std::vector<slong> pc(nc);std::iota(pc.begin(),pc.end(),0);
    assert(nmod_mat_lu(pc.data(),transpose,0)==rank);nmod_mat_clear(transpose);
    std::vector<int> projected(nc,-1);for(int i=0;i<rank;++i)projected[pc[i]]=i;
    fmpz_mat_t A,B,X;fmpz_mat_init(A,rank,rank);fmpz_mat_init(B,rank,nt);fmpz_mat_init(X,rank,nt);
    for(int j=0;j<rank;++j)for(auto term:rows[generators[j]].entries)if(projected[term.column]>=0)
        fmpz_set_si(fmpz_mat_entry(A,projected[term.column],j),term.num);
    for(int k=0;k<nt;++k)for(auto term:rows[nr+k].entries)if(projected[term.column]>=0)
        fmpz_set_si(fmpz_mat_entry(B,projected[term.column],k),term.num);
    fmpz_t den;fmpz_init(den);
    std::cerr<<"PIVOTS "<<rank<<" TARGETS "<<nt<<" SEC "<<seconds()<<'\n';
    assert(fmpz_mat_solve(X,den,A,B));assert(!fmpz_is_zero(den));fmpz_mat_clear(A);fmpz_mat_clear(B);
    std::cerr<<"SOLVED DEN_BITS "<<fmpz_bits(den)<<" SEC "<<seconds()<<'\n';
    // Check the original rectangular system, including every nonpivot column.
    fmpz_mat_t check;fmpz_mat_init(check,nc,nt);
    for(int j=0;j<rank;++j)for(int k=0;k<nt;++k){
        const fmpz* x=fmpz_mat_entry(X,j,k);
        if(fmpz_is_zero(x))continue;
        for(auto term:rows[generators[j]].entries) fmpz_addmul_si(fmpz_mat_entry(check,term.column,k),x,term.num);
    }
    for(int k=0;k<nt;++k)for(auto term:rows[nr+k].entries)fmpz_submul_si(fmpz_mat_entry(check,term.column,k),den,term.num);
    assert(fmpz_mat_is_zero(check));fmpz_mat_clear(check);
    output<<"EXACT "<<nr<<' '<<nt<<' '<<nc<<" PIVOTS "<<rank<<'\n';
    long terms=0;fmpq_t coefficient;fmpq_init(coefficient);
    for(int k=0;k<nt;++k){
        int count=0;for(int j=0;j<rank;++j)if(!fmpz_is_zero(fmpz_mat_entry(X,j,k)))++count;
        terms+=count;output<<"TARGET "<<rows[nr+k].label<<' '<<count<<'\n';
        for(int j=0;j<rank;++j)if(!fmpz_is_zero(fmpz_mat_entry(X,j,k))){
            fmpz_mul_si(fmpq_numref(coefficient),fmpz_mat_entry(X,j,k),rows[generators[j]].scale);
            fmpz_mul_si(fmpq_denref(coefficient),den,rows[nr+k].scale);fmpq_canonicalise(coefficient);
            char* value=fmpq_get_str(nullptr,10,coefficient);
            output<<rows[generators[j]].label<<' '<<value<<'\n';flint_free(value);
        }
    }
    output.flush();assert(output.good());
    std::cerr<<"FULL_QQ_IDENTITIES "<<nt<<" TERMS "<<terms<<" SEC "<<seconds()<<'\n';
    fmpq_clear(coefficient);fmpz_clear(den);fmpz_mat_clear(X);flint_cleanup();
}
