// Independent acceptance of a QQ polynomial dot product by global denominator
// clearing and fmpz_mpoly arithmetic. No discovery, division, or localization.
// CLI: one-row-native cofactor-native expected-one-polynomial-native.
#include <flint/fmpq_mpoly.h>
#include <flint/fmpz_mpoly.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cassert>
using namespace std;
fmpq_mpoly_ctx_t qc;
fmpz_mpoly_ctx_t zc;
struct Q {
    fmpq_mpoly_t p;
    Q(){fmpq_mpoly_init(p,qc);}
    Q(const Q& x):Q(){fmpq_mpoly_set(p,x.p,qc);}
    Q(Q&& x):Q(){fmpq_mpoly_swap(p,x.p,qc);}
    ~Q(){fmpq_mpoly_clear(p,qc);}
};
struct Z {
    fmpz_mpoly_t p;
    Z(){fmpz_mpoly_init(p,zc);}
    Z(const Z& x):Z(){fmpz_mpoly_set(p,x.p,zc);}
    Z(Z&& x):Z(){fmpz_mpoly_swap(p,x.p,zc);}
    ~Z(){fmpz_mpoly_clear(p,zc);}
};
void readq(istream& f,Q& q){
    long n;f>>n;assert(f&&n>=0);
    fmpq_t c;fmpq_init(c);ulong ex[3];
    for(long i=0;i<n;i++){
        string s;f>>s>>ex[0]>>ex[1]>>ex[2];
        assert(f&&fmpq_set_str(c,s.c_str(),10)==0);
        fmpq_canonicalise(c);assert(!fmpq_is_zero(c));
        fmpq_mpoly_push_term_fmpq_ui(q.p,c,ex,qc);
    }
    fmpq_mpoly_sort_terms(q.p,qc);fmpq_mpoly_combine_like_terms(q.p,qc);
    assert(fmpq_mpoly_length(q.p,qc)==n);fmpq_clear(c);
}
void denominator(fmpz_t d,const vector<Q>& v){
    fmpz_one(d);fmpz_t t;fmpz_init(t);
    for(const auto& q:v){fmpq_mpoly_get_denominator(t,q.p,qc);fmpz_lcm(d,d,t);}
    fmpz_clear(t);
}
void clear_denominator(Z& z,const Q& q,const fmpz_t d){
    fmpq_t c;fmpq_init(c);fmpz_t b;fmpz_init(b);ulong ex[3];
    for(long j=0;j<fmpq_mpoly_length(q.p,qc);j++){
        fmpq_mpoly_get_term_coeff_fmpq(c,q.p,j,qc);
        fmpq_mpoly_get_term_exp_ui(ex,q.p,j,qc);
        assert(fmpz_divisible(d,fmpq_denref(c)));
        fmpz_divexact(b,d,fmpq_denref(c));fmpz_mul(b,b,fmpq_numref(c));
        fmpz_mpoly_push_term_fmpz_ui(z.p,b,ex,zc);
    }
    fmpz_mpoly_sort_terms(z.p,zc);fmpz_mpoly_combine_like_terms(z.p,zc);
    fmpq_clear(c);fmpz_clear(b);
}
int main(int argc,char** argv){
    assert(argc==4);fmpq_mpoly_ctx_init(qc,3,ORD_LEX);fmpz_mpoly_ctx_init(zc,3,ORD_LEX);
    bool pass=false;
    {
        ifstream rf(argv[1]),lf(argv[2]),tf(argv[3]);long rows,cols,metadata,nc;
        rf>>rows>>cols>>metadata;lf>>nc;assert(rows==1&&metadata==0&&cols==nc);
        vector<Q> row(cols),cof(cols);for(auto& q:row)readq(rf,q);for(auto& q:cof)readq(lf,q);
        Q target;readq(tf,target);string extra;assert(!(rf>>extra)&&!(lf>>extra)&&!(tf>>extra));
        fmpz_t dr,dl,dt,scale;fmpz_init(dr);fmpz_init(dl);fmpz_init(dt);fmpz_init(scale);
        denominator(dr,row);denominator(dl,cof);fmpq_mpoly_get_denominator(dt,target.p,qc);
        assert(fmpz_is_one(dt));fmpz_mul(scale,dr,dl);
        Z sum,tmp,expected;clear_denominator(expected,target,scale);
        long active=0,terms=0;
        for(long j=0;j<cols;j++){
            if(fmpq_mpoly_is_zero(cof[j].p,qc))continue;
            active++;terms+=fmpq_mpoly_length(cof[j].p,qc);
            Z a,b;clear_denominator(a,row[j],dr);clear_denominator(b,cof[j],dl);
            fmpz_mpoly_mul(tmp.p,a.p,b.p,zc);fmpz_mpoly_add(sum.p,sum.p,tmp.p,zc);
        }
        pass=fmpz_mpoly_equal(sum.p,expected.p,zc);
        cout<<"{\"all_passed\":"<<(pass?"true":"false")
            <<",\"components\":"<<cols<<",\"active\":"<<active<<",\"cofactor_terms\":"<<terms
            <<",\"row_denominator_bits\":"<<fmpz_bits(dr)<<",\"cofactor_denominator_bits\":"<<fmpz_bits(dl)
            <<",\"integer_result_terms\":"<<fmpz_mpoly_length(sum.p,zc)<<"}\n";
        fmpz_clear(dr);fmpz_clear(dl);fmpz_clear(dt);fmpz_clear(scale);
    }
    fmpq_mpoly_ctx_clear(qc);fmpz_mpoly_ctx_clear(zc);flint_cleanup();return pass?0:3;
}
