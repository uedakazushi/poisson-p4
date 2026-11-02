#include "s61_binary_deeper_quartic_field_20261101.hpp"
vector<Rat> known;
void support(ostream&out,const fmpq_mpoly_t p){
    Rat a,q;fmpq_mpoly_set(a.n,p,ctx);assert(!a.zero());out<<'[';
    for(int i=0;i<int(known.size());i++){
        int power=0;while(fmpq_mpoly_divides(q.n,a.n,known[i].n,ctx)){fmpq_mpoly_swap(q.n,a.n,ctx);power++;}
        if(i)out<<',';out<<power;
    }
    assert(fmpq_mpoly_is_fmpq(a.n,ctx));out<<']';
}
int main(int argc,char**argv){
    assert(argc==4);fmpq_mpoly_ctx_init(ctx,2,ORD_LEX);
    {
        ifstream fac(argv[1]),in(argv[2]);ofstream out(argv[3]);int n;fac>>n;assert(n==13);
        known.resize(n);for(auto&f:known)readpoly(fac,f.n);
        int nr,nc;in>>nr>>nc;assert(nr==169&&nc==46);
        out<<"{\"all_passed\":true,\"rref_nonzero_coefficients\":"<<nr<<",\"relation_nonzero_coefficients\":"<<nc<<",\"denominator_factor_exponents\":[";
        for(int i=0;i<nr+nc;i++){
            Rat c;readpoly(in,c.n);readpoly(in,c.d);c.reduce();if(i)out<<',';support(out,c.d);
        }
        out<<"],\"specialized_quadratic_contents\":[";in>>n;assert(n==2);
        for(int i=0;i<2;i++){
            int id;in>>id;assert(id==(i?134:92));Rat c;readpoly(in,c.n);if(i)out<<',';
            out<<"{\"row\":"<<id<<",\"factor_exponents\":";support(out,c.n);out<<'}';
        }
        out<<"]}\n";cerr<<"215 coefficient denominators and two raw A2 contents: D-supported PASS\n";
    }
    known.clear();fmpq_mpoly_ctx_clear(ctx);flint_cleanup();
}
