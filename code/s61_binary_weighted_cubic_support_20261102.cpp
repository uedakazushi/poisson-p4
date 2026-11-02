#include "s61_binary_deeper_quartic_field_20261101.hpp"
vector<Rat> factors;
void support(ostream&out,const fmpq_mpoly_t p){
    assert(!fmpq_mpoly_is_zero(p,ctx));Rat left,q;fmpq_mpoly_set(left.n,p,ctx);
    out<<'[';
    for(int i=0;i<int(factors.size());i++){
        int power=0;
        while(fmpq_mpoly_divides(q.n,left.n,factors[i].n,ctx)){fmpq_mpoly_swap(q.n,left.n,ctx);power++;}
        if(i)out<<',';out<<power;
    }
    assert(fmpq_mpoly_is_fmpq(left.n,ctx));out<<']';
}
int main(int argc,char**argv){
    assert(argc==6);fmpq_mpoly_ctx_init(ctx,2,ORD_LEX);
    {
        ifstream fac(argv[1]),piv(argv[2]),minor(argv[3]),dual(argv[4]);ofstream out(argv[5]);
        int nf;fac>>nf;assert(nf==30);factors.resize(nf);for(auto&f:factors)readpoly(fac,f.n);
        out<<"{\"all_passed\":true,\"factor_count\":30,\"quadratic_pivots\":[";
        int row,col,count=0;
        while(piv>>row>>col){
            Rat c;readpoly(piv,c.n);readpoly(piv,c.d);if(count)out<<',';
            out<<"{\"row\":"<<row<<",\"column\":"<<col<<",\"numerator_powers\":";support(out,c.n);
            out<<",\"denominator_powers\":";support(out,c.d);out<<'}';count++;
        }
        assert(count==48);out<<"],\"dual_denominators\":[";
        int nd,nc;dual>>nd>>nc;assert(nd==2&&nc==84);
        for(int i=0;i<2;i++){
            int freecol,n;dual>>freecol>>n;if(i)out<<',';out<<"{\"free_column\":"<<freecol<<",\"entries\":[";
            for(int j=0;j<n;j++){
                int column;Rat c;dual>>column;readpoly(dual,c.n);readpoly(dual,c.d);if(j)out<<',';
                out<<"{\"column\":"<<column<<",\"powers\":";support(out,c.d);out<<'}';
            }
            out<<"]}";
        }
        out<<"],\"pairing_denominators\":[";int a,b;minor>>a>>b;assert(a==2&&b==2);
        for(int i=0;i<4;i++){
            Rat c;minor>>a>>b;readpoly(minor,c.n);readpoly(minor,c.d);if(i)out<<',';support(out,c.d);
        }
        string token;minor>>token;assert(token=="DETERMINANT");Rat determinant;readpoly(minor,determinant.n);readpoly(minor,determinant.d);
        out<<"],\"pairing_determinant_numerator_powers\":";support(out,determinant.n);
        out<<",\"pairing_determinant_denominator_powers\":";support(out,determinant.d);out<<"}\n";
        cerr<<"SUPPORT:48 pivots, two duals and pairing determinant PASS\n";
    }
    factors.clear();fmpq_mpoly_ctx_clear(ctx);flint_cleanup();
}
