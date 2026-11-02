// Exact positive-weight recursion over QQ[core, normals, rho].
#include <algorithm>
#include <cassert>
#include <flint/fmpq_mpoly.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
fmpq_mpoly_ctx_t ctx;
int variables;
struct Polynomial {
    fmpq_mpoly_struct v[1];
    Polynomial() { fmpq_mpoly_init(v,ctx); }
    Polynomial(const Polynomial& other):Polynomial() { fmpq_mpoly_set(v,other.v,ctx); }
    Polynomial(Polynomial&& other):Polynomial() { fmpq_mpoly_swap(v,other.v,ctx); }
    Polynomial& operator=(Polynomial other) { fmpq_mpoly_swap(v,other.v,ctx);return *this; }
    ~Polynomial() { fmpq_mpoly_clear(v,ctx); }
    bool zero() const { return fmpq_mpoly_is_zero(v,ctx); }
};
Polynomial product(const Polynomial& a,const Polynomial& b) {
    Polynomial result;fmpq_mpoly_mul(result.v,a.v,b.v,ctx);return result;
}
Polynomial generator(int j) { Polynomial p;fmpq_mpoly_gen(p.v,j,ctx);return p; }
void add(Polynomial& a,const Polynomial& b) { fmpq_mpoly_add(a.v,a.v,b.v,ctx); }

Polynomial read_core(istream& in) {
    int size;in>>size;assert(in&&size>=0);
    Polynomial p;fmpq_t c;fmpq_init(c);vector<ulong> exponent(variables,0);
    for(int i=0;i<size;i++) {
        string value;in>>value>>exponent[0];assert(in);
        assert(fmpq_set_str(c,value.c_str(),10)==0);fmpq_canonicalise(c);
        fmpq_mpoly_push_term_fmpq_ui(p.v,c,exponent.data(),ctx);
    }
    fmpq_mpoly_sort_terms(p.v,ctx);fmpq_mpoly_combine_like_terms(p.v,ctx);fmpq_clear(c);
    return p;
}
long write_row(ostream& out,const Polynomial& p,const vector<int>& weights,int weight) {
    long size=fmpq_mpoly_length(p.v,ctx);out<<weight<<" "<<size<<"\n";
    fmpq_t c;fmpq_init(c);vector<ulong> exponent(variables);
    for(long i=0;i<size;i++) {
        fmpq_mpoly_get_term_exp_ui(exponent.data(),p.v,i,ctx);
        int sum=0;for(int j=1;j<variables;j++)sum+=exponent[j]*weights[j-1];assert(sum==weight);
        fmpq_mpoly_get_term_coeff_fmpq(c,p.v,i,ctx);
        char* z=fmpq_get_str(nullptr,10,c);out<<z;flint_free(z);
        for(ulong d:exponent)out<<" "<<d;out<<"\n";
    }
    fmpq_clear(c);return size;
}
struct Term { int first,second;Polynomial coefficient; };
int main(int argc,char**argv) {
    assert(argc==3);ifstream in(argv[1]);ofstream out(argv[2]);assert(in&&out);
    int n,rank;in>>n>>rank;variables=n+1;
    fmpq_mpoly_ctx_init(ctx,variables,ORD_LEX);
    {
        vector<int> weights(45),normal_weights(n),selected(rank),row_weight(135);
        for(int& i:weights)in>>i;for(int& i:normal_weights)in>>i;for(int& i:selected)in>>i;
        vector<vector<Polynomial>> J(135,vector<Polynomial>(45)),T(45,vector<Polynomial>(n-1)),U(45,vector<Polynomial>(rank));
        for(auto* matrix:{&J,&T,&U})for(auto& row:*matrix)for(auto& v:row)v=read_core(in);
        vector<vector<Term>> equations(135);
        for(auto& row:equations) {
            int count;in>>count;assert(in&&count>=0);
            for(int k=0;k<count;k++) {
                int i,j,d;in>>i>>j>>d;assert((j>=0||d>0)&&d>=0);
                Polynomial coefficient=read_core(in),p;
                assert(fmpq_mpoly_pow_ui(p.v,generator(n).v,d,ctx));
                row.push_back({i,j,product(coefficient,p)});
            }
        }
        for(int& i:row_weight)in>>i;assert(in);
        vector<Polynomial> coordinates(45);
        for(int i=0;i<45;i++)for(int j=0;j<n-1;j++)if(!T[i][j].zero()) {
            assert(weights[i]==normal_weights[j]);
            add(coordinates[i],product(T[i][j],generator(j+1)));
        }
        auto nonlinear=[&](int row) {
            Polynomial result;
            for(const auto& t:equations[row]) {
                if(coordinates[t.first].zero())continue;
                if(t.second<0)add(result,product(t.coefficient,coordinates[t.first]));
                else if(!coordinates[t.second].zero())
                    add(result,product(t.coefficient,product(coordinates[t.first],coordinates[t.second])));
            }
            return result;
        };
        for(int d=1;d<=*max_element(weights.begin(),weights.end());d++) {
            vector<Polynomial> residuals(rank);
            for(int j=0;j<rank;j++)if(row_weight[selected[j]]==d)residuals[j]=nonlinear(selected[j]);
            for(int i=0;i<45;i++)if(weights[i]==d)for(int j=0;j<rank;j++)if(!U[i][j].zero()) {
                assert(row_weight[selected[j]]==d);auto term=product(U[i][j],residuals[j]);
                fmpq_mpoly_sub(coordinates[i].v,coordinates[i].v,term.v,ctx);
            }
        }
        out<<n<<" 136\n";long total=0;
        for(int i=0;i<135;i++) {
            Polynomial residual=nonlinear(i);
            for(int j=0;j<45;j++)if(!J[i][j].zero()&&!coordinates[j].zero())add(residual,product(J[i][j],coordinates[j]));
            if(find(selected.begin(),selected.end(),i)!=selected.end())assert(residual.zero());
            total+=write_row(out,residual,normal_weights,row_weight[i]);
        }
        long observable=write_row(out,coordinates[44],normal_weights,weights[44]);
        cout<<total<<" "<<observable<<"\n";assert(out);
    }
    fmpq_mpoly_ctx_clear(ctx);
}
