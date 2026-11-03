// Exact rank over QQ(t1,...,tn), with all row operations invertible.
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <gmp.h>
#include <flint/fmpz_mpoly_q.h>
using namespace std;
fmpz_mpoly_ctx_t ctx;
int variables;
struct Rational {
    fmpz_mpoly_q_struct v[1];
    Rational(){fmpz_mpoly_q_init(v,ctx);}
    Rational(const Rational& p):Rational(){fmpz_mpoly_q_set(v,p.v,ctx);}
    Rational(Rational&& p):Rational(){fmpz_mpoly_q_swap(v,p.v,ctx);}
    Rational& operator=(Rational p){fmpz_mpoly_q_swap(v,p.v,ctx);return *this;}
    ~Rational(){fmpz_mpoly_q_clear(v,ctx);}
    bool zero()const{return fmpz_mpoly_q_is_zero(v,ctx);}
};
void read_poly(istream& in,fmpz_mpoly_t p){
    long size;in>>size;assert(in&&size>=0);fmpz_mpoly_zero(p,ctx);
    fmpz_t z;fmpz_init(z);vector<ulong> powers(variables);
    for(long j=0;j<size;j++){
        string coefficient;in>>coefficient;for(auto& p:powers)in>>p;assert(in);
        assert(fmpz_set_str(z,coefficient.c_str(),10)==0);
        fmpz_mpoly_push_term_fmpz_ui(p,z,powers.data(),ctx);
    }
    fmpz_mpoly_sort_terms(p,ctx);fmpz_mpoly_combine_like_terms(p,ctx);fmpz_clear(z);
}
int main(int argc,char**argv){
    assert(argc==2);ifstream in(argv[1]);int rows,columns;in>>variables>>rows>>columns;
    assert(in&&variables>0&&rows>0&&columns>0);fmpz_mpoly_ctx_init(ctx,variables,ORD_DEGLEX);
    auto start=chrono::steady_clock::now();
    {
        vector<vector<Rational>> matrix(rows,vector<Rational>(columns));
        for(auto& row:matrix)for(auto& q:row){
            read_poly(in,fmpz_mpoly_q_numref(q.v));read_poly(in,fmpz_mpoly_q_denref(q.v));
            assert(!fmpz_mpoly_is_zero(fmpz_mpoly_q_denref(q.v),ctx));fmpz_mpoly_q_canonicalise(q.v,ctx);
        }
        int rank=0;vector<int> pivots;Rational factor,product;
        for(int column=0;column<columns&&rank<rows;column++){
            int pivot=rank;while(pivot<rows&&matrix[pivot][column].zero())pivot++;
            if(pivot==rows)continue;
            if(pivot!=rank)swap(matrix[pivot],matrix[rank]);
            // No row rescaling is needed: triangular elimination gives the
            // exact rank and avoids expanding a complete reduced normal form.
            for(int row=rank+1;row<rows;row++)if(!matrix[row][column].zero()){
                fmpz_mpoly_q_div(factor.v,matrix[row][column].v,matrix[rank][column].v,ctx);
                for(int c=column+1;c<columns;c++){
                    fmpz_mpoly_q_mul(product.v,factor.v,matrix[rank][c].v,ctx);
                    fmpz_mpoly_q_sub(matrix[row][c].v,matrix[row][c].v,product.v,ctx);
                }
                fmpz_mpoly_q_zero(matrix[row][column].v,ctx);
            }
            pivots.push_back(column);rank++;cout<<"PIVOT "<<column<<" RANK "<<rank<<endl;
        }
        for(int row=rank;row<rows;row++)for(int c=0;c<columns;c++)assert(matrix[row][c].zero());
        cout<<"EXACT_FUNCTION_FIELD_RANK "<<rank<<" ROWS "<<rows<<" COLUMNS "<<columns<<endl;
        cout<<"SECONDS "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;
    }
    fmpz_mpoly_ctx_clear(ctx);flint_cleanup();
}
