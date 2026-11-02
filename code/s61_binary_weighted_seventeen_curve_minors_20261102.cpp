#include "s61_binary_deeper_quartic_field_20261101.hpp"
#include <flint/fmpq_mpoly_factor.h>
#include <flint/nmod_mat.h>
#include <flint/nmod.h>
#include <flint/ulong_extras.h>
#include <array>
#include <chrono>
#include <numeric>
using Mon=array<int,9>;
int weights[9]={1,1,1,1,1,2,1,2,1};
void monomials(vector<Mon>&a,int d,int pos=0,Mon m={}){
    if(pos==9){if(!d)a.push_back(m);return;}
    for(int j=0;j*weights[pos]<=d;j++){m[pos]=j;monomials(a,d-j*weights[pos],pos+1,m);}
}
ulong evaluate(const fmpq_mpoly_t f,ulong a,ulong b,nmod_t mod){
    ulong x[2]={a,b};
    ulong v=fmpz_mpoly_evaluate_all_nmod(f->zpoly,x,ctx->zctx,mod);
    ulong num=fmpz_fdiv_ui(fmpq_numref(f->content),mod.n);
    ulong den=fmpz_fdiv_ui(fmpq_denref(f->content),mod.n);assert(den);
    return nmod_mul(v,nmod_mul(num,n_invmod(den,mod.n),mod),mod);
}
void array_json(ostream&out,const vector<int>&a){out<<'[';for(size_t i=0;i<a.size();i++){if(i)out<<',';out<<a[i];}out<<']';}
void minor_certificate(ostream&out,const nmod_mat_t M,int expected){
    int nr=nmod_mat_nrows(M),nc=nmod_mat_ncols(M),rank=0;
    nmod_mat_t E;nmod_mat_init_set(E,M);
    vector<int>labels(nr),rows,cols;iota(labels.begin(),labels.end(),0);
    for(int c=0;c<nc&&rank<nr;c++){
        int choice=rank;while(choice<nr&&!nmod_mat_entry(E,choice,c))choice++;
        if(choice==nr)continue;
        if(choice!=rank){for(int j=c;j<nc;j++)swap(nmod_mat_entry(E,rank,j),nmod_mat_entry(E,choice,j));swap(labels[rank],labels[choice]);}
        rows.push_back(labels[rank]);cols.push_back(c);
        ulong inv=n_invmod(nmod_mat_entry(E,rank,c),M->mod.n);
        for(int r=rank+1;r<nr;r++)if(nmod_mat_entry(E,r,c)){
            ulong scale=nmod_mul(nmod_mat_entry(E,r,c),inv,M->mod);
            for(int j=c;j<nc;j++)nmod_mat_entry(E,r,j)=nmod_sub(nmod_mat_entry(E,r,j),nmod_mul(scale,nmod_mat_entry(E,rank,j),M->mod),M->mod);
        }
        rank++;
    }
    assert(rank==expected&&nmod_mat_rank(M)==expected);
    nmod_mat_t square;nmod_mat_init(square,rank,rank,M->mod.n);
    for(int i=0;i<rank;i++)for(int j=0;j<rank;j++)nmod_mat_entry(square,i,j)=nmod_mat_entry(M,rows[i],cols[j]);
    ulong determinant=nmod_mat_det(square);assert(determinant);
    out<<"{\"rank\":"<<rank<<",\"matrix_rows\":"<<nr<<",\"matrix_columns\":"<<nc<<",\"minor_rows\":";
    array_json(out,rows);out<<",\"minor_columns\":";array_json(out,cols);out<<",\"determinant\":"<<determinant<<'}';
    nmod_mat_clear(square);nmod_mat_clear(E);
}
int main(int argc,char**argv){
    assert(argc==5);fmpq_mpoly_ctx_init(ctx,2,ORD_LEX);
    {
        ifstream in(argv[1]),pt(argv[2]),fs(argv[3]);ofstream out(argv[4]);
        int nr;in>>nr;assert(nr==61);vector<int>id(nr),wt(nr);vector<map<Mon,Rat>>rows(nr);
        for(int i=0;i<nr;i++){
            int n;in>>id[i]>>wt[i]>>n;
            for(int j=0;j<n;j++){Mon m;for(auto&a:m)in>>a;readpoly(in,rows[i][m].n);}
        }
        int nf,nk;fs>>nf;assert(nf==17);vector<Rat>factors(nf);for(auto&f:factors)readpoly(fs,f.n);
        fs>>nk;assert(nk==13);vector<Rat>known(nk);for(auto&f:known)readpoly(fs,f.n);
        for(int i=0;i<nf;i++){
            fmpq_mpoly_factor_t ff;fmpq_mpoly_factor_init(ff,ctx);assert(fmpq_mpoly_factor(ff,factors[i].n,ctx));
            assert(ff->num==1&&fmpq_mpoly_factor_get_exp_si(ff,0,ctx)==1);fmpq_mpoly_factor_clear(ff,ctx);
            for(int j=0;j<nk+i;j++){
                Rat g;assert(fmpq_mpoly_gcd(g.n,factors[i].n,j<nk?known[j].n:factors[j-nk].n,ctx));
                assert(fmpq_mpoly_is_fmpq(g.n,ctx));
            }
        }
        int np;pt>>np;assert(np==17);auto start=chrono::steady_clock::now();out<<"{\"all_passed\":true,\"QQ_irreducible_distinct_and_coprime_to_D\":true,\"cases\":[";
        for(int k=0;k<np;k++){
            int label;ulong prime,a,b;pt>>label>>prime>>a>>b;assert(label==k&&n_is_prime(prime));
            nmod_t mod;nmod_init(&mod,prime);assert(evaluate(factors[k].n,a,b,mod)==0);
            Rat de,dd;fmpq_mpoly_derivative(de.n,factors[k].n,0,ctx);fmpq_mpoly_derivative(dd.n,factors[k].n,1,ctx);
            ulong ge=evaluate(de.n,a,b,mod),gd=evaluate(dd.n,a,b,mod);assert(ge||gd);
            for(auto&f:known)assert(evaluate(f.n,a,b,mod));
            vector<map<Mon,ulong>>values(nr);for(int i=0;i<nr;i++)for(auto&[m,c]:rows[i])values[i][m]=evaluate(c.n,a,b,mod);
            if(k)out<<',';out<<"{\"id\":"<<label<<",\"prime\":"<<prime<<",\"point\":["<<a<<','<<b<<"],\"gradient\":["<<ge<<','<<gd<<"],\"ranks\":{";
            for(int degree:{2,3,4}){
                vector<Mon>cols;monomials(cols,degree);map<Mon,int>index;for(int j=0;j<int(cols.size());j++)index[cols[j]]=j;
                vector<pair<int,Mon>>source;
                for(int i=0;i<nr;i++)if(wt[i]<=degree){vector<Mon>shifts;monomials(shifts,degree-wt[i]);for(auto&s:shifts)source.emplace_back(i,s);}
                nmod_mat_t M;nmod_mat_init(M,source.size(),cols.size(),prime);
                for(int r=0;r<int(source.size());r++){auto[i,s]=source[r];for(auto&[m,c]:values[i]){Mon mon=m;for(int j=0;j<9;j++)mon[j]+=s[j];nmod_mat_entry(M,r,index.at(mon))=c;}}
                if(degree!=2)out<<',';out<<'"'<<degree<<"\":";minor_certificate(out,M,degree==2?7:degree==3?64:269);nmod_mat_clear(M);
            }
            out<<"}}";out.flush();cerr<<"curve "<<k<<" PASS\n";
        }
        out<<"]}\n";cerr<<"seventeen PASS "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';
    }
    fmpq_mpoly_ctx_clear(ctx);flint_cleanup();
}
