// One exact QQ[H,L,s0,q0,t0] constrained83 source block, no field solve.
#include <flint/fmpq_mpoly.h>
#include <flint/fmpz_mpoly.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
fmpq_mpoly_ctx_t cx;
long live_terms=0,peak_live_terms=0,max_single_terms=0,operations=0;
void guard(){
    peak_live_terms=max(peak_live_terms,live_terms);
    if(live_terms>150000) throw runtime_error("150000 live expanded terms exceeded");
}
struct P {
    fmpq_mpoly_t p;
    long terms=0;
    P(){fmpq_mpoly_init(p,cx);}
    void sync(){long n=fmpq_mpoly_length(p,cx);live_terms+=n-terms;terms=n;max_single_terms=max(max_single_terms,n);guard();}
    P(const P& a):P(){fmpq_mpoly_set(p,a.p,cx);sync();}
    P(P&& a):P(){fmpq_mpoly_swap(p,a.p,cx);terms=a.terms;a.terms=0;}
    P& operator=(P a){fmpq_mpoly_swap(p,a.p,cx);swap(terms,a.terms);return *this;}
    ~P(){live_terms-=terms;fmpq_mpoly_clear(p,cx);}
};
bool zero(const P&a){return fmpq_mpoly_is_zero(a.p,cx);}
bool same(const P&a,const P&b){return fmpq_mpoly_equal(a.p,b.p,cx);}
P sum(const P&a,const P&b){P z;fmpq_mpoly_add(z.p,a.p,b.p,cx);z.sync();operations++;return z;}
P product(const P&a,const P&b){P z;fmpq_mpoly_mul(z.p,a.p,b.p,cx);z.sync();operations++;return z;}
void accumulate(P& a,const P&b){if(!zero(b))a=sum(a,b);}
P read_poly(istream&in){
    P z;long n;in>>n;assert(n>=0);fmpq_t c;fmpq_init(c);ulong e[5];
    for(long j=0;j<n;j++){
        string s;in>>s;for(auto&v:e)in>>v;
        assert(fmpq_set_str(c,s.c_str(),10)==0);fmpq_canonicalise(c);
        fmpq_mpoly_set_coeff_fmpq_ui(z.p,c,e,cx);
    }
    fmpq_clear(c);z.sync();assert(z.terms==n);return z;
}
string qtext(const fmpq_t q){char*s=fmpq_get_str(nullptr,10,q);string z(s);flint_free(s);return z;}
string ztext(const fmpz_t q){char*s=fmpz_get_str(nullptr,10,q);string z(s);flint_free(s);return z;}
void write_poly(ostream&out,const P& a){
    out<<"{\"content\":\""<<qtext(a.p->content)<<"\",\"terms\":[";
    ulong e[5];fmpz_t c;fmpz_init(c);
    for(long j=0;j<a.terms;j++){
        if(j)out<<',';
        fmpz_mpoly_get_term_coeff_fmpz(c,a.p->zpoly,j,cx->zctx);
        fmpz_mpoly_get_term_exp_ui(e,a.p->zpoly,j,cx->zctx);
        out<<"[\""<<ztext(c)<<"\"";for(auto v:e)out<<','<<v;out<<']';
    }
    fmpz_clear(c);out<<"]}";
}
struct Term {vector<int> q;P coefficient;};
using Row=vector<P>;
using Matrix=vector<Row>;
int main(int argc,char**argv){
    assert(argc==3);auto started=chrono::steady_clock::now();
    fmpq_mpoly_ctx_init(cx,5,ORD_DEGLEX);
    int status=0;
    {
        try{
            ifstream in(argv[1]);int nv,nr,nq,nn,nl;in>>nv>>nr>>nq>>nn>>nl;
            assert(nv==5&&nr==134&&nq==45&&nn==20&&nl==51);
            Row center(45),d(45);
            for(auto& p:center)p=read_poly(in);for(auto& p:d)p=read_poly(in);
            Matrix N(45,Row(20)),L(51,Row(20));
            for(auto&row:N)for(auto&p:row)p=read_poly(in);
            for(auto&row:L)for(auto&p:row)p=read_poly(in);
            vector<vector<Term>> equations(134);
            for(auto&eq:equations){
                int k;in>>k;assert(k>=0);
                for(int j=0;j<k;j++){
                    Term t;int n;in>>n;assert(0<=n&&n<=2);t.q.resize(n);
                    for(auto&i:t.q){in>>i;assert(0<=i&&i<45);}
                    t.coefficient=read_poly(in);eq.push_back(move(t));
                }
            }
            string trailing;assert(in&&!(in>>trailing));
            const long input_terms=live_terms;
            vector<int> normal;for(int i=0;i<45;i++)if(i!=3&&i!=6&&i!=13)normal.push_back(i);
            Matrix top(134,Row(83));
            long fc_terms=0,jn_terms=0,dn_terms=0,raw_terms=0,equalities=0;
            long j_nonzero=0,d_nonzero=0,raw_nonzero=0,fc_nonzero=0;
            for(int i=0;i<134;i++){
                Row jacobian(45),polar(45);
                P fc;
                for(const auto& t:equations[i]){
                    const P& p=t.coefficient;
                    if(t.q.empty())accumulate(fc,p);
                    else if(t.q.size()==1){
                        int a=t.q[0];accumulate(fc,product(p,center[a]));accumulate(jacobian[a],p);
                    }else{
                        int a=t.q[0],b=t.q[1];
                        accumulate(fc,product(p,product(center[a],center[b])));
                        accumulate(jacobian[a],product(p,center[b]));
                        accumulate(jacobian[b],product(p,center[a]));
                        accumulate(polar[a],product(p,d[b]));
                        accumulate(polar[b],product(p,d[a]));
                    }
                }
                top[i][0]=move(fc);
                for(int k=0;k<20;k++){
                    P jn,dn;
                    for(int a=0;a<45;a++)if(!zero(N[a][k])){
                        accumulate(jn,product(jacobian[a],N[a][k]));
                        accumulate(dn,product(polar[a],N[a][k]));
                    }
                    // Independently assemble the alpha and beta coefficients
                    // of F(C+N*alpha) and F(C+d*sigma+N*beta).
                    P direct_j,direct_d;
                    for(const auto&t:equations[i]){
                        const P&p=t.coefficient;
                        if(t.q.size()==1)accumulate(direct_j,product(p,N[t.q[0]][k]));
                        else if(t.q.size()==2){
                            int a=t.q[0],b=t.q[1];
                            accumulate(direct_j,product(p,sum(product(N[a][k],center[b]),product(N[b][k],center[a]))));
                            accumulate(direct_d,product(p,sum(product(N[a][k],d[b]),product(N[b][k],d[a]))));
                        }
                    }
                    assert(same(jn,direct_j)&&same(dn,direct_d));equalities+=2;
                    top[i][1+k]=move(jn);top[i][21+k]=move(dn);
                }
                for(int k=0;k<42;k++){
                    int a=normal[k];P direct;
                    for(const auto&t:equations[i])if(t.q.size()==2){
                        if(t.q[0]==a)accumulate(direct,product(t.coefficient,d[t.q[1]]));
                        if(t.q[1]==a)accumulate(direct,product(t.coefficient,d[t.q[0]]));
                    }
                    assert(same(direct,polar[a]));equalities++;
                    top[i][41+k]=move(direct);
                }
                fc_terms+=top[i][0].terms;fc_nonzero+=!zero(top[i][0]);
                for(int k=0;k<20;k++){
                    jn_terms+=top[i][1+k].terms;j_nonzero+=!zero(top[i][1+k]);
                    dn_terms+=top[i][21+k].terms;d_nonzero+=!zero(top[i][21+k]);
                }
                for(int k=0;k<42;k++){raw_terms+=top[i][41+k].terms;raw_nonzero+=!zero(top[i][41+k]);}
                cerr<<"ROW_COMPLETE "<<i<<" LIVE_TERMS "<<live_terms<<'\n';
            }
            long lterms=0;for(const auto&row:L)for(const auto&p:row)lterms+=p.terms;
            long nnz=0,total=0;ofstream out(argv[2]);
            out<<"{\"all_passed\":true,\"parameter_variables\":[\"H\",\"L\",\"s0\",\"q0\",\"t0\"],\"shape\":[236,83],\"sparse_entries\":[";
            auto emit=[&](int i,int j,const P&p){
                if(zero(p))return;if(nnz)out<<',';out<<'['<<i<<','<<j<<',';write_poly(out,p);out<<']';nnz++;total+=p.terms;
            };
            for(int i=0;i<134;i++)for(int j=0;j<83;j++)emit(i,j,top[i][j]);
            for(int i=0;i<51;i++)for(int j=0;j<20;j++){emit(134+i,1+j,L[i][j]);emit(185+i,21+j,L[i][j]);}
            out<<"],\"input_terms\":"<<input_terms<<",\"total_matrix_terms\":"<<total<<",\"nonzero_entries\":"<<nnz;
            out<<",\"block_terms\":{\"Fcenter\":"<<fc_terms<<",\"JcenterN\":"<<jn_terms<<",\"polar_d_N\":"<<dn_terms<<",\"raw_polar42\":"<<raw_terms<<",\"two_L\":"<<2*lterms<<"}";
            out<<",\"block_nonzero_entries\":{\"Fcenter\":"<<fc_nonzero<<",\"JcenterN\":"<<j_nonzero<<",\"polar_d_N\":"<<d_nonzero<<",\"raw_polar42\":"<<raw_nonzero<<"}";
            out<<",\"all_QQ_coefficient_identity_checks\":"<<equalities<<",\"peak_live_expanded_terms\":"<<peak_live_terms<<",\"max_single_polynomial_terms\":"<<max_single_terms<<",\"polynomial_operations\":"<<operations;
            out<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"}\n";
            out.close();assert(out.good());
            ifstream sizein(argv[2],ios::ate|ios::binary);assert(sizein.tellg()<32*1024*1024);
            assert(input_terms==29028&&total<=53315&&equalities==134*(40+42));
        }catch(const exception& e){cerr<<"BOUNDED_STOP "<<e.what()<<'\n';status=2;}
    }
    assert(live_terms==0);fmpq_mpoly_ctx_clear(cx);flint_cleanup();return status;
}
