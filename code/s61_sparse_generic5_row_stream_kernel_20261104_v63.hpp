// NEW continuation from36: first FULL-safe frontier candidate, then shortlist16.
#include <flint/fmpq_mpoly.h>
#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
using namespace std;
using Exp=array<long,5>;
fmpq_mpoly_ctx_t ctx;
long live=0,peak=0,maxforecast=0;
long raw_input_bytes=0;
const long CAP=150000,BYTES=32*1024*1024;
auto epoch=chrono::steady_clock::now();
double seconds(){return chrono::duration<double>(chrono::steady_clock::now()-epoch).count();}
void guard(){peak=max(peak,live);assert(live<=CAP);}
struct P{
    fmpq_mpoly_t v;Exp shift{};long n=0;
    P(){fmpq_mpoly_init(v,ctx);}
    void sync(){long k=fmpq_mpoly_length(v,ctx);live+=k-n;n=k;guard();}
    P(const P&a):P(){shift=a.shift;fmpq_mpoly_set(v,a.v,ctx);sync();}
    P(P&&a):P(){fmpq_mpoly_swap(v,a.v,ctx);shift=a.shift;n=a.n;a.n=0;}
    P&operator=(P a){fmpq_mpoly_swap(v,a.v,ctx);swap(shift,a.shift);swap(n,a.n);return *this;}
    ~P(){live-=n;fmpq_mpoly_clear(v,ctx);}
};
bool nz(const P&a){return a.n!=0;}
Exp exponent(const P&a,long i){
    ulong e[5];fmpq_mpoly_get_term_exp_ui(e,a.v,i,ctx);Exp z;
    for(int j=0;j<5;j++){assert(e[j]<=ulong(numeric_limits<long>::max()));z[j]=a.shift[j]+long(e[j]);}return z;
}
vector<Exp> support(const P&a){vector<Exp>s;for(long i=0;i<a.n;i++)s.push_back(exponent(a,i));return s;}
P aligned(const P&a,const Exp&base){
    P z(a);ulong e[5];
    for(long i=0;i<z.n;i++){
        fmpq_mpoly_get_term_exp_ui(e,z.v,i,ctx);
        for(int j=0;j<5;j++){assert(a.shift[j]>=base[j]);e[j]+=ulong(a.shift[j]-base[j]);}
        fmpq_mpoly_set_term_exp_ui(z.v,i,e,ctx);
    }
    z.shift=base;return z;
}
P add(const P&a,const P&b,bool minus=false){
    Exp base;for(int j=0;j<5;j++)base[j]=min(a.shift[j],b.shift[j]);
    P aa=aligned(a,base),bb=aligned(b,base),z;z.shift=base;
    if(minus)fmpq_mpoly_sub(z.v,aa.v,bb.v,ctx);else fmpq_mpoly_add(z.v,aa.v,bb.v,ctx);
    z.sync();return z;
}
P mul(const P&a,const P&b){P z;for(int j=0;j<5;j++)z.shift[j]=a.shift[j]+b.shift[j];fmpq_mpoly_mul(z.v,a.v,b.v,ctx);z.sync();return z;}
bool equal(const P&a,const P&b){
    if(a.n!=b.n)return false;fmpq_t x,y;fmpq_init(x);fmpq_init(y);bool ok=true;
    for(long i=0;i<a.n&&ok;i++){
        ok=exponent(a,i)==exponent(b,i);
        fmpq_mpoly_get_term_coeff_fmpq(x,a.v,i,ctx);fmpq_mpoly_get_term_coeff_fmpq(y,b.v,i,ctx);
        ok=ok&&fmpq_equal(x,y);
    }
    fmpq_clear(x);fmpq_clear(y);return ok;
}
string qs(const fmpq_t q){char*c=fmpq_get_str(nullptr,10,q);string s(c);flint_free(c);return s;}
string zs(const fmpz_t q){char*c=fmpz_get_str(nullptr,10,q);string s(c);flint_free(c);return s;}
void put(ostream&out,const P&a){
    out<<"{\"content\":\""<<qs(a.v->content)<<"\",\"shift\":[";
    for(int j=0;j<5;j++){if(j)out<<',';out<<a.shift[j];}out<<"],\"terms\":[";
    fmpz_t c;fmpz_init(c);ulong e[5];
    for(long i=0;i<a.n;i++){
        if(i)out<<',';fmpz_mpoly_get_term_coeff_fmpz(c,a.v->zpoly,i,ctx->zctx);
        fmpq_mpoly_get_term_exp_ui(e,a.v,i,ctx);out<<"[\""<<zs(c)<<"\"";for(auto z:e)out<<','<<z;out<<']';
    }
    fmpz_clear(c);out<<"]}";
}
P read(istream&in,long n){
    P z;fmpq_t q;fmpq_init(q);ulong e[5];
    for(long k=0;k<n;k++){string s;in>>s;for(auto&x:e)in>>x;assert(fmpq_set_str(q,s.c_str(),10)==0);fmpq_canonicalise(q);fmpq_mpoly_set_coeff_fmpq_ui(z.v,q,e,ctx);}
    fmpq_clear(q);z.sync();assert(z.n==n);return z;
}
