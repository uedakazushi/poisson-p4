// Deterministic QQ-constant / parameter-monomial unit prefix, no field GCD.
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
#include <string>
#include <tuple>
#include <vector>
using namespace std;
using Exp=array<long,5>;
fmpq_mpoly_ctx_t ctx;
long live=0,peak=0,maxforecast=0;
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
struct Height{long N,D;};
Height height(const P&a){
    if(!a.n)return {0,0};long b=0;fmpz_t c;fmpz_init(c);
    for(long i=0;i<a.n;i++){fmpz_mpoly_get_term_coeff_fmpz(c,a.v->zpoly,i,ctx->zctx);b=max(b,long(fmpz_bits(c)));}
    fmpz_clear(c);return {b+long(fmpz_bits(fmpq_numref(a.v->content))),long(fmpz_bits(fmpq_denref(a.v->content)))};
}
long digits(long bits){assert(bits>=0&&bits<100000000);return (bits*30103+99999)/100000+2;}
long cellbytes(long n,Height h){return n?160+n*(digits(h.N)+digits(h.D)+5*22+48):0;}
long bytes_on_disk(const string&dir){long n=0;for(auto&p:filesystem::directory_iterator(dir))if(p.is_regular_file())n+=p.file_size();return n;}
using Matrix=vector<vector<P>>;
struct Update{int row,col;long oldn,k,b,liveupper,peakbound;Height h;};
struct Forecast{
    bool safe=true;string reason;long peakbound=0,outbound=0,live_end=0,updates=0;
    vector<vector<Update>> byrow=vector<vector<Update>>(236);
};
struct Pick{int r=-1,c=-1,phase=-1;long fill=0,pairs=0,h=0;};
Pick choose(const Matrix&A,const vector<bool>&ar,const vector<bool>&ac){
    vector<bitset<83>> bits(236);vector<long>rn(236),cn(83);
    for(int i=0;i<236;i++)if(ar[i])for(int j=0;j<83;j++)if(ac[j]&&nz(A[i][j])){bits[i].set(j);rn[i]+=A[i][j].n;cn[j]+=A[i][j].n;}
    Pick best;tuple<int,long,long,long,int,int> score;bool found=false;
    for(int i=0;i<236;i++)if(ar[i])for(int j=0;j<83;j++)if(ac[j]&&A[i][j].n==1){
        Exp e=exponent(A[i][j],0);int phase=any_of(e.begin(),e.end(),[](long x){return x!=0;})?1:0;
        long fill=0;for(int k=0;k<236;k++)if(k!=i&&ar[k]&&nz(A[k][j]))fill+=(bits[i]&~bits[k]).count();
        long pairs=(rn[i]-1)*(cn[j]-1);Height hh=height(A[i][j]);long h=hh.N+hh.D;
        auto s=make_tuple(phase,fill,pairs,h,i,j);
        if(!found||s<score){found=true;score=s;best={i,j,phase,fill,pairs,h};}
    }
    return best;
}
Forecast forecast(const Matrix&A,const vector<bool>&ar,const vector<bool>&ac,const Pick&p,const string&dir){
    Forecast f;long upper=live;fmpq_t q;fmpq_init(q);fmpq_mpoly_get_term_coeff_fmpq(q,A[p.r][p.c].v,0,ctx);
    Height unit={long(fmpz_bits(fmpq_numref(q))),long(fmpz_bits(fmpq_denref(q)))};
    Exp pe=exponent(A[p.r][p.c],0);
    long receipt=4096,residual=8192;vector<vector<long>> costs(236,vector<long>(83));
    for(int i=0;i<236;i++)if(ar[i])for(int j=0;j<83;j++)if(ac[j])costs[i][j]=cellbytes(A[i][j].n,height(A[i][j]));
    for(int j=0;j<83;j++)if(ac[j]&&j!=p.c&&nz(A[p.r][j])){
        auto h=height(A[p.r][j]);h.N+=unit.D;h.D+=unit.N;receipt+=cellbytes(A[p.r][j].n,h)+128;
    }
    for(int i=0;i<236&&f.safe;i++)if(ar[i]&&i!=p.r&&nz(A[i][p.c])){
        const P&factor=A[i][p.c];auto fs=support(factor);Height fh=height(factor);
        receipt+=cellbytes(factor.n,fh)+512;
        for(int j=0;j<83&&f.safe;j++)if(ac[j]&&j!=p.c&&nz(A[p.r][j])){
            if(seconds()>47){f.safe=false;f.reason="time_forecast_before_pivot";break;}
            auto ps=support(A[p.r][j]);long a=A[i][j].n;auto old=support(A[i][j]);
            set<Exp> product_support,united(old.begin(),old.end());
            for(auto x:fs)for(auto y:ps){
                Exp e;for(int v=0;v<5;v++)e[v]=x[v]+y[v]-pe[v];product_support.insert(e);united.insert(e);
                long support_scratch=fs.size()+ps.size()+old.size()+product_support.size()+united.size();
                f.peakbound=max(f.peakbound,live+support_scratch+16);
                if(f.peakbound>CAP){f.safe=false;f.reason="support_forecast_scratch_cap";break;}
            }
            if(!f.safe)break;
            long k=product_support.size(),b=united.size();long bound=upper+2*k+3*b+16;
            f.peakbound=max(f.peakbound,bound);
            if(bound>CAP){f.safe=false;f.reason="next_row_update_live_coefficient_forecast_cap";break;}
            auto ph=height(A[p.r][j]);ph.N+=unit.D;ph.D+=unit.N;
            long convolutions=min(factor.n,A[p.r][j].n),lg=0;while((1L<<lg)<convolutions)lg++;
            Height mh={fh.N+ph.N+lg+1,fh.D+ph.D},nh=mh;
            if(a){auto ah=height(A[i][j]);nh={max(ah.N+mh.D,mh.N+ah.D)+1,ah.D+mh.D};}
            costs[i][j]=cellbytes(b,nh);
            f.byrow[i].push_back({i,j,a,k,b,upper,bound,nh});f.updates++;receipt+=320;
            upper+=b-a;
        }
        upper-=factor.n;costs[i][p.c]=0;
    }
    for(int j=0;j<83;j++)upper-=A[p.r][j].n;
    for(int i=0;i<236;i++)if(ar[i]&&i!=p.r)for(int j=0;j<83;j++)if(ac[j]&&j!=p.c)residual+=costs[i][j]+32;
    f.live_end=upper;f.outbound=bytes_on_disk(dir)+receipt+residual+1024*1024;
    if(f.outbound>BYTES){f.safe=false;f.reason="predicted_complete_output_32MiB_cap";}
    fmpq_clear(q);return f;
}
void ints(ostream&out,const vector<int>&v){out<<'[';for(size_t i=0;i<v.size();i++){if(i)out<<',';out<<v[i];}out<<']';}
int main(int argc,char**argv){
    assert(argc==3);epoch=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,5,ORD_DEGLEX);int exitcode=0;
    {
        try{
            string dir=argv[2];ifstream in(argv[1]);int v,nr,nc,nn;in>>v>>nr>>nc>>nn;assert(v==5&&nr==236&&nc==83&&nn==9431);
            Matrix A(236,vector<P>(83));
            for(int k=0;k<nn;k++){int i,j;long n;in>>i>>j>>n;assert(!nz(A[i][j]));A[i][j]=read(in,n);}
            string tail;assert(in&&!(in>>tail));assert(live==49039);
            vector<bool>ar(236,true),ac(83,true),unitvars(5,false);vector<int>prow,pcol;int constants=0,monomials=0;long identities=0;
            string stop;Pick stopped;Forecast last;
            while(true){
                if(seconds()>47){stop="time_before_next_pivot";break;}
                Pick p=choose(A,ar,ac);if(p.r<0){stop="no_QQ_constant_or_parameter_monomial_entry";break;}
                Forecast f=forecast(A,ar,ac,p,dir);
                if(!f.safe){stop=f.reason;stopped=p;last=move(f);break;}
                maxforecast=max(maxforecast,f.peakbound);
                fmpq_t q;fmpq_init(q);fmpq_mpoly_get_term_coeff_fmpq(q,A[p.r][p.c].v,0,ctx);Exp pe=exponent(A[p.r][p.c],0);
                for(int j=0;j<5;j++)if(pe[j])unitvars[j]=true;
                int step=prow.size();string number=to_string(step);if(step<10)number="0"+number;
                ofstream out(dir+"/pivot-"+number+".json");
                out<<"{\"completed\":true,\"step\":"<<step<<",\"phase\":\""<<(p.phase?"parameter_monomial":"QQ_constant")<<"\",\"pivot_original_row\":"<<p.r<<",\"pivot_original_column\":"<<p.c;
                out<<",\"selection\":{\"structural_fill\":"<<p.fill<<",\"term_pair_majorant\":"<<p.pairs<<",\"height\":"<<p.h<<"},\"pivot_QQ_coefficient\":\""<<qs(q)<<"\",\"pivot_parameter_exponent\":[";
                for(int j=0;j<5;j++){if(j)out<<',';out<<pe[j];}out<<']';
                out<<",\"forecast\":{\"peak_live_bound\":"<<f.peakbound<<",\"output_bytes_upper\":"<<f.outbound<<",\"coefficient_live_end_upper\":"<<f.live_end<<"},\"normalized_pivot_row\":[";
                bool comma=false;
                for(int j=0;j<83;j++)if(ac[j]&&nz(A[p.r][j])){
                    fmpq_mpoly_scalar_div_fmpq(A[p.r][j].v,A[p.r][j].v,q,ctx);
                    for(int x=0;x<5;x++)A[p.r][j].shift[x]-=pe[x];A[p.r][j].sync();
                    if(comma)out<<',';out<<'['<<j<<',';put(out,A[p.r][j]);out<<']';comma=true;
                }
                fmpq_t one;fmpq_init(one);fmpq_mpoly_get_term_coeff_fmpq(one,A[p.r][p.c].v,0,ctx);assert(fmpq_is_one(one)&&exponent(A[p.r][p.c],0)==Exp{});fmpq_clear(one);
                out<<"],\"row_updates\":[";comma=false;
                for(int i=0;i<236;i++)if(ar[i]&&i!=p.r&&nz(A[i][p.c])){
                    if(comma)out<<',';comma=true;out<<"{\"original_row\":"<<i<<",\"factor\":";put(out,A[i][p.c]);out<<",\"entry_forecasts\":[";bool cc=false;
                    for(auto&u:f.byrow[i]){
                        int j=u.col;assert(A[i][j].n==u.oldn);
                        assert(live<=u.liveupper);
                        P prod=mul(A[i][p.c],A[p.r][j]);assert(prod.n<=u.k);
                        P next=add(A[i][j],prod,true);assert(next.n<=u.b);
                        P back=add(next,prod);assert(equal(back,A[i][j]));identities++;
                        if(cc)out<<',';cc=true;
                        out<<'['<<j<<','<<u.oldn<<','<<u.k<<','<<u.b<<','<<u.liveupper<<','<<u.peakbound<<','<<u.h.N<<','<<u.h.D<<']';
                        A[i][j]=move(next);
                    }
                    A[i][p.c]=P();out<<"]}";
                }
                out<<"],\"all_row_update_back_add_identities_checked\":true,\"actual_live_before_retire\":"<<live;
                for(int j=0;j<83;j++)A[p.r][j]=P();ar[p.r]=false;ac[p.c]=false;
                assert(live<=f.live_end);out<<",\"actual_residual_live_terms\":"<<live<<"}\n";out.close();assert(out.good());
                prow.push_back(p.r);pcol.push_back(p.c);if(p.phase)monomials++;else constants++;
                fmpq_clear(q);assert(bytes_on_disk(dir)<BYTES);
                cerr<<"PIVOT "<<step<<" ROW "<<p.r<<" COL "<<p.c<<" LIVE "<<live<<" PEAK "<<peak<<" SECONDS "<<seconds()<<'\n';
            }
            vector<int>rr,cc;for(int i=0;i<236;i++)if(ar[i])rr.push_back(i);for(int j=0;j<83;j++)if(ac[j])cc.push_back(j);
            ofstream residual(dir+"/residual-matrix.json");residual<<"{\"all_passed\":true,\"shape\":["<<rr.size()<<','<<cc.size()<<"],\"row_original_ids\":";ints(residual,rr);residual<<",\"column_original_ids\":";ints(residual,cc);residual<<",\"sparse_entries\":[";
            bool comma=false;long nzcount=0;
            for(int i:rr)for(int j:cc)if(nz(A[i][j])){if(comma)residual<<',';comma=true;residual<<'['<<i<<','<<j<<',';put(residual,A[i][j]);residual<<']';nzcount++;}
            residual<<"],\"total_terms\":"<<live<<",\"nonzero_entries\":"<<nzcount<<"}\n";residual.close();assert(residual.good());
            if(stopped.r>=0){ofstream s(dir+"/stopped-candidate.json");s<<"{\"reason\":\""<<stop<<"\",\"row\":"<<stopped.r<<",\"column\":"<<stopped.c<<",\"forecast_peak\":"<<last.peakbound<<",\"forecast_output_bytes\":"<<last.outbound<<",\"no_part_of_this_pivot_executed\":true}\n";}
            ofstream summary(dir+"/native-summary.json");summary<<"{\"all_passed\":true,\"pivot_count\":"<<prow.size()<<",\"phase_counts\":{\"QQ_constant\":"<<constants<<",\"parameter_monomial\":"<<monomials<<"},\"stop_reason\":\""<<stop<<"\",\"pivot_row_original_order\":";ints(summary,prow);summary<<",\"pivot_column_original_order\":";ints(summary,pcol);
            summary<<",\"unit_variables_inverted\":[";comma=false;vector<string>names={"H","L","s0","q0","t0"};for(int i=0;i<5;i++)if(unitvars[i]){if(comma)summary<<',';comma=true;summary<<'"'<<names[i]<<'"';}
            summary<<"],\"literal_row_update_identity_checks\":"<<identities<<",\"peak_live_terms\":"<<peak<<",\"max_forecast_terms\":"<<maxforecast<<",\"native_seconds\":"<<seconds()<<"}\n";summary.close();assert(summary.good());assert(bytes_on_disk(dir)<BYTES);
        }catch(const exception&e){cerr<<"STOP_EXCEPTION "<<e.what()<<'\n';exitcode=2;}
    }
    assert(live==0);fmpq_mpoly_ctx_clear(ctx);flint_cleanup();return exitcode;
}
