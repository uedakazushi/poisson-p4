// NEW continuation: full forecast of at most16 cheap-fill candidates, not global.
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
struct Height{long N,D;};
struct Stats{long zbits=0,bytes=0;};
using StatMatrix=vector<vector<Stats>>;
StatMatrix statistics(const vector<vector<P>>&A){
    StatMatrix s(236,vector<Stats>(83));fmpz_t c;fmpz_init(c);
    for(int i=0;i<236;i++)for(int j=0;j<83;j++)if(A[i][j].n){
        for(long k=0;k<A[i][j].n;k++){fmpz_mpoly_get_term_coeff_fmpz(c,A[i][j].v->zpoly,k,ctx->zctx);s[i][j].zbits=max(s[i][j].zbits,long(fmpz_bits(c)));}
        ostringstream text;put(text,A[i][j]);s[i][j].bytes=text.str().size();
    }
    fmpz_clear(c);return s;
}
long mpzdigits(const fmpz_t a){return fmpz_is_zero(a)?1:fmpz_sizeinbase(a,10);}
long shiftbytes(const Exp&e){long b=2+4;for(long x:e)b+=to_string(x).size();return b;}
long scaledbytes(const P&a,const Stats&s,const fmpq_t unit,const Exp&pe){
    fmpq_t c;fmpq_init(c);fmpq_div(c,a.v->content,unit);Exp sh=a.shift;for(int v=0;v<5;v++)sh[v]-=pe[v];
    long b=s.bytes-qs(a.v->content).size()+qs(c).size()-shiftbytes(a.shift)+shiftbytes(sh);fmpq_clear(c);return b;
}
long tightbytes(const set<Exp>&u,const Exp&base,const fmpz_t M,const fmpz_t D){
    if(u.empty()||fmpz_is_zero(M))return 0;
    ostringstream h;h<<"{\"content\":\"\",\"shift\":[";for(int v=0;v<5;v++){if(v)h<<',';h<<base[v];}h<<"],\"terms\":[]}";
    long ebytes=0;for(const auto&e:u)for(int v=0;v<5;v++){assert(e[v]>=base[v]);ebytes+=to_string(e[v]-base[v]).size();}
    // If H_i=g Z_i, |H_i|<=M, then digits(g)+sum digits(Z_i)
    // <=n*digits(M)+1. Denominator is stored once, exponents exactly.
    return h.str().size()+long(u.size())*(mpzdigits(M)+11)+mpzdigits(D)+16+ebytes;
}
long numerator_bound_bytes(const P*old,long oldbits,const P&factor,long fbits,const P&pivot,long pbits,const fmpq_t unit,const set<Exp>&u,const Exp&base,Height&height){
    fmpq_t cp;fmpq_init(cp);fmpq_div(cp,pivot.v->content,unit);fmpq_mul(cp,cp,factor.v->content);
    fmpz_t D,ma,mf,mp,a,b,t,M,da;fmpz_init(D);fmpz_init(ma);fmpz_init(mf);fmpz_init(mp);fmpz_init(a);fmpz_init(b);fmpz_init(t);fmpz_init(M);fmpz_init(da);
    if(old&&old->n)fmpz_set(da,fmpq_denref(old->v->content));else fmpz_one(da);
    fmpz_lcm(D,da,fmpq_denref(cp));
    fmpz_one(ma);fmpz_mul_2exp(ma,ma,oldbits);fmpz_sub_ui(ma,ma,1);
    fmpz_one(mf);fmpz_mul_2exp(mf,mf,fbits);fmpz_sub_ui(mf,mf,1);
    fmpz_one(mp);fmpz_mul_2exp(mp,mp,pbits);fmpz_sub_ui(mp,mp,1);
    if(old&&old->n){fmpz_divexact(a,D,da);fmpz_mul(a,a,fmpq_numref(old->v->content));fmpz_abs(a,a);fmpz_mul(a,a,ma);}else fmpz_zero(a);
    fmpz_divexact(b,D,fmpq_denref(cp));fmpz_mul(b,b,fmpq_numref(cp));fmpz_abs(b,b);
    fmpz_mul(b,b,mf);fmpz_mul(b,b,mp);fmpz_mul_ui(b,b,min(factor.n,pivot.n));fmpz_add(M,a,b);
    height={long(fmpz_bits(M)),long(fmpz_bits(D))};long answer=tightbytes(u,base,M,D);
    fmpq_clear(cp);fmpz_clear(D);fmpz_clear(ma);fmpz_clear(mf);fmpz_clear(mp);fmpz_clear(a);fmpz_clear(b);fmpz_clear(t);fmpz_clear(M);fmpz_clear(da);return answer;
}
long bytes_on_disk(const string&dir){long n=0;for(auto&p:filesystem::directory_iterator(dir))if(p.is_regular_file())n+=p.file_size();return n;}
using Matrix=vector<vector<P>>;
struct Update{int row,col;long oldn,k,b,liveupper,peakbound;Height h;};
struct Forecast{
    bool safe=true;string reason;long peakbound=0,outbound=0,live_end=0,updates=0;
    vector<vector<Update>> byrow=vector<vector<Update>>(236);
};
struct Pick{int r=-1,c=-1,phase=-1;long fill=0,pairs=0,h=0;};
bool checked_exp(__int128 x,long&z){const __int128 limit=(__int128(1)<<60);if(x < -limit||x > limit)return false;z=long(x);return true;}
vector<Pick> candidates(const Matrix&A,const vector<bool>&ar,const vector<bool>&ac){
    vector<bitset<83>> bits(236);vector<long>rn(236),cn(83);
    for(int i=0;i<236;i++)if(ar[i])for(int j=0;j<83;j++)if(ac[j]&&nz(A[i][j])){bits[i].set(j);rn[i]+=A[i][j].n;cn[j]+=A[i][j].n;}
    vector<Pick> answer;
    for(int i=0;i<236;i++)if(ar[i])for(int j=0;j<83;j++)if(ac[j]&&A[i][j].n==1){
        Exp e=exponent(A[i][j],0);int phase=any_of(e.begin(),e.end(),[](long x){return x!=0;})?1:0;
        long fill=0;for(int k=0;k<236;k++)if(k!=i&&ar[k]&&nz(A[k][j]))fill+=(bits[i]&~bits[k]).count();
        long pairs=(rn[i]-1)*(cn[j]-1);answer.push_back({i,j,phase,fill,pairs,0});
    }
    sort(answer.begin(),answer.end(),[](const Pick&a,const Pick&b){return make_tuple(a.fill,a.pairs,a.r,a.c)<make_tuple(b.fill,b.pairs,b.r,b.c);});return answer;
}
Forecast forecast(const Matrix&A,const StatMatrix&stats,const vector<bool>&ar,const vector<bool>&ac,const Pick&p,const string&dir,long incumbent=CAP+1,int bestr=999,int bestc=999,long selectionbytes=0){
    Forecast f;long upper=live;fmpq_t q;fmpq_init(q);fmpq_mpoly_get_term_coeff_fmpq(q,A[p.r][p.c].v,0,ctx);
    Exp pe=exponent(A[p.r][p.c],0);
    long receipt=4096,residual=8192;vector<vector<long>> costs(236,vector<long>(83));
    for(int i=0;i<236;i++)if(ar[i])for(int j=0;j<83;j++)if(ac[j])costs[i][j]=stats[i][j].bytes;
    for(int j=0;j<83;j++)if(ac[j]&&j!=p.c&&nz(A[p.r][j])){
        for(int v=0;v<5;v++){long test;if(!checked_exp(__int128(A[p.r][j].shift[v])-pe[v],test)){f.safe=false;f.reason="exponent_representation_guard";fmpq_clear(q);return f;}}
        receipt+=scaledbytes(A[p.r][j],stats[p.r][j],q,pe)+128;
    }
    for(int i=0;i<236&&f.safe;i++)if(ar[i]&&i!=p.r&&nz(A[i][p.c])){
        const P&factor=A[i][p.c];auto fs=support(factor);receipt+=stats[i][p.c].bytes+512;
        for(int j=0;j<83&&f.safe;j++)if(ac[j]&&j!=p.c&&nz(A[p.r][j])){
            if(seconds()>39){f.safe=false;f.reason="time_forecast_before_pivot";break;}
            auto ps=support(A[p.r][j]);long a=A[i][j].n;auto old=support(A[i][j]);
            set<Exp> product_support,united(old.begin(),old.end());
            for(auto x:fs)for(auto y:ps){
                Exp e;for(int v=0;v<5;v++)if(!checked_exp(__int128(x[v])+y[v]-pe[v],e[v])){f.safe=false;f.reason="exponent_representation_guard";}
                if(!f.safe)break;product_support.insert(e);united.insert(e);
                long support_scratch=fs.size()+ps.size()+old.size()+product_support.size()+united.size();
                f.peakbound=max(f.peakbound,live+support_scratch+64);
                if(f.peakbound>CAP){f.safe=false;f.reason="support_forecast_scratch_cap";break;}
            }
            if(!f.safe)break;
            long k=product_support.size(),b=united.size();long bound=upper+2*k+3*b+64;
            f.peakbound=max(f.peakbound,bound);
            if(bound>CAP){f.safe=false;f.reason="next_row_update_live_coefficient_forecast_cap";break;}
            if(f.peakbound>incumbent||(f.peakbound==incumbent&&make_pair(p.r,p.c)>=make_pair(bestr,bestc))){f.safe=false;f.reason="dominated_by_proved_incumbent_forecast";break;}
            Exp base;for(int v=0;v<5;v++){long prodshift;if(!checked_exp(__int128(factor.shift[v])+A[p.r][j].shift[v]-pe[v],prodshift)){f.safe=false;f.reason="exponent_representation_guard";break;}base[v]=min(A[i][j].shift[v],prodshift);}
            if(!f.safe)break;
            for(const auto&e:united)for(int v=0;v<5;v++)if(__int128(e[v])-base[v] < 0 || __int128(e[v])-base[v] > (__int128(1)<<62)){f.safe=false;f.reason="exponent_representation_guard";}
            if(!f.safe)break;
            Height nh;costs[i][j]=numerator_bound_bytes(a?&A[i][j]:nullptr,stats[i][j].zbits,factor,stats[i][p.c].zbits,A[p.r][j],stats[p.r][j].zbits,q,united,base,nh);
            f.byrow[i].push_back({i,j,a,k,b,upper,bound,nh});f.updates++;receipt+=320;
            upper+=b-a;
        }
        upper-=factor.n;costs[i][p.c]=0;
    }
    for(int j=0;j<83;j++)upper-=A[p.r][j].n;
    for(int i=0;i<236;i++)if(ar[i]&&i!=p.r)for(int j=0;j<83;j++)if(ac[j]&&j!=p.c)residual+=costs[i][j]+32;
    long remaining_rows=count(ar.begin(),ar.end(),true)-1,remaining_cols=count(ac.begin(),ac.end(),true)-1;
    long next_frontier_reserve=128*remaining_rows*remaining_cols;
    f.live_end=upper;f.outbound=raw_input_bytes+bytes_on_disk(dir)+selectionbytes+receipt+residual+1024*1024+next_frontier_reserve;
    if(f.safe&&f.outbound>BYTES){f.safe=false;f.reason="predicted_complete_output_32MiB_cap";}
    fmpq_clear(q);return f;
}
void ints(ostream&out,const vector<int>&v){out<<'[';for(size_t i=0;i<v.size();i++){if(i)out<<',';out<<v[i];}out<<']';}
int main(int argc,char**argv){
    assert(argc==3);epoch=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,5,ORD_DEGLEX);int exitcode=0;
    {
        try{
            string dir=argv[2];raw_input_bytes=filesystem::file_size(argv[1]);ifstream in(argv[1]);int v,nr,nc,nar,nac,nn;in>>v>>nr>>nc>>nar>>nac>>nn;assert(v==5&&nr==236&&nc==83&&nar==206&&nac==53&&nn==6078);
            vector<bool>ar(236,false),ac(83,false),unitvars(5,false);for(int k=0;k<nar;k++){int i;in>>i;ar[i]=true;}for(int k=0;k<nac;k++){int j;in>>j;ac[j]=true;}
            Matrix A(236,vector<P>(83));
            for(int k=0;k<nn;k++){int i,j;long n;in>>i>>j>>n;Exp sh;for(auto&e:sh)in>>e;assert(ar[i]&&ac[j]&&!nz(A[i][j]));A[i][j]=read(in,n);A[i][j].shift=sh;}
            string tail;assert(in&&!(in>>tail));assert(live==85459);
            vector<int>prow,pcol;int constants=0,monomials=0;long identities=0;
            string stop;Pick stopped;Forecast last;
            while(true){
                if(seconds()>39){stop="time_before_next_pivot";break;}
                StatMatrix stats=statistics(A);Pick p;Forecast f;
                vector<tuple<int,int,string,long,long>> selection_records;
                string selection_mode;
                if(prow.empty()){
                    auto first_options=candidates(A,ar,ac);p=*find_if(first_options.begin(),first_options.end(),[](const Pick&x){return x.r==12&&x.c==67;});assert(A[12][67].n==1);
                    assert(exponent(A[12][67],0)==Exp({0,0,0,0,1}));
                    fmpq_t test,want;fmpq_init(test);fmpq_init(want);fmpq_set_si(want,-32,405);fmpq_mpoly_get_term_coeff_fmpq(test,A[12][67].v,0,ctx);assert(fmpq_equal(test,want));fmpq_clear(test);fmpq_clear(want);
                    f=forecast(A,stats,ar,ac,p,dir,CAP+1,999,999,4096);
                    selection_mode="best_FULL_safe_saved_examined_record_NOT_global_minimum";
                    selection_records.emplace_back(p.r,p.c,f.safe?"safe":"unsafe",f.peakbound,f.outbound);
                    if(!f.safe){stop=f.reason;stopped=p;last=move(f);break;}
                }else{
                    auto options=candidates(A,ar,ac);long all_eligible=options.size();if(options.size()>16)options.resize(16);
                    if(options.empty()){stop="no_QQ_constant_or_parameter_monomial_entry";break;}
                    long best=CAP+1;bool found=false,incomplete=false;long selection_bytes=4096+128*options.size();
                    for(const auto&candidate:options){
                        if(seconds()>39){incomplete=true;break;}
                        Forecast g=forecast(A,stats,ar,ac,candidate,dir,best,p.r<0?999:p.r,p.c<0?999:p.c,selection_bytes);
                        string status=g.safe?"safe":g.reason=="dominated_by_proved_incumbent_forecast"?"dominated":g.reason=="predicted_complete_output_32MiB_cap"?"bytes":g.reason=="time_forecast_before_pivot"?"time":"terms";
                        selection_records.emplace_back(candidate.r,candidate.c,status,g.peakbound,g.outbound);
                        if(status=="time"){incomplete=true;break;}
                        if(g.safe&&(!found||make_tuple(g.peakbound,candidate.r,candidate.c)<make_tuple(best,p.r,p.c))){best=g.peakbound;p=candidate;f=move(g);found=true;}
                    }
                    selection_mode="best_safe_among_at_most16_cheapfill_shortlisted_NOT_global";
                    if(incomplete||!found){
                        stop=incomplete?"shortlist_selection_not_completed_before_time_guard":"no_safe_pivot_among_at_most16_examined_NOT_all_unsafe";
                        ofstream front(dir+"/selection-frontier.json");front<<"{\"no_pivot_executed\":true,\"eligible_candidates\":"<<all_eligible<<",\"shortlist_candidates\":"<<options.size()<<",\"records\":[";
                        bool co=false;for(auto&rec:selection_records){if(co)front<<',';co=true;front<<'['<<get<0>(rec)<<','<<get<1>(rec)<<",\""<<get<2>(rec)<<"\","<<get<3>(rec)<<','<<get<4>(rec)<<']';}front<<"]}\n";break;
                    }
                }
                maxforecast=max(maxforecast,f.peakbound);
                fmpq_t q;fmpq_init(q);fmpq_mpoly_get_term_coeff_fmpq(q,A[p.r][p.c].v,0,ctx);Exp pe=exponent(A[p.r][p.c],0);
                for(int j=0;j<5;j++)if(pe[j])unitvars[j]=true;
                int step=prow.size();string number=to_string(step);if(step<10)number="0"+number;
                {
                    ofstream sel(dir+"/selection-"+number+".json");sel<<"{\"mode\":\""<<selection_mode<<"\",\"selected\":["<<p.r<<','<<p.c<<"],\"not_global_minimum_claim\":true,\"all_shortlisted_considered_or_dominated_before_arithmetic\":"<<(step?"true":"false")<<",\"records\":[";
                    bool co=false;for(auto&rec:selection_records){if(co)sel<<',';co=true;sel<<'['<<get<0>(rec)<<','<<get<1>(rec)<<",\""<<get<2>(rec)<<"\","<<get<3>(rec)<<','<<get<4>(rec)<<']';}sel<<"]}\n";sel.close();assert(sel.good());
                }
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
                fmpq_clear(q);assert(raw_input_bytes+bytes_on_disk(dir)<BYTES);
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
            summary<<"],\"literal_row_update_identity_checks\":"<<identities<<",\"peak_live_terms\":"<<peak<<",\"max_forecast_terms\":"<<maxforecast<<",\"source_initial_terms\":85459,\"prior30_steps_replayed\":0,\"shortlist_limit\":16,\"not_global_minimum_claim\":true,\"numeric_bound_scratch_reserve\":64,\"raw_input_bytes_included_in32MiB\":"<<raw_input_bytes<<",\"future_frontier_RAW_reserve_used\":true,\"native_soft_seconds\":39,\"native_seconds\":"<<seconds()<<"}\n";summary.close();assert(summary.good());assert(raw_input_bytes+bytes_on_disk(dir)<BYTES);
        }catch(const exception&e){cerr<<"STOP_EXCEPTION "<<e.what()<<'\n';exitcode=2;}
    }
    assert(live==0);fmpq_mpoly_ctx_clear(ctx);flint_cleanup();return exitcode;
}
