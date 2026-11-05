// ONE fixed two-matrix adjugate worker. No rank/pivot search or field gcd.
#include "s61_sparse_generic5_row_stream_kernel_20261104_v63.hpp"
#include <fstream>
#include <memory>
#include <stdexcept>
const long META=25000,LIMIT=85813,RAWLIMIT=10622833;
namespace fs=std::filesystem;
void need(bool b,const string&s){if(!b)throw runtime_error(s);}
long prospect=0,checks=0,nodes_done=0;fs::path outdir,workdir;long actual_input=0,external=0;
struct Buf:streambuf{
    ofstream f;long n=0,cap;Buf(const fs::path&p,long c):f(p,ios::binary),cap(c){need(bool(f),"write open");}
    streamsize xsputn(const char*s,streamsize k)override{need(n+k<=cap,"file RAW upper");f.write(s,k);n+=k;need(bool(f),"write");return k;}
    int overflow(int c)override{if(c==traits_type::eof())return traits_type::not_eof(c);char a=char(c);xsputn(&a,1);return c;}
    int sync()override{f.flush();return f?0:-1;}
};
struct Output{
    fs::path p,tmp;Buf b;ostream f;Output(const fs::path&q,long cap):p(q),tmp(q.string()+".partial"),b(tmp,cap),f(&b){need(!fs::exists(p),"immutable output");f.exceptions(ios::badbit|ios::failbit);}
    void commit(){f.flush();b.f.close();need(bool(b.f),"close");fs::rename(tmp,p);}
};
long tree(const fs::path&p){long n=0;for(auto&x:fs::recursive_directory_iterator(p))if(x.is_regular_file())n+=x.file_size();return n;}
void raw(){need(actual_input+external+tree(outdir)+tree(workdir)<=RAWLIMIT,"whole RAW preflight guard");}
void owner(){need(live+META<=LIMIT&&peak+META<=LIMIT,"owned incl fixed metadata");}
void forecast(long base,long k,long b){long p=base+2*k+3*b+256+META;prospect=max(prospect,p);need(p<=LIMIT,"whole prospective scratch before operation");}
using Matrix=vector<vector<P>>;
P inputpoly(istream&in){long n;in>>n;P p;fmpq_t q;fmpq_init(q);ulong e[5];for(long k=0;k<n;k++){string s;in>>s;for(auto&v:e)in>>v;need(bool(in)&&fmpq_set_str(q,s.c_str(),10)==0,"literal cleared coefficient");need(fmpz_is_one(fmpq_denref(q)),"cleared integer coefficient");fmpq_mpoly_set_coeff_fmpq_ui(p.v,q,e,ctx);p.sync();owner();}fmpq_clear(q);need(p.n==n,"literal input support");return p;}
void savepoly(const string&name,const P&p,long cap){Output f(outdir/name,cap);put(f.f,p);f.f<<'\n';f.commit();raw();}
P fixedminor(istream&in,const Matrix&A,Output&ev,int block,int output){
    int n;in>>n;vector<int>rr(n),cc(n);for(auto&r:rr)in>>r;for(auto&c:cc)in>>c;int nn;in>>nn;need(nn==(1<<n)-1,"DP count");
    vector<P>d(1<<n);fmpq_mpoly_one(d[0].v,ctx);d[0].sync();owner();
    for(int mask=1;mask<=nn;mask++){
        int m;long count;int nt;in>>m>>count>>nt;need(m==mask&&nt==__builtin_popcount(unsigned(mask)),"fixed mask");P acc;
        for(int h=0;h<nt;h++){
            int c,prev;long na,nb,k,b;in>>c>>prev>>na>>nb>>k>>b;int row=rr[__builtin_popcount(unsigned(mask))-1];
            auto it=find(cc.begin(),cc.end(),c);need(it!=cc.end(),"column in minor");int pos=it-cc.begin();need(prev==(mask^(1<<pos)),"minor source anchor");
            need(A[row][c].n<=na&&d[prev].n<=nb,"source supports");forecast(live-acc.n,k,b);
            P p=mul(A[row][c],d[prev]);need(p.n<=k,"product support");int parity=(__builtin_popcount(unsigned(mask))-1+__builtin_popcount(unsigned(mask&((1<<pos)-1))))&1;
            P next=add(acc,p,parity);need(next.n<=b,"sum support");owner();acc=std::move(next);checks++;
        }
        need(acc.n<=count,"node support");d[mask]=std::move(acc);nodes_done++;ev.f<<'['<<block<<','<<output<<','<<mask<<','<<d[mask].n<<",true]\n";owner();
    }
    return std::move(d.back());
}
int main(int argc,char**argv){
    need(argc==4,"arguments");outdir=argv[2];workdir=argv[3];ifstream in(argv[1]);need(bool(in),"control");in>>actual_input>>external;
    fmpq_mpoly_ctx_init(ctx,5,ORD_DEGLEX);
    {
        vector<Matrix>A(2),B(2);for(int b=0;b<2;b++){A[b].resize(5);B[b].resize(5);for(int i=0;i<5;i++){A[b][i].resize(5);B[b][i].resize(2);for(auto&p:A[b][i])p=inputpoly(in);for(auto&p:B[b][i])p=inputpoly(in);}}
        need(live<=1150,"literal1150 source envelope");long source_terms=live;raw();
        for(int b=0;b<2;b++)for(int i=0;i<5;i++)for(int j=0;j<7;j++)savepoly(string(b?"beta":"alpha")+"-cleared-"+to_string(i)+"-"+to_string(j)+".json",j<5?A[b][i][j]:B[b][i][j-5],20000);
        vector<vector<P>>outputs(2);for(auto&v:outputs)v.resize(26);
        Output events(outdir/"fixed-DP-events.jsonl",200000);
        vector<vector<long>>caps(2,vector<long>(26));
        for(int block=0;block<2;block++)for(int j=0;j<26;j++){
            need(seconds()<45,"soft45s complete-minor boundary");int bi,oi;long cap;in>>bi>>oi>>cap;need(bi==block&&oi==j,"fixed output order");caps[block][j]=cap;
            P value=fixedminor(in,A[block],events,block,j);
            if(j){int r=(j-1)/5,c=(j-1)%5;if((r+c)&1)fmpq_mpoly_neg(value.v,value.v,ctx);}
            need(!j?value.n>0:true,"actual determinant nonzero");outputs[block][j]=std::move(value);owner();
            string label;if(!j)label="det";else{int r=(j-1)/5,c=(j-1)%5;label="adj-"+to_string(c)+"-"+to_string(r);}
            savepoly(string(block?"beta-":"alpha-")+label+".json",outputs[block][j],cap);
            cout<<"block "<<block<<" output "<<j<<" terms "<<outputs[block][j].n<<" live "<<live<<" seconds "<<seconds()<<endl;
        }
        events.commit();raw();
        long left_right=0;
        Output receipts(outdir/"all100-adjugate-identities.jsonl",100000);
        for(int block=0;block<2;block++)for(int order=0;order<2;order++)for(int i=0;i<5;i++)for(int j=0;j<5;j++){
            int bi,ord,ii,jj;in>>bi>>ord>>ii>>jj;need(bi==block&&ord==order&&ii==i&&jj==j,"fixed identity order");P acc;
            for(int k=0;k<5;k++){
                int kk;long na,nb,prod_n,b;in>>kk>>na>>nb>>prod_n>>b;need(kk==k,"identity component");
                const P&a=order?A[block][k][j]:A[block][i][k];const P&c=order?outputs[block][1+k*5+i]:outputs[block][1+j*5+k];
                need(a.n<=na&&c.n<=nb,"fixed identity supports");forecast(live-acc.n,prod_n,b);P p=mul(a,c);need(p.n<=prod_n,"identity product support");P next=add(acc,p);need(next.n<=b,"identity support union");owner();acc=std::move(next);
            }
            need(i==j?equal(acc,outputs[block][0]):!nz(acc),"literal HAdj or AdjH equals determinant I");left_right++;receipts.f<<'['<<block<<','<<order<<','<<i<<','<<j<<",true]\n";
        }
        string extra;need(bool(in)&&!(in>>extra),"full fixed control consumed");receipts.commit();raw();owner();
        long output_terms=live-source_terms;need(output_terms<=43847&&left_right==100,"all outputs and identities");
        Output summary(outdir/"native-summary.json",20000);summary.f<<"{\"all_passed\":true,\"actual_determinants_nonzero\":true,\"all52_outputs\":52,\"all100_left_right_identities\":"<<left_right<<",\"actual_output_terms\":"<<output_terms<<",\"source_terms\":"<<source_terms<<",\"DP_nodes\":"<<nodes_done<<",\"DP_product_adds\":"<<checks<<",\"coefficient_owned_peak\":"<<peak<<",\"fixed_metadata_charge\":"<<META<<",\"whole_owned_upper_checked\":"<<peak+META<<",\"prospective_owned_upper\":"<<prospect<<",\"seconds\":"<<seconds()<<",\"no_rank_search_GCD_point_prime_or_top_substitution\":true}\n";summary.commit();raw();
    }
    need(!live,"release all coefficients");fmpq_mpoly_ctx_clear(ctx);return 0;
}
