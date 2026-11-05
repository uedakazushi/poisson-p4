// ONE16-cell exact coefficient job. No full index, field gcd, point or prime.
#include "s61_sparse_generic5_row_stream_kernel_20261104_v63.hpp"
#include <openssl/evp.h>
#include <iomanip>
#include <map>
#include <memory>
#include <stdexcept>
namespace fs=std::filesystem;
const long META=25000;long LIMIT=0,RAWLIMIT=0,PAIRLIMIT=0;
void need(bool b,const string&s){if(!b)throw runtime_error(s);}
#include "s61_sparse_generic5_primitive_binary_native_20261105_v64.hpp"
fs::path outdir,workdir;long inputbytes=0,external=0,prospect=0,pairs=0,backadds=0,cellcount=0,rowsdone=0;
string matrixpath,matrixsha;
struct Buf:streambuf{
    ofstream f;long n=0,cap;Buf(const fs::path&p,long c):f(p,ios::binary),cap(c){need(bool(f),"output open");}
    streamsize xsputn(const char*s,streamsize k)override{need(n+k<=cap,"per-file RAW guard");f.write(s,k);need(bool(f),"write");n+=k;return k;}
    int overflow(int c)override{if(c==traits_type::eof())return traits_type::not_eof(c);char a=char(c);xsputn(&a,1);return c;}
    int sync()override{f.flush();return f?0:-1;}
};
struct Output{
    fs::path path,tmp;Buf b;ostream f;Output(const fs::path&p,long cap):path(p),tmp(p.string()+".partial"),b(tmp,cap),f(&b){need(!fs::exists(path),"immutable output");f.exceptions(ios::badbit|ios::failbit);}
    void commit(){f.flush();b.f.close();need(bool(b.f),"close");fs::rename(tmp,path);}
};
string digest(const string&s){unsigned char md[EVP_MAX_MD_SIZE];unsigned n=0;need(EVP_Digest(s.data(),s.size(),md,&n,EVP_sha256(),nullptr)==1,"SHA");ostringstream o;for(unsigned i=0;i<n;i++)o<<hex<<setw(2)<<setfill('0')<<int(md[i]);return o.str();}
string filehash(const fs::path&p){ifstream f(p,ios::binary);need(bool(f),"hash open");unique_ptr<EVP_MD_CTX,decltype(&EVP_MD_CTX_free)>c(EVP_MD_CTX_new(),EVP_MD_CTX_free);need(c&&EVP_DigestInit_ex(c.get(),EVP_sha256(),nullptr)==1,"hash init");char b[65536];while(f){f.read(b,sizeof(b));if(f.gcount())need(EVP_DigestUpdate(c.get(),b,f.gcount())==1,"hash update");}unsigned char md[EVP_MAX_MD_SIZE];unsigned n;need(EVP_DigestFinal_ex(c.get(),md,&n)==1,"hash final");ostringstream o;for(unsigned i=0;i<n;i++)o<<hex<<setw(2)<<setfill('0')<<int(md[i]);return o.str();}
// Canonical actual-exponent QQ coefficient hash; unaffected by internal shifts.
string valuehash(const P&p){
    unique_ptr<EVP_MD_CTX,decltype(&EVP_MD_CTX_free)>c(EVP_MD_CTX_new(),EVP_MD_CTX_free);need(c&&EVP_DigestInit_ex(c.get(),EVP_sha256(),nullptr)==1,"value hash init");
    auto emit=[&](const string&s){need(EVP_DigestUpdate(c.get(),s.data(),s.size())==1,"value hash update");};emit(to_string(p.n)+"\n");fmpq_t q;fmpq_init(q);
    for(long i=0;i<p.n;i++){fmpq_mpoly_get_term_coeff_fmpq(q,p.v,i,ctx);emit(qs(q));auto e=exponent(p,i);for(auto v:e)emit(" "+to_string(v));emit("\n");}
    fmpq_clear(q);unsigned char md[EVP_MAX_MD_SIZE];unsigned n;need(EVP_DigestFinal_ex(c.get(),md,&n)==1,"value hash final");ostringstream o;for(unsigned i=0;i<n;i++)o<<hex<<setw(2)<<setfill('0')<<int(md[i]);return o.str();
}
long tree(const fs::path&p){long n=0;for(auto&e:fs::recursive_directory_iterator(p))if(e.is_regular_file())n+=e.file_size();return n;}
void raw(){need(inputbytes+external+tree(outdir)+tree(workdir)<=RAWLIMIT,"whole RAW guard");}
void owner(){need(live+META<=LIMIT&&peak+META<=LIMIT,"owned+thin metadata guard");}
void forecast(long base,long k,long b,long bound){long p=base+2*k+3*b+256+META;prospect=max(prospect,p);need(p<=bound&&p<=LIMIT,"prospective source bound before multiplication");}
P literal(istream&in){P p;string scalar,h;in>>scalar;for(auto&v:p.shift)in>>v;long n;in>>n>>h;fmpq_t q,c;fmpq_init(q);fmpq_init(c);need(fmpq_set_str(c,scalar.c_str(),10)==0,"literal scalar");ulong e[5];
    for(long k=0;k<n;k++){string z;in>>z;for(auto&v:e)in>>v;need(fmpq_set_str(q,z.c_str(),10)==0,"literal coefficient");fmpq_mul(q,q,c);fmpq_mpoly_set_coeff_fmpq_ui(p.v,q,e,ctx);p.sync();owner();}fmpq_clear(q);fmpq_clear(c);need(p.n==n&&valuehash(p)==h,"actual cached numerator/P scalar-monomial source binding");return p;}
struct Ref{int col;long off,bytes,n;string sha;};
int main(int argc,char**argv){need(argc==4,"arguments");outdir=argv[2];workdir=argv[3];ifstream in(argv[1]);need(bool(in),"control");in>>inputbytes>>external>>LIMIT>>RAWLIMIT>>PAIRLIMIT>>matrixpath>>matrixsha;need(LIMIT<=150000&&RAWLIMIT<=32*1024*1024&&LIMIT>META,"authenticated source caps");fmpq_mpoly_ctx_init(ctx,5,ORD_DEGLEX);
    {
        vector<P>global;global.reserve(15);for(int i=0;i<15;i++)global.push_back(literal(in));need(live==4603,"4313 cached numerators+290 P");need(filehash(matrixpath)==matrixsha,"literal current matrix SHA");ifstream src(matrixpath,ios::binary);string magic(8,0);src.read(magic.data(),8);need(magic=="Q5PRIM1\n","source binary magic");
        int nr;in>>nr;need(nr>0&&nr<=4,"bounded at most four source rows");Output events(outdir/"all-source-product-back-adds.jsonl",200000);Output results(outdir/"sixteen-evaluated-cell-records.jsonl",100000);
        for(int ri=0;ri<nr;ri++){
            need(seconds()<45,"soft45s complete-cell frontier");int row,nrefs;in>>row>>nrefs;need(row>=0&&row<134,"selected actual top row");vector<P>current(83);vector<Ref>refs(nrefs);
            for(auto&r:refs){in>>r.col>>r.off>>r.bytes>>r.n>>r.sha;need(live+4*r.n+256+META<=LIMIT,"source parser forecast");src.clear();src.seekg(r.off);string payload(size_t(r.bytes),0);src.read(payload.data(),r.bytes);need(src.gcount()==r.bytes&&digest(payload)==r.sha,"thin literal source cell SHA");current[r.col]=binary_decode(payload,r.n);owner();}
            long base=live;int nc;in>>nc;need(nc==4,"only four changed columns");
            for(int ci=0;ci<nc;ci++){
                int col,ns;long cap,nupper;string Ds,Ms;Exp expectedbase;in>>col>>cap>>nupper>>Ds>>Ms;for(auto&v:expectedbase)in>>v;in>>ns;need(col==vector<int>{7,13,27,33}[ci],"fixed free column");P acc;
                for(int step=0;step<ns;step++){
                    int sourcecol,g;long na,nb,k,b,bound;in>>sourcecol>>g>>na>>nb>>k>>b>>bound;need(g>=0&&g<15&&current[sourcecol].n==na&&global[g].n==nb,"actual polynomial source routing");
                    forecast(live-acc.n,k,b,bound);string before=valuehash(acc);P prod=mul(current[sourcecol],global[g]);need(prod.n<=k,"product support");pairs+=na*nb;need(pairs<=PAIRLIMIT,"fixed coefficient pair limit");
                    P next=add(acc,prod);need(next.n<=b,"union support");acc=std::move(next);owner();
                    {P back=add(acc,prod,true);need(valuehash(back)==before,"literal product back-add QQ identity");backadds++;owner();}
                    events.f<<'['<<row<<','<<col<<','<<step<<','<<sourcecol<<','<<g<<','<<acc.n<<",true]\n";
                }
                need(acc.n<=nupper&&acc.shift==expectedbase,"final support/shift");string expected=valuehash(acc);long n=acc.n;string payload;
                if(n)payload=binary_encode(acc);acc=P();need(live==base,"release output before binary decoder");need(live+4*n+256+META<=LIMIT,"binary roundtrip prospective bound");
                string filename="row-"+to_string(row)+"-col-"+to_string(col)+".bin";Output bf(outdir/filename,cap);bf.f<<"Q5CELL1\n";bf.f.put(n?char(1):char(0));if(n)bf.f.write(payload.data(),payload.size());bf.commit();raw();
                P back;if(n)back=binary_decode(payload,n);need(valuehash(back)==expected,"new QQ cell binary roundtrip");owner();
                long count;in>>count;need(count==nupper,"support certificate size");
                fmpz_t D,M,z;fmpz_init(D);fmpz_init(M);fmpz_init(z);need(fmpz_set_str(D,Ds.c_str(),10)==0&&fmpz_set_str(M,Ms.c_str(),10)==0,"integer height bounds");fmpq_t q;fmpq_init(q);
                auto greater=[](const Exp&a,const Exp&b){long da=0,db=0;for(int j=0;j<5;j++){da+=a[j];db+=b[j];}return da!=db?da>db:a>b;};
                long matched=0;Exp previous{};bool first=true;
                for(long k=0;k<count;k++){Exp e;for(auto&v:e)in>>v;need(first||greater(previous,e),"strict DEGLEX support order");previous=e;first=false;
                    if(matched==back.n)continue;auto actual=exponent(back,matched);if(e!=actual){need(greater(e,actual),"actual coefficient belongs to streamed support");continue;}
                    fmpq_mpoly_get_term_coeff_fmpq(q,back.v,matched,ctx);need(fmpz_divisible(D,fmpq_denref(q)),"exact known scalar common denominator");fmpz_divexact(z,D,fmpq_denref(q));fmpz_mul(z,z,fmpq_numref(q));fmpz_abs(z,z);need(fmpz_cmp(z,M)<=0,"actual coefficient integer height");matched++;
                }need(matched==back.n,"all actual terms matched streamed certificate");
                fmpq_clear(q);fmpz_clear(D);fmpz_clear(M);fmpz_clear(z);string h=filehash(outdir/filename);
                results.f<<"{\"row\":"<<row<<",\"column\":"<<col<<",\"file\":\""<<filename<<"\",\"sha256\":\""<<h<<"\",\"terms\":"<<n<<",\"actual_QQ_value_sha256\":\""<<expected<<"\",\"all_source_product_identities\":true,\"binary_roundtrip\":true,\"height_support_verified\":true}\n";
                cellcount++;cout<<"cell "<<cellcount<<" row "<<row<<" col "<<col<<" terms "<<n<<" owned "<<peak+META<<" seconds "<<seconds()<<endl;
            }
            rowsdone++;raw();
        }
        string extra;need(bool(in)&&!(in>>extra),"fixed thin control consumed");events.commit();results.commit();need(filehash(matrixpath)==matrixsha,"source matrix unchanged");owner();
        Output report(outdir/"native-summary.json",20000);report.f<<"{\"all_passed\":true,\"rows\":"<<rowsdone<<",\"cells\":"<<cellcount<<",\"all_source_polynomial_identities\":true,\"all_binary_roundtrips\":true,\"back_adds\":"<<backadds<<",\"coefficient_pairs\":"<<pairs<<",\"owned_coefficients_peak\":"<<peak<<",\"fixed_metadata\":"<<META<<",\"whole_owned_peak\":"<<peak+META<<",\"prospective_whole_peak\":"<<prospect<<",\"seconds\":"<<seconds()<<",\"no_full111_or_inverse_claim\":true}\n";report.commit();raw();
    }
    need(!live,"all coefficients released");fmpq_mpoly_ctx_clear(ctx);return 0;
}
