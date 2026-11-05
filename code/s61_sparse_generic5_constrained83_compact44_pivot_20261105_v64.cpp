// ONE single-file streamed pivot from compact generation39. Never owns a full coefficient matrix.
#include "s61_sparse_generic5_row_stream_kernel_20261104_v63.hpp"
#include <openssl/evp.h>
#include <map>
#include <memory>
#include <stdexcept>
#include <iomanip>
namespace fs=std::filesystem;
long owned_bound=14734,metadata_charge=142,read_bytes=0,hash_bytes=0;
long input_bytes=7432267,index_bytes=0,other_external_bytes=0;
long rows_cap=22473147,norm_cap=45738,other_cap=3500000;
long parse_peak=0,forecast_peak=0,support_peak=0,backchecks=0,product_pairs=0;
// Count every retained integer index/plan/receipt scalar as well as the
// algebraic coefficient owners. Hash characters are byte/memory metadata.
long index_scalar_slots=0;
fs::path outdir,workdir;string sourcepath,indexpath,gatepath;
void need(bool b,const string&s){if(!b)throw runtime_error(s);}
#include "s61_sparse_generic5_primitive_binary_native_20261105_v64.hpp"
long verified_binary_cells=0;
string digest(const string&s){
    unsigned char md[EVP_MAX_MD_SIZE];unsigned n=0;
    need(EVP_Digest(s.data(),s.size(),md,&n,EVP_sha256(),nullptr)==1,"SHA payload");
    ostringstream o;for(unsigned j=0;j<n;j++)o<<hex<<setw(2)<<setfill('0')<<int(md[j]);return o.str();
}
string filehash(const fs::path&p){
    ifstream f(p,ios::binary);need(bool(f),"SHA open");
    unique_ptr<EVP_MD_CTX,decltype(&EVP_MD_CTX_free)> c(EVP_MD_CTX_new(),EVP_MD_CTX_free);
    need(bool(c)&&EVP_DigestInit_ex(c.get(),EVP_sha256(),nullptr)==1,"SHA init");
    char b[65536];while(f){f.read(b,sizeof(b));auto n=f.gcount();if(n){hash_bytes+=n;need(EVP_DigestUpdate(c.get(),b,n)==1,"SHA update");}}
    unsigned char md[EVP_MAX_MD_SIZE];unsigned n=0;need(EVP_DigestFinal_ex(c.get(),md,&n)==1,"SHA final");
    ostringstream o;for(unsigned j=0;j<n;j++)o<<hex<<setw(2)<<setfill('0')<<int(md[j]);return o.str();
}
struct Ref{int row,col;long offset,bytes,terms;string sha;};
struct Update{int col;long a,k,b,resident,peak,mbits,dbits;};
struct Plan{int row,factor;long terms,newterms,newbytes;vector<Update>updates;};
vector<Ref> refs;map<int,vector<int>> rowrefs;vector<int> allrows,allcols;vector<Plan> plans;
const vector<int>& rowindices(int r){static const vector<int> empty;auto it=rowrefs.find(r);return it==rowrefs.end()?empty:it->second;}
// A bounded stream does not construct an intermediate row string.
struct BoundBuf:streambuf{
    ofstream f;long written=0,cap;
    BoundBuf(const fs::path&p,long c):f(p,ios::binary),cap(c){need(bool(f),"output open");}
    streamsize xsputn(const char*s,streamsize n)override{need(written+n<=cap,"per-file RAW guard");f.write(s,n);need(bool(f),"output write");written+=n;return n;}
    int overflow(int c)override{if(c==traits_type::eof())return traits_type::not_eof(c);char x=char(c);xsputn(&x,1);return c;}
    int sync()override{f.flush();return f?0:-1;}
};
struct Output{
    fs::path path,tmp;BoundBuf b;ostream o;
    Output(const fs::path&p,long cap):path(p),tmp(p.string()+".partial"),b(tmp,cap),o(&b){o.exceptions(ios::badbit|ios::failbit);need(!fs::exists(path),"immutable output exists");}
    void commit(){o.flush();b.f.close();need(bool(b.f),"close output");fs::rename(tmp,path);}
};
long treebytes(const fs::path&p){long z=0;for(auto&e:fs::recursive_directory_iterator(p))if(e.is_regular_file())z+=e.file_size();return z;}
void ledger(){
    long rb=fs::exists(outdir/"rows")?treebytes(outdir/"rows"):0;
    long nb=0;for(auto&e:fs::directory_iterator(outdir))if(e.is_regular_file()&&e.path().filename().string().rfind("normalized-pivot-row",0)==0)nb+=e.file_size();
    long ob=treebytes(outdir)+treebytes(workdir)+other_external_bytes-rb-nb;
    need(rb<=rows_cap&&nb<=norm_cap&&ob<=other_cap,"whole RAW category guard");
    need(input_bytes+index_bytes+rb+nb+ob<=33451152&&input_bytes+index_bytes+rb+nb+ob<=BYTES,"whole RAW total guard");
}
ifstream source;
P load(const Ref&r){
    long parsed=live+4*r.terms+256;parse_peak=max(parse_peak,parsed);need(parsed<=owned_bound&&parsed<=CAP,"binary parser/value-roundtrip owned guard");
    source.clear();source.seekg(r.offset);string bytes(size_t(r.bytes),'\0');source.read(bytes.data(),r.bytes);need(source.gcount()==r.bytes,"binary cell read");read_bytes+=r.bytes;
    need(digest(bytes)==r.sha,"binary source cell SHA");P p=binary_decode(bytes,r.terms);verified_binary_cells++;return p;
}
using Row=vector<P>;
Row loadrow(int r){Row a(83);for(int k:rowindices(r)){auto&t=refs[k];a[t.col]=load(t);}return a;}
long terms(const Row&a){long n=0;for(auto&p:a)n+=p.n;return n;}
void rowfile(const fs::path&p,const Row&a,int row,long cap,bool normalized){
    Output f(p,cap);f.o<<"{\"original_row\":"<<row<<",\"pivot_original_column\":68,\"normalized\":"<<(normalized?"true":"false")<<",\"entries\":[";
    bool first=true;for(int j:allcols)if(nz(a[j])){if(!first)f.o<<',';first=false;f.o<<'['<<j<<',';put(f.o,a[j]);f.o<<']';}f.o<<"],\"total_terms\":"<<terms(a)<<"}\n";f.commit();
}
struct Saved{string path,sha;long terms,bytes;};map<int,Saved> saved;
void manifest(const string&name,bool complete,const string&sourcehash,const string&indexhash,const string&normhash){
    Output f(outdir/name,180000);auto&o=f.o;
    o<<"{\"complete_pivot\":"<<(complete?"true":"false")<<",\"shape\":["<<(complete?191:192)<<','<<(complete?38:39)<<"],\"source_file\":\""<<sourcepath<<"\",\"source_sha256\":\""<<sourcehash<<"\",\"index_sha256\":\""<<indexhash<<"\",\"pivot\":[32,68],\"normalized_pivot_sha256\":\""<<normhash<<"\",\"rows\":[";
    bool first=true;long logical=424697;for(int r:allrows){if(complete&&r==32)continue;if(!first)o<<',';first=false;o<<"{\"original_row\":"<<r;
        if(saved.count(r)){auto&t=saved.at(r);o<<",\"literal_row_file\":\""<<t.path<<"\",\"sha256\":\""<<t.sha<<"\",\"terms\":"<<t.terms;long old=0;for(int k:rowindices(r))old+=refs[k].terms;logical+=t.terms-old;}
        else if(!normhash.empty()&&r==32){o<<",\"normalized_pivot_row\":\"normalized-pivot-row.json\",\"sha256\":\""<<normhash<<"\",\"terms\":1315";}
        else{o<<",\"indexed_source_cell_indices\":[";bool z=true;for(int k:rowindices(r)){if(!z)o<<',';z=false;o<<k;}o<<']';}o<<'}';}
    if(complete)logical-=1315;o<<"],\"column_original_ids\":[";first=true;for(int j:allcols)if(!complete||j!=68){if(!first)o<<',';first=false;o<<j;}o<<"],\"logical_total_terms\":"<<logical<<",\"all_coefficients_preserved\":true}\n";f.commit();
}
int main(int argc,char**argv){
    need(argc==4,"arguments");outdir=argv[2];workdir=argv[3];fmpq_mpoly_ctx_init(ctx,5,ORD_DEGLEX);
    string sh,ih,gh;int nr,nc,nrefs,nplans;ifstream in(argv[1]);
    in>>sourcepath>>sh>>indexpath>>ih>>gatepath>>gh>>other_external_bytes;
    in>>nr;allrows.resize(nr);for(auto&v:allrows)in>>v;in>>nc;allcols.resize(nc);for(auto&v:allcols)in>>v;
    in>>nrefs;refs.resize(nrefs);for(int k=0;k<nrefs;k++){auto&r=refs[k];in>>r.row>>r.col>>r.offset>>r.bytes>>r.terms>>r.sha;rowrefs[r.row].push_back(k);}
    in>>nplans;plans.resize(nplans);for(auto&p:plans){int n;in>>p.row>>p.terms>>p.newterms>>p.newbytes>>p.factor>>n;p.updates.resize(n);for(auto&u:p.updates)in>>u.col>>u.a>>u.k>>u.b>>u.resident>>u.peak>>u.mbits>>u.dbits;}
    need(bool(in)&&nr==192&&nc==39&&nrefs==4457&&nplans==93,"control input");
    long nu=0;for(auto&p:plans)nu+=p.updates.size();
    index_scalar_slots=5*nrefs+nrefs+long(rowrefs.size())+nr+nc+5*nplans+8*nu+3*nplans+512;
    need(index_scalar_slots+owned_bound<=CAP,"index/plan + algebraic scratch owned guard");
    need(filehash(sourcepath)==sh&&filehash(indexpath)==ih&&filehash(gatepath)==gh,"whole source/index/gate SHA");
    fs::create_directory(outdir/"rows");fs::create_directory(outdir/"receipts");manifest("generation-0.json",false,sh,ih,"");ledger();
    source.open(sourcepath,ios::binary);need(bool(source),"open source");{string magic(8,0);source.read(magic.data(),8);need(magic=="Q5PRIM1\n","binary file magic");}live=metadata_charge;guard();
    string normhash;bool complete=false;long completed=0,logical=424697;
    {
        Row pivot=loadrow(32);need(terms(pivot)==1315,"pivot row369");auto&p=pivot[68];need(p.n==1&&exponent(p,0)==Exp{1,0,0,0,0},"selected monomial");
        fmpq_t q;fmpq_init(q);fmpq_set_si(q,2048,1215);fmpq_t qc;fmpq_init(qc);fmpq_mpoly_get_term_coeff_fmpq(qc,p.v,0,ctx);need(fmpq_equal(q,qc),"unit QQ2048/1215");Exp pe=exponent(p,0);P unit(p);
        for(int j:allcols)if(nz(pivot[j])){P old(pivot[j]);for(int v=0;v<5;v++)pivot[j].shift[v]-=pe[v];fmpq_mpoly_scalar_div_fmpq(pivot[j].v,pivot[j].v,q,ctx);pivot[j].sync();P back=mul(pivot[j],unit);need(equal(back,old),"normalization back-add");backchecks++;}
        fmpq_clear(q);fmpq_clear(qc);unit=P();need(pivot[68].n==1&&exponent(pivot[68],0)==Exp{},"normalized1");
        rowfile(outdir/"normalized-pivot-row.json",pivot,32,norm_cap,true);normhash=filehash(outdir/"normalized-pivot-row.json");ledger();
        {Output ec(outdir/"event-context.json",4096);ec.o<<"{\"schema\":\"[row,column,new_terms,back_add_true]\",\"pivot\":[32,68],\"normalized_pivot_file\":\"normalized-pivot-row.json\",\"normalized_pivot_sha256\":\""<<normhash<<"\",\"source_gate_sha256\":\""<<gh<<"\",\"source_index_sha256\":\""<<ih<<"\",\"reconstruction\":\"source(row,col)-source(row,68)*normalized_pivot(col)\",\"factor_source_cell_indices\":[";bool first=true;for(auto&p:plans){if(!first)ec.o<<',';first=false;ec.o<<'['<<p.row<<','<<p.factor<<']';}ec.o<<"]}\n";ec.commit();ledger();}
        ofstream events(outdir/"arithmetic-events.jsonl",ios::binary);need(bool(events),"events open");
        for(const auto&plan:plans){
            if(seconds()>42)break;need(live==metadata_charge+1315,"released previous row");
            Row current=loadrow(plan.row);need(terms(current)==plan.terms,"source row total");need(refs[plan.factor].row==plan.row&&refs[plan.factor].col==68,"factor ref");P&factor=current[68];need(factor.n==refs[plan.factor].terms,"factor nonzero");
            for(const auto&u:plan.updates){
                need(current[u.col].n==u.a&&nz(pivot[u.col]),"update source counts");
                long prospective=live+2*u.k+3*u.b+96;forecast_peak=max(forecast_peak,prospective);need(prospective<=u.peak&&prospective<=owned_bound&&prospective<=CAP,"whole update forecast");
                {
                    auto fsu=support(factor),psu=support(pivot[u.col]),oldsu=support(current[u.col]);set<Exp>prod,sun;
                    for(auto&a:fsu)for(auto&b:psu){Exp z;for(int v=0;v<5;v++)z[v]=a[v]+b[v];prod.insert(z);}
                    sun.insert(oldsu.begin(),oldsu.end());sun.insert(prod.begin(),prod.end());
                    long sp=live+fsu.size()+psu.size()+oldsu.size()+prod.size()+sun.size()+96;support_peak=max(support_peak,sp);
                    need(prod.size()==size_t(u.k)&&sun.size()==size_t(u.b)&&sp<=prospective,"exact support union forecast");product_pairs+=factor.n*pivot[u.col].n;
                }
                P prod=mul(factor,pivot[u.col]);P next=add(current[u.col],prod,true);
                {P back=add(next,prod);need(equal(back,current[u.col]),"literal back-add");backchecks++;}
                need(next.n<=u.b&&peak+96<=owned_bound,"actual owned scratch");current[u.col]=std::move(next);
                events<<'['<<plan.row<<','<<u.col<<','<<current[u.col].n<<",true]\n";need(events.tellp()<=1048576,"compact events byte cap");events.flush();
            }
            long nf=factor.n;factor=P();need(terms(current)<=plan.newterms,"whole new row support bound");
            string filename="rows/row-"+to_string(plan.row)+".json";rowfile(outdir/filename,current,plan.row,plan.newbytes,false);
            Saved s{filename,filehash(outdir/filename),terms(current),long(fs::file_size(outdir/filename))};saved[plan.row]=s;completed++;logical+=s.terms-plan.terms;
            Output receipt(outdir/("receipts/row-"+to_string(plan.row)+".json"),1024);receipt.o<<"{\"original_row\":"<<plan.row<<",\"factor_source_cell_index\":"<<plan.factor<<",\"all_cell_back_add_identities_checked\":true,\"pivot_column_cleared\":true,\"literal_row_sha256\":\""<<s.sha<<"\",\"terms\":"<<s.terms<<",\"RAW_bytes\":"<<s.bytes<<"}\n";receipt.commit();ledger();
            cout<<"row "<<completed<<"/93 original "<<plan.row<<" terms "<<s.terms<<" owned "<<peak<<" seconds "<<seconds()<<endl;
        }
        events.close();complete=completed==93;manifest(complete?"generation-1.json":"exact-prefix-generation.json",complete,sh,ih,normhash);ledger();
    }
    need(live==metadata_charge,"row ownership fully released");long owned=max({peak+96,parse_peak,forecast_peak,support_peak});need(owned<=owned_bound&&owned<=CAP,"final owned bound");
    if(complete)logical-=1315;
    need(filehash(sourcepath)==sh&&filehash(indexpath)==ih&&filehash(gatepath)==gh,"source/index/gate unchanged");
    Output summary(outdir/"native-summary.json",12000);summary.o<<"{\"all_passed\":true,\"completed_pivot\":"<<(complete?"true":"false")<<",\"updated_rows\":"<<completed<<",\"pivot\":[32,68],\"shape\":["<<(complete?191:192)<<','<<(complete?38:39)<<"],\"logical_input_terms\":424697,\"logical_after_terms\":"<<logical<<",\"logical_forecast_upper\":467576,\"owned_coefficients_peak\":"<<peak<<",\"owned_including_scratch_peak\":"<<owned<<",\"index_plan_receipt_integer_slots_upper\":"<<index_scalar_slots<<",\"whole_owned_scalar_and_scratch_upper\":"<<owned+index_scalar_slots<<",\"parse_phase_peak\":"<<parse_peak<<",\"support_phase_peak\":"<<support_peak<<",\"update_forecast_peak\":"<<forecast_peak<<",\"metadata_scalar_charge\":142,\"back_add_identities_checked\":"<<backchecks<<",\"coefficient_pair_operations\":"<<product_pairs<<",\"compact_event_schema\":true,\"native_binary_value_roundtrips\":"<<verified_binary_cells<<",\"source_payload_read_bytes\":"<<read_bytes<<",\"hash_read_bytes\":"<<hash_bytes<<",\"seconds\":"<<setprecision(12)<<seconds()<<",\"no_full_coefficient_matrix\":true,\"no_second_pivot\":true,\"old44_steps_not_replayed\":true,\"stop_reason\":\""<<(complete?"one_authorized_pivot_complete":"soft42s_exact_prefix_saved")<<"\"}\n";summary.commit();ledger();
    live-=metadata_charge;fmpq_mpoly_ctx_clear(ctx);return 0;
}
