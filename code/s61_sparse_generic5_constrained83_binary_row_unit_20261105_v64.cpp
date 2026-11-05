// ONE generic source-gated binary-row unit pivot. No full coefficient matrix is owned.
#include "s61_sparse_generic5_row_stream_kernel_20261104_v63.hpp"
#include <openssl/evp.h>
#include <map>
#include <memory>
#include <stdexcept>
#include <iomanip>
namespace fs=std::filesystem;
long owned_bound=0,metadata_charge=0,read_bytes=0,hash_bytes=0;
int PR,PC,PRIOR;long logical_input,pivot_terms,logical_forecast,raw_total;string unit_q;Exp unit_exp;
long input_bytes=0,index_bytes=0,other_external_bytes=0;
long rows_cap=0,norm_cap=0,other_cap=0;
long parse_peak=0,forecast_peak=0,support_peak=0,backchecks=0,product_pairs=0,writer_peak=0,writer_roundtrips=0;
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
struct BinaryBound{int col;long bytes,terms,indexbytes;};
struct Plan{int row,factor;long terms,newterms,newbytes,binbytes,indexbytes,writerpeak;vector<Update>updates;vector<BinaryBound>binarycells;};
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
    need(input_bytes+index_bytes+rb+nb+ob<=raw_total&&input_bytes+index_bytes+rb+nb+ob<=BYTES,"whole RAW total guard");
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
    Output f(p,cap);f.o<<"{\"original_row\":"<<row<<",\"pivot_original_column\":"<<PC<<",\"normalized\":"<<(normalized?"true":"false")<<",\"entries\":[";
    bool first=true;for(int j:allcols)if(nz(a[j])){if(!first)f.o<<',';first=false;f.o<<'['<<j<<',';put(f.o,a[j]);f.o<<']';}f.o<<"],\"total_terms\":"<<terms(a)<<"}\n";f.commit();
}
struct Saved{string path,sha,indexpath,indexsha;long terms,bytes,indexbytes;};map<int,Saved> saved;

Saved binaryrow(const Row&a,const Plan&p){
    string rel="rows/row-"+to_string(p.row)+".bin",idxrel="rows/row-"+to_string(p.row)+".index.json";
    Output blob(outdir/rel,p.binbytes),index(outdir/idxrel,p.indexbytes);blob.o<<"Q5PRIM1\n";
    index.o<<"{\"encoding\":\"primitive-int-exponent-v1\",\"original_row\":"<<p.row<<",\"pivot_original_column\":"<<PC<<",\"source_binary_file\":\"row-"<<p.row<<".bin\",\"entries\":[";
    long covered=0;set<int> coveredcols;for(const auto&b:p.binarycells){need(coveredcols.insert(b.col).second,"unique binary output column");covered+=a[b.col].n;}need(covered==terms(a),"all computed nonzero QQ cells are indexed");
    bool first=true;
    for(const auto&b:p.binarycells){
        const P&cell=a[b.col];need(cell.n<=b.terms,"binary cell support bound");if(!cell.n)continue;
        long prospective=live+4*cell.n+256;writer_peak=max(writer_peak,prospective);need(prospective<=p.writerpeak&&prospective<=owned_bound,"writer full owned scratch");
        string payload=binary_encode(cell);need(long(payload.size())<=b.bytes,"binary per-cell bytes");long offset=blob.b.written;
        blob.o.write(payload.data(),payload.size());
        {P back=binary_decode(payload,cell.n);need(equal(back,cell),"stored binary cell equals computed QQ cell");writer_roundtrips++;}
        if(!first)index.o<<',';first=false;long before=index.b.written;
        index.o<<'['<<b.col<<','<<offset<<','<<payload.size()<<",\""<<digest(payload)<<"\","<<cell.n<<',';
        if(cell.n==1){fmpq_t q;fmpq_init(q);fmpq_mpoly_get_term_coeff_fmpq(q,cell.v,0,ctx);index.o<<'"'<<qs(q)<<"\",[";auto e=exponent(cell,0);for(int v=0;v<5;v++){if(v)index.o<<',';index.o<<e[v];}index.o<<']';fmpq_clear(q);}
        else index.o<<"null,null";
        index.o<<']';need(index.b.written-before<=b.indexbytes,"paired-index entry bound");
    }
    blob.commit();string h=filehash(outdir/rel);index.o<<"],\"total_terms\":"<<terms(a)<<",\"source_binary_sha256\":\""<<h<<"\"}\n";index.commit();
    Saved s{rel,h,idxrel,filehash(outdir/idxrel),terms(a),long(fs::file_size(outdir/rel)),long(fs::file_size(outdir/idxrel))};
    need(s.bytes<=p.binbytes&&s.indexbytes<=p.indexbytes&&s.bytes+s.indexbytes<=p.newbytes,"whole paired binary-row bound");return s;
}
void manifest(const string&name,bool complete,const string&sourcehash,const string&indexhash,const string&normhash){
    Output f(outdir/name,180000);auto&o=f.o;
    o<<"{\"complete_pivot\":"<<(complete?"true":"false")<<",\"shape\":["<<(allrows.size()-(complete?1:0))<<','<<(allcols.size()-(complete?1:0))<<"],\"source_file\":\""<<sourcepath<<"\",\"source_sha256\":\""<<sourcehash<<"\",\"index_sha256\":\""<<indexhash<<"\",\"pivot\":["<<PR<<','<<PC<<"],\"normalized_pivot_sha256\":\""<<normhash<<"\",\"rows\":[";
    bool first=true;long logical=logical_input;for(int r:allrows){if(complete&&r==PR)continue;if(!first)o<<',';first=false;o<<"{\"original_row\":"<<r;
        if(saved.count(r)){auto&t=saved.at(r);o<<",\"literal_binary_row_file\":\""<<t.path<<"\",\"sha256\":\""<<t.sha<<"\",\"paired_index_file\":\""<<t.indexpath<<"\",\"paired_index_sha256\":\""<<t.indexsha<<"\",\"terms\":"<<t.terms;long old=0;for(int k:rowindices(r))old+=refs[k].terms;logical+=t.terms-old;}
        else if(!normhash.empty()&&r==PR){o<<",\"normalized_pivot_row\":\"normalized-pivot-row.json\",\"sha256\":\""<<normhash<<"\",\"terms\":"<<pivot_terms<<"";}
        else{o<<",\"indexed_source_cell_indices\":[";bool z=true;for(int k:rowindices(r)){if(!z)o<<',';z=false;o<<k;}o<<']';}o<<'}';}
    if(complete)logical-=pivot_terms;o<<"],\"column_original_ids\":[";first=true;for(int j:allcols)if(!complete||j!=PC){if(!first)o<<',';first=false;o<<j;}o<<"],\"logical_total_terms\":"<<logical<<",\"all_coefficients_preserved\":true}\n";f.commit();
}
int main(int argc,char**argv){
    need(argc==4,"arguments");outdir=argv[2];workdir=argv[3];fmpq_mpoly_ctx_init(ctx,5,ORD_DEGLEX);
    string sh,ih,gh,proofpath,proofsha,boundspath,boundssha;int nr,nc,nrefs,nplans;ifstream in(argv[1]);
    in>>sourcepath>>sh>>indexpath>>ih>>gatepath>>gh>>proofpath>>proofsha>>boundspath>>boundssha>>other_external_bytes;
    in>>PR>>PC>>PRIOR>>logical_input>>pivot_terms>>logical_forecast>>owned_bound>>metadata_charge>>input_bytes>>rows_cap>>norm_cap>>other_cap>>raw_total>>unit_q;for(auto&v:unit_exp)in>>v;
    in>>nr;allrows.resize(nr);for(auto&v:allrows)in>>v;in>>nc;allcols.resize(nc);for(auto&v:allcols)in>>v;
    in>>nrefs;refs.resize(nrefs);for(int k=0;k<nrefs;k++){auto&r=refs[k];in>>r.row>>r.col>>r.offset>>r.bytes>>r.terms>>r.sha;rowrefs[r.row].push_back(k);}
    in>>nplans;plans.resize(nplans);for(auto&p:plans){int n;in>>p.row>>p.terms>>p.newterms>>p.newbytes>>p.factor>>n;p.updates.resize(n);for(auto&u:p.updates)in>>u.col>>u.a>>u.k>>u.b>>u.resident>>u.peak>>u.mbits>>u.dbits;int bn;in>>p.binbytes>>p.indexbytes>>p.writerpeak>>bn;p.binarycells.resize(bn);for(auto&b:p.binarycells)in>>b.col>>b.bytes>>b.terms>>b.indexbytes;}
    need(bool(in)&&nr==236-PRIOR&&nc==83-PRIOR&&nrefs>0&&nplans>0,"control input");
    long nu=0,nbc=0;for(auto&p:plans){nu+=p.updates.size();nbc+=p.binarycells.size();}
    index_scalar_slots=5*nrefs+nrefs+long(rowrefs.size())+nr+nc+5*nplans+8*nu+3*nplans+512+4*nbc+4*nplans+512;need(index_scalar_slots>0,"all binary-writer metadata slots");
    need(index_scalar_slots+owned_bound<=CAP,"index/plan + algebraic scratch owned guard");
    need(filehash(sourcepath)==sh&&filehash(indexpath)==ih&&filehash(gatepath)==gh&&filehash(proofpath)==proofsha&&filehash(boundspath)==boundssha,"whole source/index/gate/binary-proof SHA");
    fs::create_directory(outdir/"rows");fs::create_directory(outdir/"receipts");manifest("generation-0.json",false,sh,ih,"");ledger();
    source.open(sourcepath,ios::binary);need(bool(source),"open source");{string magic(8,0);source.read(magic.data(),8);need(magic=="Q5PRIM1\n","binary file magic");}live=metadata_charge;guard();
    string normhash;bool complete=false;long completed=0,logical=logical_input;
    {
        Row pivot=loadrow(PR);need(terms(pivot)==pivot_terms,"pivot row369");auto&p=pivot[PC];need(p.n==1&&exponent(p,0)==unit_exp,"selected monomial");
        fmpq_t q;fmpq_init(q);need(fmpq_set_str(q,unit_q.c_str(),10)==0,"parse unit QQ");fmpq_t qc;fmpq_init(qc);fmpq_mpoly_get_term_coeff_fmpq(qc,p.v,0,ctx);need(fmpq_equal(q,qc),"source parameter-monomial unit QQ");Exp pe=exponent(p,0);P unit(p);
        for(int j:allcols)if(nz(pivot[j])){P old(pivot[j]);for(int v=0;v<5;v++)pivot[j].shift[v]-=pe[v];fmpq_mpoly_scalar_div_fmpq(pivot[j].v,pivot[j].v,q,ctx);pivot[j].sync();P back=mul(pivot[j],unit);need(equal(back,old),"normalization back-add");backchecks++;}
        fmpq_clear(q);fmpq_clear(qc);unit=P();need(pivot[PC].n==1&&exponent(pivot[PC],0)==Exp{},"normalized1");
        rowfile(outdir/"normalized-pivot-row.json",pivot,PR,norm_cap,true);normhash=filehash(outdir/"normalized-pivot-row.json");ledger();
        {Output ec(outdir/"event-context.json",4096);ec.o<<"{\"schema\":\"[row,column,new_terms,back_add_true]\",\"pivot\":["<<PR<<','<<PC<<"],\"normalized_pivot_file\":\"normalized-pivot-row.json\",\"normalized_pivot_sha256\":\""<<normhash<<"\",\"source_gate_sha256\":\""<<gh<<"\",\"source_index_sha256\":\""<<ih<<"\",\"reconstruction\":\"source(row,col)-source(row,pivot_column)*normalized_pivot(col)\",\"factor_source_cell_indices\":[";bool first=true;for(auto&p:plans){if(!first)ec.o<<',';first=false;ec.o<<'['<<p.row<<','<<p.factor<<']';}ec.o<<"]}\n";ec.commit();ledger();}
        ofstream events(outdir/"arithmetic-events.jsonl",ios::binary);need(bool(events),"events open");
        for(const auto&plan:plans){
            if(seconds()>42)break;need(live==metadata_charge+pivot_terms,"released previous row");
            Row current=loadrow(plan.row);need(terms(current)==plan.terms,"source row total");need(refs[plan.factor].row==plan.row&&refs[plan.factor].col==PC,"factor ref");P&factor=current[PC];need(factor.n==refs[plan.factor].terms,"factor nonzero");
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
            Saved s=binaryrow(current,plan);saved[plan.row]=s;completed++;logical+=s.terms-plan.terms;
            Output receipt(outdir/("receipts/row-"+to_string(plan.row)+".json"),2048);receipt.o<<"{\"original_row\":"<<plan.row<<",\"factor_source_cell_index\":"<<plan.factor<<",\"all_cell_back_add_identities_checked\":true,\"pivot_column_cleared\":true,\"literal_row_sha256\":\""<<s.sha<<"\",\"paired_index_sha256\":\""<<s.indexsha<<"\",\"binary_cell_value_roundtrips\":true,\"terms\":"<<s.terms<<",\"RAW_bytes\":"<<s.bytes+s.indexbytes<<"}\n";receipt.commit();ledger();
            cout<<"row "<<completed<<'/'<<plans.size()<<" original "<<plan.row<<" terms "<<s.terms<<" owned "<<peak<<" seconds "<<seconds()<<endl;
        }
        events.close();complete=completed==long(plans.size());manifest(complete?"generation-1.json":"exact-prefix-generation.json",complete,sh,ih,normhash);ledger();
    }
    need(live==metadata_charge,"row ownership fully released");long owned=max({peak+96,parse_peak,forecast_peak,support_peak,writer_peak});need(owned<=owned_bound&&owned<=CAP,"final owned bound");
    if(complete)logical-=pivot_terms;
    need(filehash(sourcepath)==sh&&filehash(indexpath)==ih&&filehash(gatepath)==gh&&filehash(proofpath)==proofsha&&filehash(boundspath)==boundssha,"source and all gate files unchanged");
    Output summary(outdir/"native-summary.json",12000);summary.o<<"{\"all_passed\":true,\"completed_pivot\":"<<(complete?"true":"false")<<",\"updated_rows\":"<<completed<<",\"pivot\":["<<PR<<','<<PC<<"],\"shape\":["<<(allrows.size()-(complete?1:0))<<','<<(allcols.size()-(complete?1:0))<<"],\"logical_input_terms\":"<<logical_input<<",\"logical_after_terms\":"<<logical<<",\"logical_forecast_upper\":"<<logical_forecast<<",\"owned_coefficients_peak\":"<<peak<<",\"owned_including_scratch_peak\":"<<owned<<",\"index_plan_receipt_integer_slots_upper\":"<<index_scalar_slots<<",\"whole_owned_scalar_and_scratch_upper\":"<<owned+index_scalar_slots<<",\"parse_phase_peak\":"<<parse_peak<<",\"support_phase_peak\":"<<support_peak<<",\"update_forecast_peak\":"<<forecast_peak<<",\"metadata_scalar_charge\":"<<metadata_charge<<",\"back_add_identities_checked\":"<<backchecks<<",\"coefficient_pair_operations\":"<<product_pairs<<",\"compact_event_schema\":true,\"new_binary_row_output\":true,\"writer_value_roundtrips\":"<<writer_roundtrips<<",\"writer_owned_peak\":"<<writer_peak<<",\"native_binary_value_roundtrips\":"<<verified_binary_cells<<",\"source_payload_read_bytes\":"<<read_bytes<<",\"hash_read_bytes\":"<<hash_bytes<<",\"seconds\":"<<setprecision(12)<<seconds()<<",\"no_full_coefficient_matrix\":true,\"no_second_pivot\":true,\"old"<<PRIOR<<"_steps_not_replayed\":true,\"stop_reason\":\""<<(complete?"one_authorized_pivot_complete":"soft42s_exact_prefix_saved")<<"\"}\n";summary.commit();ledger();
    live-=metadata_charge;fmpq_mpoly_ctx_clear(ctx);return 0;
}
