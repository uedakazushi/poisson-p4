// READ ONLY: exact source/support/height forecasts, no row operation calls.
#define main dormant_shortlist_arithmetic_main
#include "s61_sparse_generic5_constrained83_shortlist_continuation_20261104_v63.cpp"
#undef main
int main(int argc,char**argv){
    assert(argc==3);epoch=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,5,ORD_DEGLEX);int status=0;
    {
        try{
            string dir=argv[2];raw_input_bytes=filesystem::file_size(argv[1]);ifstream in(argv[1]);int nv,nr,nc,nar,nac,nn;
            in>>nv>>nr>>nc>>nar>>nac>>nn;assert(nv==5&&nr==236&&nc==83&&nar==200&&nac==47&&nn==5483);
            vector<bool>ar(236,false),ac(83,false);for(int k=0;k<nar;k++){int i;in>>i;ar[i]=true;}for(int k=0;k<nac;k++){int j;in>>j;ac[j]=true;}
            Matrix A(236,vector<P>(83));for(int k=0;k<nn;k++){int i,j;long n;Exp sh;in>>i>>j>>n;for(auto&e:sh)in>>e;A[i][j]=read(in,n);A[i][j].shift=sh;}
            int ne;in>>ne;assert(ne==16);set<pair<int,int>>excluded;for(int k=0;k<ne;k++){int i,j;in>>i>>j;excluded.emplace(i,j);}string tail;assert(in&&!(in>>tail));assert(live==143477);
            auto options=candidates(A,ar,ac);assert(options.size()==328);set<pair<int,int>>first16;for(int k=0;k<16;k++)first16.emplace(options[k].r,options[k].c);assert(first16==excluded);
            vector<Pick>remaining;for(const auto&p:options)if(!excluded.count({p.r,p.c}))remaining.push_back(p);assert(remaining.size()==312);
            StatMatrix stats=statistics(A);long completed=0,unsafe=0;bool found=false;string reason;int next=0;
            for(int k=0;k<int(remaining.size());k++){
                next=k;if(seconds()>32){reason="source_soft_time_guard_before_candidate";break;}
                if(raw_input_bytes+bytes_on_disk(dir)+2*1024*1024+8192>BYTES){reason="source_RAW_guard_before_candidate";break;}
                const Pick&p=remaining[k];Forecast f=forecast(A,stats,ar,ac,p,dir,CAP+1,999,999,4096);
                if(f.reason=="time_forecast_before_pivot"){reason="source_time_guard_in_incomplete_candidate";break;}
                string number=to_string(k);while(number.size()<3)number="0"+number;
                ofstream out(dir+"/forecast-"+number+".json");
                out<<"{\"candidate_cursor\":"<<k<<",\"original_row\":"<<p.r<<",\"original_column\":"<<p.c<<",\"completed_source_forecast\":true,\"FULL_safe\":"<<(f.safe?"true":"false")<<",\"reason\":\""<<(f.safe?"all_full_row_updates_fit":f.reason)<<"\",\"peak_forecast\":"<<f.peakbound<<",\"complete_RAW_forecast\":"<<f.outbound<<",\"final_terms_upper\":"<<f.live_end<<",\"all_completed_update_forecasts\":[";
                bool co=false;for(int i=0;i<236;i++)if(!f.byrow[i].empty()){
                    if(co)out<<',';co=true;out<<'['<<i<<",[";bool cj=false;for(const auto&u:f.byrow[i]){if(cj)out<<',';cj=true;out<<'['<<u.col<<','<<u.oldn<<','<<u.k<<','<<u.b<<','<<u.liveupper<<','<<u.peakbound<<','<<u.h.N<<','<<u.h.D<<']';}out<<"]]";
                }
                out<<"],\"no_coefficient_update_executed\":true}\n";out.close();assert(out.good());
                completed++;if(!f.safe)unsafe++;next=k+1;
                ofstream cursor(dir+"/cursor-"+number+".json");cursor<<"{\"completed_remaining_candidates\":"<<completed<<",\"next_cursor\":"<<next<<",\"remaining_total\":312,\"matrix_unchanged\":true}\n";cursor.close();
                assert(raw_input_bytes+bytes_on_disk(dir)<BYTES);assert(live==143477);
                if(f.safe){
                    fmpq_t q;fmpq_init(q);fmpq_mpoly_get_term_coeff_fmpq(q,A[p.r][p.c].v,0,ctx);Exp e=exponent(A[p.r][p.c],0);
                    ofstream chosen(dir+"/chosen-candidate.txt");chosen<<p.r<<' '<<p.c<<' '<<qs(q);for(long x:e)chosen<<' '<<x;chosen<<'\n';chosen.close();
                    fmpq_clear(q);found=true;reason="first_FULL_safe_remaining_candidate";break;
                }
                if(k+1==int(remaining.size()))reason="all312_remaining_rigorously_unsafe_under_frozen_forecast";
            }
            ofstream final(dir+"/source-frontier-summary.json");
            final<<"{\"all_passed\":true,\"source_only\":true,\"source_owned_terms_unchanged\":143477,\"peak_owned_polynomial_terms\":"<<peak<<",\"excluded_previous_exact16\":true,\"eligible_total\":328,\"remaining_total\":312,\"completed_remaining\":"<<completed<<",\"rigorous_unsafe_remaining\":"<<unsafe<<",\"next_cursor\":"<<next<<",\"FULL_safe_found\":"<<(found?"true":"false")<<",\"all328_unsafe_under_frozen_forecast\":"<<(!found&&completed==312?"true":"false")<<",\"stop_reason\":\""<<reason<<"\",\"native_seconds\":"<<seconds()<<"}\n";final.close();assert(final.good());assert(raw_input_bytes+bytes_on_disk(dir)<BYTES);
        }catch(const exception&e){cerr<<"SOURCE_STOP "<<e.what()<<'\n';status=2;}
    }
    assert(live==0);fmpq_mpoly_ctx_clear(ctx);flint_cleanup();return status;
}
