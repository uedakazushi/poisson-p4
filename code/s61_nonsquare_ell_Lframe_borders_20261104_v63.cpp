// Only the new 317-term L-frame input; never reconstruct the old dense19 stage.
#define main frozen_Lframe_preflight_main
#include "s61_nonsquare_ell_Lframe_preflight_20261104_v63.cpp"
#undef main
#include <cstdio>
#include <climits>
#include <stdexcept>
struct Stop:runtime_error{using runtime_error::runtime_error;};
chrono::steady_clock::time_point epoch;ofstream operations,progress;
Mat inputmat(const string&path,int nr,int nc){ifstream f(path);int r,c;f>>r>>c;assert(r==nr&&c==nc);Mat A(r,vector<P>(c));for(auto&row:A)for(auto&p:row)readp(f,p);return A;}
double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-epoch).count();}
void guard(const vector<P>&row){long total=0;for(auto&p:row){long n=fmpq_mpoly_length(p.p,ctx);if(n>18000)throw Stop("single_entry18000terms");total+=n;}if(total>150000)throw Stop("single_row150000terms");if(elapsed()>53)throw Stop("native53s_boundary");}
void normalize(vector<P>&row,int id){vector<int>ids;for(int j=0;j<(int)row.size();j++)if(!fmpq_mpoly_is_zero(row[j].p,ctx))ids.push_back(j);if(ids.empty())return;sort(ids.begin(),ids.end(),[&](int x,int y){return fmpq_mpoly_length(row[x].p,ctx)<fmpq_mpoly_length(row[y].p,ctx);});P g=row[ids[0]],tmp;for(int j:ids){if(fmpq_mpoly_is_fmpq(g.p,ctx))break;assert(fmpq_mpoly_gcd(tmp.p,g.p,row[j].p,ctx));g=tmp;}if(!fmpq_mpoly_is_fmpq(g.p,ctx)){for(int j:ids){assert(fmpq_mpoly_divides(tmp.p,row[j].p,g.p,ctx));row[j]=tmp;}}else g=one();fmpq_t c;fmpq_init(c);fmpq_mpoly_get_term_coeff_fmpq(c,row[ids[0]].p,0,ctx);for(int j:ids)fmpq_mpoly_scalar_div_fmpq(row[j].p,row[j].p,c,ctx);operations<<"NORMALIZE "<<id<<'\n';writep(operations,g);char*s=fmpq_get_str(nullptr,10,c);operations<<s<<'\n';flint_free(s);fmpq_clear(c);guard(row);}
void eliminate(Mat&A,int target,int pivot,int column){if(fmpq_mpoly_is_zero(A[target][column].p,ctx))return;P g,u,v,x,y;assert(fmpq_mpoly_gcd(g.p,A[pivot][column].p,A[target][column].p,ctx));assert(fmpq_mpoly_divides(u.p,A[pivot][column].p,g.p,ctx));assert(fmpq_mpoly_divides(v.p,A[target][column].p,g.p,ctx));vector<P>row(A[target].size());for(int j=0;j<(int)row.size();j++){if(j==column)continue;fmpq_mpoly_mul(x.p,u.p,A[target][j].p,ctx);fmpq_mpoly_mul(y.p,v.p,A[pivot][j].p,ctx);fmpq_mpoly_sub(row[j].p,x.p,y.p,ctx);guard(row);}operations<<"ELIMINATE "<<target<<' '<<pivot<<' '<<column<<'\n';writep(operations,u);writep(operations,v);A[target]=move(row);normalize(A[target],target);}
pair<int,int>choose(const Mat&A,const vector<int>&rows,const vector<int>&cols){long best=LONG_MAX;pair<int,int>p{-1,-1};for(int r:rows)for(int c:cols){long n=fmpq_mpoly_length(A[r][c].p,ctx);if(n&&n<best){best=n;p={r,c};}}assert(p.first>=0);return p;}
void removeid(vector<int>&v,int i){auto p=find(v.begin(),v.end(),i);assert(p!=v.end());v.erase(p);}
void checkpoint(const string&out,const string&stage,const Mat&A,const vector<pair<int,int>>&piv,const vector<int>&remaining,int count){string tmp=out+"/current-matrix.tmp";writem(tmp,A);assert(!rename(tmp.c_str(),(out+"/current-matrix.txt").c_str()));ofstream f(out+"/current-stage.json");f<<"{\"stage\":\""<<stage<<"\",\"completed_stage_steps\":"<<count<<",\"seconds\":"<<elapsed()<<",\"pivots\":[";for(int i=0;i<(int)piv.size();i++)f<<(i?",":"")<<'['<<piv[i].first<<','<<piv[i].second<<']';f<<"],\"remaining_columns\":[";for(int i=0;i<(int)remaining.size();i++)f<<(i?",":"")<<remaining[i];f<<"]}\n";operations.flush();progress<<"CHECKPOINT "<<stage<<' '<<count<<" sec="<<elapsed()<<'\n';progress.flush();}
int main(int argc,char**argv){assert(argc==3);fmpq_mpoly_ctx_init(ctx,3,ORD_DEGREVLEX);{
 epoch=chrono::steady_clock::now();string source=argv[1],out=argv[2];operations.open(out+"/operations.txt");progress.open(out+"/progress.log");Mat Q=inputmat(source+"/quad10x26.txt",10,26),G=inputmat(source+"/genuine16x31.txt",16,31),A(26,vector<P>(31));for(int r=0;r<2;r++)for(int j=0;j<26;j++)A[r][j]=Q[8+r][j];for(int r=0;r<8;r++)for(int j=0;j<26;j++)A[r+2][j]=Q[r][j];for(int r=0;r<16;r++)A[r+10]=G[r];Q.clear();G.clear();
 vector<pair<int,int>>piv;vector<int>available;for(int j=0;j<31;j++)available.push_back(j);string stage="initial",status="blocked_no_identity_conclusion",stop;int completed=0;vector<int>heavy{11,12,18,20,21},remaining11{10,13,14,15,16,17,19,22,23,24,25};
 try{
  for(int r=0;r<26;r++)normalize(A[r],r);
  vector<int>qr{0,1,2,3,4,5,6,7,8,9};stage="quad10_echelon";
  for(int k=0;k<10;k++){vector<int>rr,cc;if(k<2){for(int r:qr)if(r<2)rr.push_back(r);for(int c:available)if(c>=20&&c<26)cc.push_back(c);}else{rr=qr;for(int c:available)if(c<20)cc.push_back(c);}auto[r,c]=choose(A,rr,cc);piv.push_back({r,c});removeid(qr,r);removeid(available,c);for(int t:qr)eliminate(A,t,r,c);completed=k+1;checkpoint(out,stage,A,piv,available,completed);}
  Mat qe;for(auto[r,c]:piv)qe.push_back(A[r]);writem(out+"/completed-quad10-echelon.txt",qe);
  stage="genuine16_reduce_quad10";for(int k=0;k<10;k++){auto[r,c]=piv[k];for(int t=10;t<26;t++)eliminate(A,t,r,c);completed=k+1;checkpoint(out,stage,A,piv,available,completed);}
  Mat qg;for(int r=10;r<26;r++){vector<P>v;for(int c:available)v.push_back(A[r][c]);qg.push_back(move(v));}writem(out+"/completed-afterquad16x21.txt",qg);
  stage="heavy5";vector<int>hr=heavy,gr;for(int r=10;r<26;r++)gr.push_back(r);
  for(int k=0;k<5;k++){vector<int>hc;for(int c:available)if(c>=26)hc.push_back(c);auto[r,c]=choose(A,hr,hc);piv.push_back({r,c});removeid(hr,r);removeid(gr,r);removeid(available,c);for(int t:gr)eliminate(A,t,r,c);completed=k+1;checkpoint(out,stage,A,piv,available,completed);}
  Mat r11;for(int r:remaining11){vector<P>v;for(int c:available)v.push_back(A[r][c]);r11.push_back(move(v));}assert(available.size()==16);writem(out+"/completed-R11x16.txt",r11);
  stage="ten_cubic_pivots_and_six_borders";vector<int>cr=remaining11;removeid(cr,16);
  for(int k=0;k<10;k++){auto[r,c]=choose(A,cr,available);piv.push_back({r,c});removeid(cr,r);removeid(available,c);for(int t:cr)eliminate(A,t,r,c);eliminate(A,16,r,c);completed=k+1;checkpoint(out,stage,A,piv,available,completed);if(terms(A[16])==0){status="all_six_zero_exact";break;}}
  if(status!="all_six_zero_exact"){assert(available.size()==6);status=terms(A[16])?"at_least_one_border_nonzero_exact":"all_six_zero_exact";}
  ofstream f(out+"/six-remainders.txt");f<<status<<'\n'<<available.size()<<'\n';for(int c:available){f<<c<<'\n';writep(f,A[16][c]);}
 }catch(const Stop&s){stop=s.what();progress<<"STOP "<<stop<<'\n';progress.flush();}
 ofstream f(out+"/native.json");f<<"{\"status\":\""<<status<<"\",\"last_stage\":\""<<stage<<"\",\"completed_stage_steps\":"<<completed<<",\"stop_reason\":\""<<stop<<"\",\"seconds\":"<<elapsed()<<",\"target_genuine_row\":6,\"source_rows_unchanged_global_ids\":true,\"only_generic_nonzero_factor_inversions\":true}\n";
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();return 0;}
