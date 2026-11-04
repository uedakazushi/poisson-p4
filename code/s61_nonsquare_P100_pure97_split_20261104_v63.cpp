// Saved-row arithmetic only. No source/oracle evaluation or Gaussian.
#define main frozen_q5_row_main_NOT_CALLED
#include "s61_nonsquare_P100_q5_pure98_row_20261104_v63.cpp"
#undef main
#include <zlib.h>
#include <streambuf>
long compressed_before=0;
struct GBuf:streambuf{
 gzFile f;char b[4096];
 explicit GBuf(const string&path){f=gzopen(path.c_str(),"wb9");assert(f);setp(b,b+sizeof b);}
 void flushbuf(){long n=pptr()-pbase();if(n){if(compressed_before+gzoffset(f)+131072>28*1024*1024)throw runtime_error("28MiB_compressed_output_guard");assert(gzwrite(f,b,n)==n);pbump(-n);}}
 int overflow(int c)override{flushbuf();if(c!=EOF){*pptr()=c;pbump(1);}return c;}
 int sync()override{flushbuf();return 0;}
 void close(){if(f){flushbuf();assert(gzflush(f,Z_FINISH)==Z_OK);compressed_before+=gzoffset(f);assert(gzclose(f)==Z_OK);f=nullptr;}}
 ~GBuf(){if(f)gzclose(f);}
};
bool certifiedunitden(const R&r){P p=r.d;for(int i=0;i<100;i++){P g;fmpq_poly_gcd(g.v,p.v,unitproduct.v);guard(g);if(fmpq_poly_degree(g.v)==0)break;p=quotient(p,g);}return fmpq_poly_degree(p.v)==0;}
int main(int argc,char**argv){assert(argc==3);began=chrono::steady_clock::now();string folder=argv[2];try{
 ifstream f(argv[1]);f>>input_reserve;vector<R>r(98);for(auto&v:r)v=inputr(f);R c0=inputr(f),c1=inputr(f);
 vector<vector<R>>V(7,vector<R>(98));for(auto&row:V)for(auto&v:row)v=inputr(f);string extra;assert(!(f>>extra));
 z=shiftp(one(),1);unitproduct=mulp(mulp(mulp(mulp(z,plusp(z,pc("20"))),plusp(z,pc("60"))),plusp(scalar(z,"203"),pc("2700"))),plusp(z,pc("12")));
 stage="literal_saved_q5_dual_and_GL2";assert(eq(add(mul(r[0],c0),mul(r[1],c1)),rc("1")));
 vector<vector<R>>T={{c0,neg(r[1])},{c1,r[0]}},Ti={{r[0],r[1]},{neg(c1),c0}};
 for(int i=0;i<2;i++)for(int j=0;j<2;j++){R a,b;for(int k=0;k<2;k++){a=add(a,mul(T[i][k],Ti[k][j]));b=add(b,mul(Ti[i][k],T[k][j]));}assert(eq(a,rc(i==j?"1":"0"))&&eq(b,rc(i==j?"1":"0")));}
 assert(eq(sub(mul(c0,r[0]),mul(neg(r[1]),c1)),rc("1")));
 R q0=add(mul(r[0],c0),mul(r[1],c1)),q1=add(mul(r[0],neg(r[1])),mul(r[1],r[0]));assert(eq(q0,rc("1"))&&zero(q1));
 for(int j=2;j<98;j++)assert(zero(sub(r[j],mul(q0,r[j]))));
 vector<int>nz,un;vector<long>rowterms;
 for(int i=0;i<7;i++){
  stage="saved_G18_row_"+to_string(i);auto&v=V[i];R a=add(mul(v[0],c0),mul(v[1],c1));vector<R>w;w.push_back(add(mul(v[0],neg(r[1])),mul(v[1],r[0])));
  for(int j=2;j<98;j++)w.push_back(sub(v[j],mul(a,r[j])));
  assert(eq(v[0],sub(mul(a,r[0]),mul(w[0],c1))));assert(eq(v[1],add(mul(a,r[1]),mul(w[0],c0))));
  for(int j=2;j<98;j++)assert(eq(v[j],add(w[j-1],mul(a,r[j]))));
  int n=0;long terms=0;for(auto&u:w){n+=!zero(u);assert(certifiedunitden(u));terms+=rterms(u);}nz.push_back(n);rowterms.push_back(terms);
  {GBuf gb(folder+"/pure97-row"+to_string(i)+".json.gz");ostream o(&gb);o<<"{\"source319_index\":"<<174+i<<",\"source_q5_cofactor\":";jr(o,a);o<<",\"full98_source_reconstruction_verified\":true,\"values97\":[";for(int j=0;j<97;j++){if(j)o<<',';jr(o,w[j]);}o<<"]}\n";o.flush();gb.close();}
  vector<R>().swap(V[i]);
 }
 ofstream o(folder+"/audit.json");o<<"{\"all_passed\":true,\"literal_GL2_det1_and_both_inverse_products_verified\":true,\"all98_q5_source_identity_verified\":true,\"all686_G18_source_reconstruction_identities_verified\":true,\"known_curve_denominators_only\":true,\"seven_nonzero_counts\":[";
 for(int i=0;i<7;i++){if(i)o<<',';o<<nz[i];}o<<"],\"seven_RF_termcounts\":[";for(int i=0;i<7;i++){if(i)o<<',';o<<rowterms[i];}o<<"],\"peak_live_RF_terms_with_inputreserve\":"<<peak_live+input_reserve<<",\"serialized_terms\":"<<stored<<",\"compressed_row_bytes\":"<<compressed_before<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";
 }catch(const exception&e){ofstream o(folder+"/frontier.json");o<<"{\"all_passed\":false,\"stage\":\""<<stage<<"\",\"reason\":\""<<e.what()<<"\"}\n";return 2;}flint_cleanup();}
