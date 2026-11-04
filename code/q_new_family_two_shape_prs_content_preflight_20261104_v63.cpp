// LIGHT authorized preflight: strip only integer content and r,r-1,T,Z0.
// No pseudo-remainder, norm, resultant, or polynomial gcd is evaluated.
#define main unused_completed_prs_main
#include "q_new_family_two_shape_prs_review_20261104_v63.cpp"
#undef main
U read_array(const string&path){ifstream f(path);int nv,n;f>>nv>>n;assert(nv==5&&n>0);U a(n);fmpz_t c;fmpz_init(c);ulong e[5];for(auto&p:a){long k;f>>k;assert(k>0);for(long i=0;i<k;i++){string t;f>>t;for(auto&v:e)f>>v;assert(f&&fmpz_set_str(c,t.c_str(),10)==0);fmpz_mpoly_push_term_fmpz_ui(p.a,c,e,ctx);}fmpz_mpoly_sort_terms(p.a,ctx);fmpz_mpoly_combine_like_terms(p.a,ctx);assert(terms(p)==k);}string extra;assert(!(f>>extra));fmpz_clear(c);return a;}
string strz(const fmpz_t a){char*t=fmpz_get_str(nullptr,10,a);string s(t);flint_free(t);return s;}
void stats(ostream&o,const U&a){long degs[5]={0,0,0,0,0},bits=0;for(auto&p:a){slong d[5];fmpz_mpoly_degrees_si(d,p.a,ctx);for(int i=0;i<5;i++)degs[i]=max(degs[i],long(d[i]));fmpz_t c;fmpz_init(c);for(long j=0;j<terms(p);j++){fmpz_mpoly_get_term_coeff_fmpz(c,p.a,j,ctx);bits=max(bits,long(fmpz_bits(c)));}fmpz_clear(c);}o<<"{\"s_degree\":"<<deg(a)<<",\"coefficient_term_counts_ascending\":[";for(unsigned i=0;i<a.size();i++){if(i)o<<',';o<<terms(a[i]);}o<<"],\"total_terms\":"<<allterms(a)<<",\"max_coefficient_bits\":"<<bits<<",\"coefficient_multidegree_reversedLEX\":[";for(int i=0;i<5;i++){if(i)o<<',';o<<degs[i];}o<<"]}";}
void strip(const string&name,U a,const vector<P>&fac,const vector<string>&names,ostream&report){
 U original=a;fmpz_t content,c;fmpz_init(content);fmpz_init(c);fmpz_zero(content);
 for(auto&p:a)for(long j=0;j<terms(p);j++){fmpz_mpoly_get_term_coeff_fmpz(c,p.a,j,ctx);fmpz_gcd(content,content,c);}
 assert(!fmpz_is_zero(content));P F;fmpz_mpoly_set_fmpz(F.a,content,ctx);a=divexact(a,F);
 vector<int> powers;
 for(unsigned k=0;k<fac.size();k++){
  int n=0;
  while(true){U q(a.size());bool ok=true;for(unsigned i=0;i<a.size();i++)if(!fmpz_mpoly_divides(q[i].a,a[i].a,fac[k].a,ctx)){ok=false;break;}
   if(!ok)break;for(unsigned i=0;i<a.size();i++)assert(equal(q[i]*fac[k],a[i]));a=move(q);F=F*fac[k];++n;assert(n<500);
  }
  powers.push_back(n);front(name+"_stripped_"+names[k],allterms(a));
 }
 for(unsigned i=0;i<a.size();i++)assert(equal(F*a[i],original[i]));
 saveu(name+"-primitive",a);saveu(name+"-known-factor",U{F});
 report<<"\""<<name<<"\":{\"integer_content\":\""<<strz(content)<<"\",\"known_factor_powers\":{";
 for(unsigned i=0;i<names.size();i++){if(i)report<<',';report<<'"'<<names[i]<<"\":"<<powers[i];}
 report<<"},\"all_coefficient_global_recomposition_identities_passed\":true,\"no_further_COMMON_known_factor_divides_all_coefficients\":true,\"primitive_statistics\":";stats(report,a);report<<"}";
 fmpz_clear(c);fmpz_clear(content);
}
int main(int argc,char**argv){assert(argc==4);output_dir=argv[3];fmpz_mpoly_ctx_init(ctx,5,ORD_LEX);int status=0;auto start=chrono::steady_clock::now();{
 try{P r,R,Z;fmpz_mpoly_gen(r.a,4,ctx);fmpz_mpoly_gen(R.a,3,ctx);fmpz_mpoly_gen(Z.a,2,ctx);P rm1=r-P(1),T=rm1*R-P(2)*r;
  vector<P> fac={r,Z,rm1,T};vector<string> names={"r","Z0","r-1","T"};
  U c2=read_array(argv[1]),c1=read_array(argv[2]);assert(deg(c2)==2&&deg(c1)==1);ofstream o(output_dir+"/native-report.json");o<<"{\"all_passed\":true,\"coefficient_order\":[\"Z2\",\"Z1\",\"Z0\",\"R\",\"r\"],";
  strip("C2",move(c2),fac,names,o);o<<',';strip("C1",move(c1),fac,names,o);
  o<<",\"only_integer_content_and_named_factors_removed\":true,\"no_norm_PRS_resultant_polynomial_gcd_points_primes_or_source26\":true,\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
  front("both_known_content_strips_and_global_identities_complete");
 }catch(const exception&e){cerr<<e.what()<<endl;status=2;}
 }fmpz_mpoly_ctx_clear(ctx);flint_cleanup();return status;}
