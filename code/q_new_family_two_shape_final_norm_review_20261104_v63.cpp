// ONE different final step from saved primitive P2/P1; no PRS/resultant rerun.
#define main unused_completed_prs_main
#include "q_new_family_two_shape_prs_review_20261104_v63.cpp"
#undef main
U read_array(const string&path){ifstream f(path);int nv,n;f>>nv>>n;assert(nv==5&&n>0);U a(n);fmpz_t c;fmpz_init(c);ulong e[5];for(auto&p:a){long k;f>>k;assert(k>0);for(long i=0;i<k;i++){string t;f>>t;for(auto&v:e)f>>v;assert(f&&fmpz_set_str(c,t.c_str(),10)==0);fmpz_mpoly_push_term_fmpz_ui(p.a,c,e,ctx);}fmpz_mpoly_sort_terms(p.a,ctx);fmpz_mpoly_combine_like_terms(p.a,ctx);assert(terms(p)==k);}string extra;assert(!(f>>extra));fmpz_clear(c);return a;}
void stats(ostream&o,const U&a){long degs[5]={0,0,0,0,0},bits=0;for(auto&p:a){slong d[5];fmpz_mpoly_degrees_si(d,p.a,ctx);for(int i=0;i<5;i++)degs[i]=max(degs[i],long(d[i]));fmpz_t c;fmpz_init(c);for(long j=0;j<terms(p);j++){fmpz_mpoly_get_term_coeff_fmpz(c,p.a,j,ctx);bits=max(bits,long(fmpz_bits(c)));}fmpz_clear(c);}o<<"{\"s_degree\":"<<deg(a)<<",\"coefficient_term_counts_ascending\":[";for(unsigned i=0;i<a.size();i++){if(i)o<<',';o<<terms(a[i]);}o<<"],\"total_terms\":"<<allterms(a)<<",\"max_coefficient_bits\":"<<bits<<",\"coefficient_multidegree_reversedLEX\":[";for(int i=0;i<5;i++){if(i)o<<',';o<<degs[i];}o<<"]}";}
void checkpoint_poly(const string&name,const P&p){
 ofstream o(output_dir+"/"+name+"-statistics.json");U a{p};o<<"{\"sha256_reversedLEX_terms\":\""<<digest(p)<<"\",\"statistics\":";stats(o,a);o<<"}\n";
 front(name,terms(p));
}
int main(int argc,char**argv){assert(argc==5);output_dir=argv[4];fmpz_mpoly_ctx_init(ctx,5,ORD_LEX);int status=0;auto start=chrono::steady_clock::now();{
 try{
  // These are small unrelated univariate fixtures, not new source points.
  fixture(constants({1,3,-1,0,2,1}),constants({-1,2,1,1}),"monic");
  fixture(constants({1,-5,2,-1,3,2}),constants({2,4,-2,3}),"nonmonic");front("fixture_signs_rechecked");
  U a=read_array(argv[1]),b=read_array(argv[2]);assert(a.size()==3&&b.size()==2);
  assert(allterms(a)==15005&&allterms(b)==21012);front("primitive_inputs_bound");
  P norm;
  {P sq=b[0]*b[0];checkpoint_poly("b0_squared",sq);norm=a[2]*sq;}
  checkpoint_poly("norm_first_term",norm);
  {P mixed=b[0]*b[1];P term=a[1]*mixed;norm=norm-term;}
  checkpoint_poly("norm_first_two_terms",norm);
  {P sq=b[1]*b[1];P term=a[0]*sq;norm=norm+term;}
  checkpoint_poly("complete_primitive_norm",norm);
  P denominator=a[2]*a[2];checkpoint_poly("known_leading_square",denominator);
  P quotient;bool exact=fmpz_mpoly_divides(quotient.a,norm.a,denominator.a,ctx);
  if(!exact)throw runtime_error("exact division of complete Norm by a2^2 returned false");
  checkpoint_poly("exact_quotient_returned",quotient);
  {P recomposed=denominator*quotient;assert(equal(recomposed,norm));}
  checkpoint_poly("exact_division_remultiplied",quotient);saveu("primitive_final_quotient",U{quotient});
  // Clear large unused norm/denominator arrays before restoring content.
  norm=P();denominator=P();
  P r,R;fmpz_mpoly_gen(r.a,4,ctx);fmpz_mpoly_gen(R.a,3,ctx);P rm1=r-P(1),T=rm1*R-P(2)*r;
  P scalar;fmpz_mpoly_set_ui(scalar.a,1UL<<52,ctx);
  P content=scalar*pw(r,18)*pw(rm1,10)*pw(T,13);
  P c0=content*quotient;checkpoint_poly("full_C0_from_known_contents",c0);
  P raw=load_result(argv[3]);assert(equal(neg(c0),raw));front("all47114_raw_resultant_coefficients_match",terms(raw));saveu("full_C0",U{c0});
  ofstream o(output_dir+"/native-report.json");
  o<<"{\"all_passed\":true,\"different_last_step_only_saved_primitive_P2_P1\":true,\"three_product_norm_identity\":\"a2*b0^2-a1*b0*b1+a0*b1^2\",\"complete_norm_exact_division_by_a2_squared_and_remultiplication\":true,\"restored_content\":\"2^52*r^18*(r-1)^10*T^13\",\"raw_all47114_integer_coefficients_match_negative_C0\":true,\"raw_resultant_identity_independently_verified\":true,\"quotient_terms\":"<<terms(quotient)<<",\"C0_terms\":"<<terms(c0)<<",\"C0_reversedLEX_sha256\":\""<<digest(c0)<<"\",\"no_PRS_resultant_gcd_prime_point_source26_or_alternate_retry\":true,\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
 }catch(const exception&e){ofstream o(output_dir+"/native-failure.json");o<<"{\"all_passed\":false,\"error\":\""<<e.what()<<"\"}\n";cerr<<e.what()<<endl;status=2;}
 }fmpz_mpoly_ctx_clear(ctx);flint_cleanup();return status;}
