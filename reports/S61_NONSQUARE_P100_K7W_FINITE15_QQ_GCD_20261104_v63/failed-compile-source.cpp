// One QQ gcd/xgcd of the SAVED H15 and R14. No determinant/source recomputation.
#include <flint/fmpq_poly.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
using namespace std;
auto began=chrono::steady_clock::now();
long live_slots=0,peak_slots=0,stored_slots=0,max_bits=0;
string stage="input";
struct P {
 fmpq_poly_t p; mutable long slots=0;
 P(){fmpq_poly_init(p);} P(const P&)=delete; P&operator=(const P&)=delete;
 ~P(){live_slots-=slots;fmpq_poly_clear(p);}
 void check()const{
  long n=fmpq_poly_length(p);live_slots+=n-slots;slots=n;peak_slots=max(peak_slots,live_slots);
  if(peak_slots+31>150000)throw runtime_error("150000_live_terms_guard");
  if(chrono::duration<double>(chrono::steady_clock::now()-began).count()>49)throw runtime_error("49s_native_guard");
  fmpq_t q;fmpq_init(q);
  for(long i=0;i<n;i++){fmpq_poly_get_coeff_fmpq(q,p,i);max_bits=max(max_bits,long(max(fmpz_bits(fmpq_numref(q)),fmpz_bits(fmpq_denref(q)))));}
  fmpq_clear(q);if(max_bits>300000)throw runtime_error("300000_coefficient_bits_guard");
 }
};
void readp(istream&in,P&v){long n;in>>n;assert(n>=0&&n<=16);fmpq_t q;fmpq_init(q);for(long j=0;j<n;j++){long e;string s;in>>e>>s;assert(e>=0&&e<=15);assert(!fmpq_set_str(q,s.c_str(),10));fmpq_canonicalise(q);fmpq_poly_set_coeff_fmpq(v.p,e,q);}fmpq_clear(q);v.check();}
string qs(const fmpq_t q){char*s=fmpq_get_str(nullptr,10,q);string t(s);flint_free(s);return t;}
void primitive(const P&v,P&w,fmpq_t content){
 fmpz_t den,g,num,t;fmpz_init(den);fmpz_one(den);fmpz_init(g);fmpz_init(num);fmpz_init(t);fmpq_t q;fmpq_init(q);
 for(long i=0;i<fmpq_poly_length(v.p);i++){fmpq_poly_get_coeff_fmpq(q,v.p,i);fmpz_lcm(den,den,fmpq_denref(q));}
 for(long i=0;i<fmpq_poly_length(v.p);i++){fmpq_poly_get_coeff_fmpq(q,v.p,i);fmpz_divexact(t,den,fmpq_denref(q));fmpz_mul(num,t,fmpq_numref(q));fmpz_gcd(g,g,num);}
 assert(!fmpz_is_zero(g));fmpq_poly_get_coeff_fmpq(q,v.p,fmpq_poly_degree(v.p));if(fmpq_sgn(q)<0)fmpz_neg(g,g);
 fmpq_set_fmpz_frac(content,g,den);fmpq_canonicalise(content);fmpq_poly_scalar_div_fmpq(w.p,v.p,content);w.check();
 P back;fmpq_poly_scalar_mul_fmpq(back.p,w.p,content);back.check();assert(fmpq_poly_equal(back.p,v.p));
 fmpq_clear(q);fmpz_clear(den);fmpz_clear(g);fmpz_clear(num);fmpz_clear(t);
}
void exactdiv(const P&a,const P&b,P&q){P rem,back;fmpq_poly_divrem(q.p,rem.p,a.p,b.p);q.check();rem.check();assert(fmpq_poly_is_zero(rem.p));fmpq_poly_mul(back.p,q.p,b.p);back.check();assert(fmpq_poly_equal(back.p,a.p));}
void jp(ostream&o,const P&v){stored_slots+=fmpq_poly_length(v.p);if(stored_slots>150000)throw runtime_error("150000_serialized_terms_guard");o<<'[';bool first=true;fmpq_t q;fmpq_init(q);for(long i=0;i<fmpq_poly_length(v.p);i++){fmpq_poly_get_coeff_fmpq(q,v.p,i);if(fmpq_is_zero(q))continue;if(!first)o<<',';first=false;o<<"[["<<i<<"],\""<<qs(q)<<"\"]";}fmpq_clear(q);o<<']';}
struct CapBuf:streambuf{string data;streamsize xsputn(const char*s,streamsize n)override{if(data.size()+n>27*1024*1024)throw runtime_error("27MiB_serialization_guard");data.append(s,n);return n;}int overflow(int c)override{if(c!=EOF){char x=c;xsputn(&x,1);}return c;}};
int main(int argc,char**argv){assert(argc==3);string folder=argv[2];try{
 ifstream in(argv[1]);P H,R,PH,PR,G,S,T,QH,QR;readp(in,H);readp(in,R);string extra;assert(!(in>>extra));assert(fmpq_poly_degree(H.p)==15&&fmpq_poly_degree(R.p)==14);
 fmpq_t ch,cr;fmpq_init(ch);fmpq_init(cr);primitive(H,PH,ch);primitive(R,PR,cr);
 stage="one_exact_QQ_xgcd";fmpq_poly_xgcd(G.p,S.p,T.p,PH.p,PR.p);G.check();S.check();T.check();assert(!fmpq_poly_is_zero(G.p));
 stage="literal_Bezout_and_divisibility_checks";P a,b,c;fmpq_poly_mul(a.p,S.p,PH.p);a.check();fmpq_poly_mul(b.p,T.p,PR.p);b.check();fmpq_poly_add(c.p,a.p,b.p);c.check();assert(fmpq_poly_equal(c.p,G.p));exactdiv(H,G,QH);exactdiv(R,G,QR);
 fmpq_t lc;fmpq_init(lc);fmpq_poly_get_coeff_fmpq(lc,G.p,fmpq_poly_degree(G.p));assert(fmpq_is_one(lc));long gd=fmpq_poly_degree(G.p);assert(gd<=2);
 stage="bounded_serialization";CapBuf buf;ostream o(&buf);o.exceptions(ios::badbit|ios::failbit);
 o<<"{\"all_passed\":true,\"one_saved_H15_R14_QQ_xgcd_only\":true,\"exact_QQ_gcd_degree\":"<<gd<<",\"whole_P100_pointwise_quartic_lower150\":"<<(gd==0?"true":"false")<<",\"original_H15\":";jp(o,H);o<<",\"original_saved_R14\":";jp(o,R);o<<",\"primitive_ZZ_H\":";jp(o,PH);o<<",\"primitive_ZZ_R\":";jp(o,PR);
 o<<",\"H_equals_content_times_primitive_content\":\""<<qs(ch)<<"\",\"R_equals_content_times_primitive_content\":\""<<qs(cr)<<"\",\"exact_gcd_monic\":";jp(o,G);o<<",\"QQ_Bezout_primitiveH_cofactor\":";jp(o,S);o<<",\"QQ_Bezout_primitiveR_cofactor\":";jp(o,T);o<<",\"H_div_gcd\":";jp(o,QH);o<<",\"R_div_gcd\":";jp(o,QR);
 o<<",\"both_content_source_identities_verified\":true,\"literal_g_equals_S_primitiveH_plus_T_primitiveR_verified\":true,\"both_literal_gcd_times_quotient_equals_saved_input_verified\":true,\"exact_gcd_certificate_not_modular_inference\":true,\"no_source_NF_minor_factor_point_prime_job\":true,\"maximum_observed_coefficient_bits\":"<<max_bits<<",\"peak_live_polynomial_slots_including31_input\":"<<peak_slots+31<<",\"serialized_polynomial_slots\":"<<stored_slots<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";
 ofstream out(folder+"/native.json");out<<buf.data;fmpq_clear(lc);fmpq_clear(ch);fmpq_clear(cr);
 }catch(const exception&e){ofstream out(folder+"/frontier.json");out<<"{\"all_passed\":false,\"stage\":\""<<stage<<"\",\"reason\":\""<<e.what()<<"\",\"max_observed_bits\":"<<max_bits<<",\"peak_live_slots\":"<<peak_slots+31<<"}\n";return 2;}flint_cleanup();}
