// GCD-only of the FIVE saved primitive coefficient polynomials; no QQ Bezout.
#define main frozen_q5_row_main_NOT_CALLED
#include "s61_nonsquare_P100_q5_pure98_row_20261104_v63.cpp"
#undef main
#include <flint/nmod_poly.h>
#include <flint/ulong_extras.h>
#include <sstream>
const ulong PRIME=32003;
P readpoly(istream&f){long n;f>>n;P p;fmpq_t q;fmpq_init(q);for(long i=0;i<n;i++){long e;string c;f>>e>>c;assert(!fmpq_set_str(q,c.c_str(),10));fmpq_canonicalise(q);fmpq_poly_set_coeff_fmpq(p.v,e,q);}fmpq_clear(q);guard(p);return p;}
bool same(const P&a,const P&b){return fmpq_poly_equal(a.v,b.v);}bool zz(const P&p){return fmpq_poly_is_zero(p.v);}P subp(const P&a,const P&b){return plusp(a,negp(b));}
bool pint(const P&p){bool ok=true;fmpq_t q;fmpq_init(q);for(long i=0;i<fmpq_poly_length(p.v);i++){fmpq_poly_get_coeff_fmpq(q,p.v,i);ok&=fmpz_fdiv_ui(fmpq_denref(q),PRIME)!=0;}fmpq_clear(q);return ok;}
bool lcunit(const P&p){if(zz(p))return true;fmpq_t q;fmpq_init(q);fmpq_poly_get_coeff_fmpq(q,p.v,fmpq_poly_degree(p.v));bool ok=fmpz_fdiv_ui(fmpq_denref(q),PRIME)!=0&&fmpz_fdiv_ui(fmpq_numref(q),PRIME)!=0;fmpq_clear(q);return ok;}
struct Primitive{P p;string c="1";};
Primitive primitive(const P&a){Primitive r;if(zz(a))return r;fmpz_t d,g,n,t;fmpz_init(d);fmpz_one(d);fmpz_init(g);fmpz_init(n);fmpz_init(t);fmpq_t q,c;fmpq_init(q);fmpq_init(c);
 for(long i=0;i<fmpq_poly_length(a.v);i++){fmpq_poly_get_coeff_fmpq(q,a.v,i);fmpz_lcm(d,d,fmpq_denref(q));}
 for(long i=0;i<fmpq_poly_length(a.v);i++){fmpq_poly_get_coeff_fmpq(q,a.v,i);fmpz_divexact(t,d,fmpq_denref(q));fmpz_mul(n,t,fmpq_numref(q));fmpz_gcd(g,g,n);}
 fmpq_poly_get_coeff_fmpq(q,a.v,fmpq_poly_degree(a.v));if(fmpq_sgn(q)<0)fmpz_neg(g,g);fmpq_set_fmpz_frac(c,g,d);fmpq_canonicalise(c);char*s=fmpq_get_str(nullptr,10,c);r.c=s;flint_free(s);fmpq_poly_scalar_div_fmpq(r.p.v,a.v,c);guard(r.p);assert(same(scalar(r.p,r.c),a));
 fmpq_clear(q);fmpq_clear(c);fmpz_clear(d);fmpz_clear(g);fmpz_clear(n);fmpz_clear(t);return r;}
struct NP{nmod_poly_t v;NP(){nmod_poly_init(v,PRIME);}NP(const NP&a):NP(){nmod_poly_set(v,a.v);}NP(NP&&a):NP(){nmod_poly_swap(v,a.v);}NP&operator=(NP a){nmod_poly_swap(v,a.v);return*this;}~NP(){nmod_poly_clear(v);}};
NP modp(const P&p){NP r;fmpq_t q;fmpq_init(q);for(long i=0;i<fmpq_poly_length(p.v);i++){fmpq_poly_get_coeff_fmpq(q,p.v,i);ulong d=fmpz_fdiv_ui(fmpq_denref(q),PRIME);assert(d);nmod_poly_set_coeff_ui(r.v,i,fmpz_fdiv_ui(fmpq_numref(q),PRIME)*n_invmod(d,PRIME)%PRIME);}fmpq_clear(q);return r;}
void jnp(ostream&o,const NP&p){o<<'[';bool first=true;for(long i=0;i<nmod_poly_length(p.v);i++){ulong c=nmod_poly_get_coeff_ui(p.v,i);if(c){if(!first)o<<',';first=false;o<<'['<<i<<','<<c<<']';}}o<<']';}
string leadstr(const P&p){fmpq_t q;fmpq_init(q);fmpq_poly_get_coeff_fmpq(q,p.v,fmpq_poly_degree(p.v));char*s=fmpq_get_str(nullptr,10,q);string out=s;flint_free(s);fmpq_clear(q);return out;}
struct PRS{int family,pair,step,exponent;P A,B,Q,R;string leading,content;};
vector<PRS>prs;
P prs_pair(P A,P B,int family){int step=0;if(fmpq_poly_degree(A.v)<fmpq_poly_degree(B.v))swap(A,B);
 while(!zz(B)){P R=A,Qpoly;string lc=leadstr(B);int e=0;
  while(!zz(R)&&fmpq_poly_degree(R.v)>=fmpq_poly_degree(B.v)){long shift=fmpq_poly_degree(R.v)-fmpq_poly_degree(B.v);P t=shiftp(pc(leadstr(R)),shift);Qpoly=plusp(scalar(Qpoly,lc),t);R=subp(scalar(R,lc),mulp(t,B));e++;}
  P scale=one();for(int k=0;k<e;k++)scale=scalar(scale,lc);assert(same(mulp(scale,A),plusp(mulp(Qpoly,B),R)));Primitive rr=primitive(R);
  prs.push_back({family,0,step++,e,A,B,Qpoly,rr.p,lc,rr.c});A=B;B=rr.p;
 }
 return primitive(A).p;}
struct CapBuf:streambuf{string data;streamsize xsputn(const char*s,streamsize n)override{if(data.size()+n>27*1024*1024)throw runtime_error("27MiB_JSON_guard");data.append(s,n);return n;}int overflow(int c)override{if(c!=EOF){char q=c;xsputn(&q,1);}return c;}};
int main(int argc,char**argv){assert(argc==3);began=chrono::steady_clock::now();string folder=argv[2];try{
 ifstream f(argv[1]);f>>input_reserve;vector<P>Z(5);for(auto&p:Z)p=readpoly(f);string extra;assert(!(f>>extra));stage="QQ_gcd_ONLY";P H;
 for(auto&p:Z){P g;fmpq_poly_gcd(g.v,H.v,p.v);guard(g);H=g;}assert(!zz(H)&&fmpq_poly_degree(H.v)<=15);Primitive hp=primitive(H);vector<P>Qv(5);
 bool gate=pint(H)&&lcunit(H);for(int i=0;i<5;i++){Qv[i]=quotient(Z[i],H);assert(same(mulp(H,Qv[i]),Z[i]));gate&=pint(Qv[i])&&lcunit(Qv[i]);}
 {ofstream o(folder+"/gcd-phase.json");o<<"{\"gcd_candidate_monic\":";jp(o,H);o<<",\"all5_exact_divisions\":[";for(int i=0;i<5;i++){if(i)o<<',';jp(o,Qv[i]);}o<<"],\"cofactor_oldprime_gate\":"<<(gate?"true":"false")<<",\"exact_gcd_certificate_pending\":true}\n";}
 stage="cofactor_mod32003_certificate";vector<NP>bars(5),weights(5);NP G;bool modunit=false;
 if(gate){for(int i=0;i<5;i++){bars[i]=modp(Qv[i]);assert(nmod_poly_degree(bars[i].v)==fmpq_poly_degree(Qv[i].v));}
  for(int i=0;i<5;i++){if(nmod_poly_is_zero(G.v)&&nmod_poly_is_zero(bars[i].v))continue;NP g,s,t,l,r,ch;nmod_poly_xgcd(g.v,s.v,t.v,G.v,bars[i].v);nmod_poly_mul(l.v,s.v,G.v);nmod_poly_mul(r.v,t.v,bars[i].v);nmod_poly_add(ch.v,l.v,r.v);assert(nmod_poly_equal(ch.v,g.v));
   for(int j=0;j<i;j++){NP v;nmod_poly_mul(v.v,weights[j].v,s.v);weights[j]=move(v);}weights[i]=t;G=g;}
  NP ch;for(int i=0;i<5;i++){NP v;nmod_poly_mul(v.v,weights[i].v,bars[i].v);nmod_poly_add(ch.v,ch.v,v.v);}assert(nmod_poly_equal(ch.v,G.v));modunit=nmod_poly_is_one(G.v);
 }
 bool prspass=false;if(!modunit){stage="small_integer_primitive_PRS_fallback";P g;for(int i=0;i<5;i++){P q=primitive(Qv[i]).p;if(zz(g))g=q;else g=prs_pair(g,q,i);}assert(fmpq_poly_degree(g.v)==0);prspass=true;}
 CapBuf buf;ostream o(&buf);o.exceptions(ios::badbit|ios::failbit);o<<"{\"all_passed\":true,\"exact_gcd_monic\":";jp(o,H);o<<",\"primitive_ZZ_gcd\":";jp(o,hp.p);o<<",\"monic_equals_QQscalar_times_primitive\":\""<<hp.c<<"\",\"all5_quotients\":[";
 for(int i=0;i<5;i++){if(i)o<<',';jp(o,Qv[i]);}o<<"],\"all5_literal_gcd_times_quotient_equals_input_verified\":true,\"cofactor_QQ_den_and_LC_primeunit_gate\":"<<(gate?"true":"false")<<",\"cofactor_reduction_degrees_preserved\":"<<(gate?"true":"false")<<",\"oldprime_cofactor_gcd_is1\":"<<(modunit?"true":"false")<<",\"modular_cofactor_gcd\":";jnp(o,G);o<<",\"cofactor_modular_Bezout5\":[";
 for(int i=0;i<5;i++){if(i)o<<',';jnp(o,weights[i]);}o<<"],\"cofactor_reductions5\":[";for(int i=0;i<5;i++){if(i)o<<',';jnp(o,bars[i]);}o<<"],\"literal_modular_Bezout_verified\":"<<(gate?"true":"false")<<",\"primitive_PRS_fallback_used\":"<<(prspass?"true":"false")<<",\"primitive_PRS_transcript\":[";
 for(int i=0;i<int(prs.size());i++){if(i)o<<',';auto&s=prs[i];o<<"{\"family\":"<<s.family<<",\"step\":"<<s.step<<",\"left_A\":";jp(o,s.A);o<<",\"divisor_B\":";jp(o,s.B);o<<",\"leading_scalar\":\""<<s.leading<<"\",\"exponent\":"<<s.exponent<<",\"pseudoquotient\":";jp(o,s.Q);o<<",\"primitive_remainder\":";jp(o,s.R);o<<",\"remainder_content\":\""<<s.content<<"\"}";}
 o<<"],\"exact_gcd_justified_by_divisibility_and_coprime_cofactors\":true,\"gcd_degree\":"<<fmpq_poly_degree(H.v)<<",\"whole_P100_lower150\":"<<(fmpq_poly_degree(H.v)==0?"true":"false")<<",\"no_QQ_Bezout_or_factorsearch\":true,\"peak_live_polynomial_terms_with_literal_input\":"<<peak_live+input_reserve<<",\"serialized_terms\":"<<stored<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";
 ofstream of(folder+"/native.json");of<<buf.data;
 }catch(const exception&e){ofstream o(folder+"/frontier.json");o<<"{\"all_passed\":false,\"stage\":\""<<stage<<"\",\"reason\":\""<<e.what()<<"\"}\n";return 2;}flint_cleanup();}
