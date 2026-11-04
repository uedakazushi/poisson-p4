// Fixed independent source substitution and fixed-order determinant only.
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
struct P{fmpq_poly_t v;P(){fmpq_poly_init(v);}P(const P&p):P(){fmpq_poly_set(v,p.v);}P(P&&p)noexcept:P(){fmpq_poly_swap(v,p.v);}P&operator=(const P&p){fmpq_poly_set(v,p.v);return *this;}P&operator=(P&&p)noexcept{fmpq_poly_swap(v,p.v);return *this;}~P(){fmpq_poly_clear(v);}};
struct Z{fmpz_poly_t v;Z(){fmpz_poly_init(v);}Z(const Z&p):Z(){fmpz_poly_set(v,p.v);}Z(Z&&p)noexcept:Z(){fmpz_poly_swap(v,p.v);}Z&operator=(const Z&p){fmpz_poly_set(v,p.v);return *this;}Z&operator=(Z&&p)noexcept{fmpz_poly_swap(v,p.v);return *this;}~Z(){fmpz_poly_clear(v);}};
void check(bool b,const char*s){if(!b)throw std::runtime_error(s);}
P one(){P z;fmpq_poly_one(z.v);return z;}
P mul(const P&a,const P&b){P z;fmpq_poly_mul(z.v,a.v,b.v);return z;}
P rational(const std::string&s){P z;fmpq_t q;fmpq_init(q);check(!fmpq_set_str(q,s.c_str(),10),"rational input");fmpq_canonicalise(q);fmpq_poly_set_coeff_fmpq(z.v,0,q);fmpq_clear(q);return z;}
P readP(){int n;check(bool(std::cin>>n),"EOF");P p;for(int j=0;j<n;j++){int e;std::string v;std::cin>>e>>v;P c=rational(v);fmpq_t q;fmpq_init(q);fmpq_poly_get_coeff_fmpq(q,c.v,0);fmpq_poly_set_coeff_fmpq(p.v,e,q);fmpq_clear(q);}return p;}
Z integer(const P&p){check(fmpz_is_one(fmpq_poly_denref(p.v)),"integer expected");Z z;fmpq_poly_get_numerator(z.v,p.v);return z;}
bool eq(const P&a,const P&b){return fmpq_poly_equal(a.v,b.v);}
struct Term{int a,e,d;P c;};
using Source=std::vector<Term>;
Source readSource(){int n;std::cin>>n;Source p;for(int j=0;j<n;j++){int a,e,d;std::string c;std::cin>>a>>e>>d>>c;p.push_back({a,e,d,rational(c)});}return p;}
std::vector<P>ap,ep,dp;
std::map<std::tuple<int,int,int>,P>cache;
P pull(const Source&s){P result;for(const auto&t:s){int shift=24-2*t.a-t.e+t.d;check(shift>=0,"source Laurent denominator bound");auto key=std::make_tuple(t.a,t.e,t.d);auto it=cache.find(key);if(it==cache.end())it=cache.emplace(key,mul(mul(ap.at(t.a),ep.at(t.e)),dp.at(t.d))).first;P v=mul(t.c,it->second);fmpq_poly_shift_left(v.v,v.v,shift);fmpq_poly_add(result.v,result.v,v.v);}return result;}
void exact_shift_right(P&p,int n){fmpq_t q;fmpq_init(q);for(int j=0;j<n;j++){fmpq_poly_get_coeff_fmpq(q,p.v,j);check(fmpq_is_zero(q),"Laurent clearing exactness");}fmpq_clear(q);fmpq_poly_shift_right(p.v,p.v,n);}
int main(){try{
 auto start=std::chrono::steady_clock::now();
 P an=readP(),en=readP(),dn=readP();ap={one()};ep={one()};dp={one()};
 for(int j=0;j<13;j++)ap.push_back(mul(ap.back(),an));for(int j=0;j<23;j++)ep.push_back(mul(ep.back(),en));for(int j=0;j<6;j++)dp.push_back(mul(dp.back(),dn));
 Source Bsrc=readSource();P b=pull(Bsrc);exact_shift_right(b,24);check(eq(b,readP()),"B literal substitution");
 std::vector<std::vector<P>> entries(17,std::vector<P>(24));long rawterms=0;
 for(int i=0;i<17;i++)for(int j=0;j<24;j++){
   int kind;std::string scale;std::cin>>kind>>scale;Source p=readSource();rawterms+=p.size();P v=pull(p);
   if(kind!=1){v=mul(v,an);exact_shift_right(v,2);}if(kind!=2)v=mul(v,b);
   check(kind>=0&&kind<=2,"denominator kind");P invscale=rational("1");fmpq_t qq;fmpq_init(qq);check(!fmpq_set_str(qq,scale.c_str(),10),"scale input");fmpq_canonicalise(qq);check(!fmpq_is_zero(qq),"nonzero denominator scale");fmpq_poly_scalar_div_fmpq(v.v,v.v,qq);fmpq_clear(qq);entries[i][j]=std::move(v);
 }
 check(rawterms==84253,"all original terms");
 std::vector<std::vector<Z>> matrix(17,std::vector<Z>(24));long terms=0,degree=0,bits=0;
 fmpz_t lcm,content,saved,c;fmpz_init(lcm);fmpz_init(content);fmpz_init(saved);fmpz_init(c);
 for(int i=0;i<17;i++){
   std::string ls,cs;std::cin>>ls>>cs;fmpz_one(lcm);fmpz_zero(content);
   for(int j=0;j<24;j++)fmpz_lcm(lcm,lcm,fmpq_poly_denref(entries[i][j].v));
   check(!fmpz_set_str(saved,ls.c_str(),10)&&fmpz_equal(lcm,saved),"literal row LCM");
   fmpq_t ql;fmpq_init(ql);fmpq_set_fmpz(ql,lcm);
   for(int j=0;j<24;j++){
     P w;fmpq_poly_scalar_mul_fmpq(w.v,entries[i][j].v,ql);matrix[i][j]=integer(w);fmpz_poly_content(c,matrix[i][j].v);fmpz_gcd(content,content,c);
   }fmpq_clear(ql);check(!fmpz_set_str(saved,cs.c_str(),10)&&fmpz_equal(content,saved)&&!fmpz_is_zero(content),"literal row content");
   for(int j=0;j<24;j++){
     fmpz_poly_scalar_divexact_fmpz(matrix[i][j].v,matrix[i][j].v,content);Z expected=integer(readP());check(fmpz_poly_equal(matrix[i][j].v,expected.v),"all408 exact source entries");
     degree=std::max(degree,fmpz_poly_degree(matrix[i][j].v));
     for(int k=0;k<fmpz_poly_length(matrix[i][j].v);k++){const fmpz* v=fmpz_poly_get_coeff_ptr(matrix[i][j].v,k);if(!fmpz_is_zero(v))terms++;bits=std::max(bits,(long)fmpz_bits(v));}
   }
 }
 check(terms==5182&&degree==34&&bits==169,"substituted matrix profile");
 std::vector<int> rows(17),cols(17);for(int&i:rows)std::cin>>i;for(int&i:cols)std::cin>>i;
 std::vector<std::vector<Z>> m(17,std::vector<Z>(17));for(int i=0;i<17;i++)for(int j=0;j<17;j++){check(rows[i]>=0&&rows[i]<17&&cols[j]>=0&&cols[j]<24,"minor indices");m[i][j]=matrix[rows[i]][cols[j]];}
 Z previous;fmpz_poly_one(previous.v);long divisions=0;
 for(int k=0;k<16;k++){
   Z pivot=m[k][k];check(!fmpz_poly_is_zero(pivot.v),"fixed chosen pivot nonzero; no search");
   for(int i=k+1;i<17;i++)for(int j=k+1;j<17;j++){
     Z numerator,tmp,quotient,remainder,replay;
     fmpz_poly_mul(numerator.v,pivot.v,m[i][j].v);fmpz_poly_mul(tmp.v,m[i][k].v,m[k][j].v);fmpz_poly_sub(numerator.v,numerator.v,tmp.v);
     fmpz_poly_divrem(quotient.v,remainder.v,numerator.v,previous.v);check(fmpz_poly_is_zero(remainder.v),"Bareiss exact division");
     fmpz_poly_mul(replay.v,quotient.v,previous.v);check(fmpz_poly_equal(replay.v,numerator.v),"every Bareiss division remultiplied");m[i][j]=std::move(quotient);divisions++;
   }previous=std::move(pivot);
 }
 Z expected=integer(readP()),fflu=integer(readP()),negative;
 check(!fmpz_poly_is_zero(m[16][16].v)&&fmpz_poly_equal(expected.v,m[16][16].v),"literal fixed17 determinant");
 fmpz_poly_neg(negative.v,fflu.v);check(fmpz_poly_equal(expected.v,fflu.v)||fmpz_poly_equal(expected.v,negative.v),"saved denominator sign");
 long detterms=0;for(int k=0;k<fmpz_poly_length(expected.v);k++)detterms+=!fmpz_is_zero(fmpz_poly_get_coeff_ptr(expected.v,k));
 check(fmpz_poly_degree(expected.v)==447&&detterms==152,"determinant support");std::string extra;check(!(std::cin>>extra),"trailing input");
 std::cout<<"{\"all_passed\":true,\"all408_source_entries\":true,\"raw_terms\":"<<rawterms<<",\"distinct_monomials\":"<<cache.size()<<",\"all17_row_clearings\":true,\"fixed17_determinant_exact\":true,\"minor_degree\":447,\"minor_terms\":152,\"Bareiss_exact_remultiplications\":"<<divisions<<",\"no_rank_search_points_primes_or_author_solver\":true,\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
 fmpz_clear(lcm);fmpz_clear(content);fmpz_clear(saved);fmpz_clear(c);return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
