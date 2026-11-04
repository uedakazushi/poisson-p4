// Verification only: exact QQ[e] arithmetic and two monic quotient reductions.
// No inverse, gcd, row reduction, rank search, or author implementation.
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <array>
#include <chrono>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

struct P {
  fmpq_poly_t v;
  P() { fmpq_poly_init(v); }
  P(const P& p) { fmpq_poly_init(v); fmpq_poly_set(v,p.v); }
  P(P&& p) noexcept { fmpq_poly_init(v); fmpq_poly_swap(v,p.v); }
  P& operator=(const P& p) { if(this!=&p) fmpq_poly_set(v,p.v); return *this; }
  P& operator=(P&& p) noexcept { if(this!=&p) fmpq_poly_swap(v,p.v); return *this; }
  ~P() { fmpq_poly_clear(v); }
};
P one(){ P z; fmpq_poly_one(z.v); return z; }
P plus(const P& x,const P& y){P z;fmpq_poly_add(z.v,x.v,y.v);return z;}
P minus(const P& x){P z;fmpq_poly_neg(z.v,x.v);return z;}
P times(const P& x,const P& y){P z;fmpq_poly_mul(z.v,x.v,y.v);return z;}
bool eq(const P& x,const P& y){return fmpq_poly_equal(x.v,y.v);}
P readP(){
  int n;if(!(std::cin>>n))throw std::runtime_error("EOF");P p;fmpq_t q;fmpq_init(q);
  for(int k=0;k<n;k++){int e;std::string c;std::cin>>e>>c;
    if(e<0||fmpq_set_str(q,c.c_str(),10))throw std::runtime_error("bad coefficient");
    fmpq_canonicalise(q);fmpq_poly_set_coeff_fmpq(p.v,e,q);
  }fmpq_clear(q);return p;
}
void check(bool x,const char* message){if(!x)throw std::runtime_error(message);}
P f6,q0,q1;
P red(const P& p){P z;fmpq_poly_rem(z.v,p.v,f6.v);return z;}
P em(const P& p,const P& q){return red(times(p,q));}
using B=std::array<P,2>;
using A=std::array<B,4>;
B bOne(){B z;z[0]=one();return z;}
B ba(const B& x,const B& y){return {plus(x[0],y[0]),plus(x[1],y[1])};}
B bn(const B& x){return {minus(x[0]),minus(x[1])};}
B scale(const B& x,const P& y){return {em(x[0],y),em(x[1],y)};}
B bm(const B& x,const B& y){P v=em(x[1],y[1]);return {
  plus(em(x[0],y[0]),minus(em(q0,v))),
  plus(plus(em(x[0],y[1]),em(x[1],y[0])),minus(em(q1,v)))};}
bool beq(const B& x,const B& y){return eq(x[0],y[0])&&eq(x[1],y[1]);}
B readB(){return {readP(),readP()};}
A readA(){A x;for(auto& b:x)b=readB();return x;}
std::array<B,4> h;
A aOne(){A x;x[0]=bOne();return x;}
A aa(const A& x,const A& y){A z;for(int j=0;j<4;j++)z[j]=ba(x[j],y[j]);return z;}
A ab(const A& x,const B& y){A z;for(int j=0;j<4;j++)z[j]=bm(x[j],y);return z;}
A am(const A& x,const A& y){
  std::array<B,7> z;
  for(int i=0;i<4;i++)for(int j=0;j<4;j++)z[i+j]=ba(z[i+j],bm(x[i],y[j]));
  for(int j=6;j>=4;j--){B v=z[j];for(int i=0;i<4;i++)z[j-4+i]=ba(z[j-4+i],bn(bm(v,h[i])));}
  return {z[0],z[1],z[2],z[3]};
}
bool aeq(const A& x,const A& y){for(int j=0;j<4;j++)if(!beq(x[j],y[j]))return false;return true;}
struct Term{int a,e,d;P c;};
using Source=std::vector<Term>;
Source readSource(){int n;std::cin>>n;Source p;for(int j=0;j<n;j++){
  int a,e,d;std::string c;std::cin>>a>>e>>d>>c;P z;fmpq_t v;fmpq_init(v);
  check(!fmpq_set_str(v,c.c_str(),10),"bad source rational");fmpq_canonicalise(v);
  fmpq_poly_set_coeff_fmpq(z.v,0,v);fmpq_clear(v);check(a>=0&&e>=0&&d>=0,"negative source exponent");
  p.push_back({a,e,d,std::move(z)});
}return p;}
std::vector<P> ep;
std::vector<B> ap;
std::vector<A> dp;
A image(const Source& p){A z;for(const auto& t:p)z=aa(z,ab(dp.at(t.d),scale(ap.at(t.a),em(t.c,ep.at(t.e)))));return z;}

int main(){try{
  auto start=std::chrono::steady_clock::now();f6=readP();long field_degree=fmpq_poly_degree(f6.v);check(field_degree==9||field_degree==11,"degree9_or11");
  Source Q=readSource(),H=readSource(),U=readSource();std::vector<Source> ps;for(int i=0;i<22;i++)ps.push_back(readSource());
  P savedq0=readP(),savedq1=readP();B ai=readB();std::array<B,5> hc,hm;
  for(auto& b:hc)b=readB();for(auto& b:hm)b=readB();
  int inverses;std::cin>>inverses;check(inverses==9,"nine inverse certificates");
  P q0_inverse;
  for(int i=0;i<inverses;i++){
    P element=readP(),inverse=readP(),bezf=readP(),beze=readP();
    check(eq(plus(times(bezf,f6),times(beze,element)),one()),"full QQ[e] Bezout");
    check(eq(em(element,inverse),one()),"quotient inverse");
    if(i==0){check(eq(element,savedq0),"q0 inverse anchor");q0_inverse=inverse;}
  }
  int ma=2,me=1,md=4;
  for(const auto& p:ps)for(const auto& t:p){ma=std::max(ma,t.a);me=std::max(me,t.e);md=std::max(md,t.d);}
  for(const auto& p:{Q,H,U})for(const auto& t:p){ma=std::max(ma,t.a);me=std::max(me,t.e);md=std::max(md,t.d);}
  P e;fmpq_poly_set_coeff_si(e.v,1,1);ep={one()};for(int j=0;j<me;j++)ep.push_back(em(ep.back(),e));
  P q2;
  for(const auto& t:Q){check(t.d==0&&t.a<=2,"source Q shape");P v=em(t.c,ep[t.e]);if(t.a==0)q0=plus(q0,v);else if(t.a==1)q1=plus(q1,v);else q2=plus(q2,v);}
  check(eq(q2,one())&&eq(q0,savedq0)&&eq(q1,savedq1),"actual Q coefficient binding");
  B av;av[1]=one();check(beq(bm(av,ai),bOne()),"a inverse");
  check(eq(ai[0],minus(em(q1,q0_inverse)))&&eq(ai[1],minus(q0_inverse)),"a inverse formula");
  ap={bOne()};for(int j=0;j<ma;j++)ap.push_back(bm(ap.back(),av));
  std::array<B,5> actualhc;
  for(const auto& t:H){check(t.d<=4,"H degree4");actualhc[t.d]=ba(actualhc[t.d],scale(ap[t.a],em(t.c,ep[t.e])));}
  B ai5=bOne();for(int j=0;j<5;j++)ai5=bm(ai5,ai);
  for(int j=0;j<5;j++){check(beq(actualhc[j],hc[j]),"original Hclear coefficient");check(beq(bm(hc[j],ai5),hm[j]),"H monic division by actual a5");if(j<4)h[j]=hm[j];}
  check(beq(hm[4],bOne()),"H monic top coefficient");
  A dd;dd[1]=bOne();dp={aOne()};for(int j=0;j<md;j++)dp.push_back(am(dp.back(),dd));
  check(aeq(image(Q),A())&&aeq(image(H),A()),"source quotient equations vanish");
  std::vector<A> images;
  for(int i=0;i<22;i++){A v=image(ps[i]);check(aeq(v,readA()),"original source image");images.push_back(std::move(v));}
  A ui=image(U);check(aeq(ui,readA()),"source U image");A target=aOne();for(int i=0;i<8;i++)target=am(target,ui);check(aeq(target,readA()),"U8 target");
  // Read saved matrix rowwise but compare every column to a reconstructed product.
  std::vector<A> columns(176);
  for(int row=0;row<8;row++)for(int col=0;col<176;col++)columns[col][row/2][row%2]=readP();
  for(int i=0;i<22;i++)for(int b=0;b<8;b++){A basis;basis[b/2][b%2]=one();check(aeq(am(images[i],basis),columns[8*i+b]),"all176 source ideal columns");}
  A dot;
  for(int i=0;i<22;i++)dot=aa(dot,am(readA(),images[i]));
  check(aeq(dot,target),"literal22 cofactor identity U8");std::string extra;check(!(std::cin>>extra),"trailing input");
  std::cout<<"{\"all_passed\":true,\"QQe_Bezout_identities\":9,\"actual_source_images\":22,\"actual_ideal_columns\":176,\"field_degree\":"<<field_degree<<",\"rational_matrix_slots\":"<<8*176*field_degree<<",\"literal_original22_cofactor_dot_equals_U8\":true,\"source_Q_H_U_checked\":true,\"no_inverse_GCD_Gaussian_rank_search\":true,\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
  return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
