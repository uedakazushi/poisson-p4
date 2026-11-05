#include <flint/fmpz.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
static long owned=0,peak=0;static long long pairs=0;static int maxbits=0;static string phase="input";static int current=-1;
static chrono::steady_clock::time_point started;static double deadline=48;
struct Z{fmpz_t x;Z(){fmpz_init(x);}~Z(){fmpz_clear(x);}Z(const Z&)=delete;Z(Z&&z){fmpz_init(x);fmpz_swap(x,z.x);}};
static void ck(bool b,const string&s){if(!b)throw runtime_error(s);}
static void guard(){peak=max(peak,owned+3000);ck(owned+3000<=73322,"owned_numeric_metadata_guard");ck(chrono::duration<double>(chrono::steady_clock::now()-started).count()<deadline,"time_guard");}
static void bits(const fmpz_t x){maxbits=max(maxbits,(int)fmpz_bits(x));ck(fmpz_bits(x)<=434,"434bit_guard");}
static string zs(const fmpz_t x){char*p=fmpz_get_str(nullptr,10,x);string s(p);flint_free(p);return s;}
static uint32_t pack(int a,int e,int d){ck(a>=-32&&a<992&&e>=-32&&e<992&&d>=-32&&d<992,"packed_key_range");return ((a+32)<<20)|((e+32)<<10)|(d+32);}
static array<int,3> unpack(uint32_t k){return {(int)(k>>20)-32,(int)((k>>10)&1023)-32,(int)(k&1023)-32};}
struct Spec{int a0=0,a1=-1,e0=0,e1=-1,d0=0,d1=-1,lo=0,hi=-1;string den="1";};
struct P{Spec b;vector<size_t> off;vector<Z> c;Z den;long charge=0;
 explicit P(const Spec&s):b(s){
  int A=max(0,b.a1-b.a0+1),E=max(0,b.e1-b.e0+1);size_t rows=(size_t)A*E;off.resize(rows+1);size_t n=0;
  for(int a=b.a0;a<=b.a1;a++)for(int e=b.e0;e<=b.e1;e++){size_t r=(a-b.a0)*E+e-b.e0;off[r]=n;int l=max(b.d0,b.lo-a-e),h=min(b.d1,b.hi-a-e);n+=max(0,h-l+1);}off[rows]=n;
  charge=n+off.size()+17;owned+=charge;guard();c.resize(n);ck(fmpz_set_str(den.x,b.den.c_str(),10)==0&&fmpz_sgn(den.x)>0,"den_input");
 }
 ~P(){owned-=charge;}
 int index(int a,int e,int d)const{if(a<b.a0||a>b.a1||e<b.e0||e>b.e1)return -1;int l=max(b.d0,b.lo-a-e),h=min(b.d1,b.hi-a-e);if(d<l||d>h)return -1;return off[(a-b.a0)*(b.e1-b.e0+1)+e-b.e0]+d-l;}
 template<class F>void each(F f)const{for(int a=b.a0;a<=b.a1;a++)for(int e=b.e0;e<=b.e1;e++){int l=max(b.d0,b.lo-a-e),h=min(b.d1,b.hi-a-e);for(int d=l;d<=h;d++){int k=index(a,e,d);if(!fmpz_is_zero(c[k].x))f(a,e,d,c[k].x);}}}
};
struct Sparse{vector<uint32_t>keys;vector<Z>co;long charge=0;~Sparse(){owned-=charge;}};
struct Def{int id;char op;int left=-1,right=-1;Spec spec;unique_ptr<Sparse> atom;};
static vector<Def> defs;static vector<unique_ptr<P>> val;static vector<int> users;
static unique_ptr<Sparse> readsparse(istream&in){int n;in>>n;ck(n>=0&&n<=11728,"sparse_input_count");auto p=make_unique<Sparse>();p->charge=2*n+16;owned+=p->charge;guard();p->keys.reserve(n);p->co.resize(n);for(int i=0;i<n;i++){int a,e,d;string z;in>>a>>e>>d>>z;p->keys.push_back(pack(a,e,d));ck(fmpz_set_str(p->co[i].x,z.c_str(),10)==0&&!fmpz_is_zero(p->co[i].x),"sparse_input_integer");bits(p->co[i].x);}return p;}
static unique_ptr<P> dense(const Spec&s,const Sparse&a){auto p=make_unique<P>(s);for(int i=0;i<(int)a.keys.size();i++){auto m=unpack(a.keys[i]);int k=p->index(m[0],m[1],m[2]);ck(k>=0&&fmpz_is_zero(p->c[k].x),"atom_offset_or_duplicate");fmpz_set(p->c[k].x,a.co[i].x);}return p;}
static void addto(P&z,const P&a,int sign){Z f,t;ck(fmpz_divisible(z.den.x,a.den.x),"addition_denominator_divexact");fmpz_divexact(f.x,z.den.x,a.den.x);a.each([&](int aa,int e,int d,const fmpz_t v){int k=z.index(aa,e,d);ck(k>=0,"addition_support");fmpz_mul(t.x,v,f.x);bits(t.x);if(sign>0)fmpz_add(z.c[k].x,z.c[k].x,t.x);else fmpz_sub(z.c[k].x,z.c[k].x,t.x);bits(z.c[k].x);});}
static unique_ptr<P> add(const Spec&s,const P*a,const P*b,int sign){auto z=make_unique<P>(s);if(a)addto(*z,*a,1);if(b)addto(*z,*b,sign);
 // Full backadd, coefficient by coefficient, with no polynomial copy.
 Z f,g,t,u,v;if(a){ck(fmpz_divisible(z->den.x,a->den.x),"backadd_denA");fmpz_divexact(f.x,z->den.x,a->den.x);}if(b){ck(fmpz_divisible(z->den.x,b->den.x),"backadd_denB");fmpz_divexact(g.x,z->den.x,b->den.x);}
 for(int aa=s.a0;aa<=s.a1;aa++)for(int e=s.e0;e<=s.e1;e++)for(int d=max(s.d0,s.lo-aa-e);d<=min(s.d1,s.hi-aa-e);d++){int k=z->index(aa,e,d);fmpz_zero(t.x);if(a){int j=a->index(aa,e,d);if(j>=0)fmpz_mul(t.x,a->c[j].x,f.x);}if(b){int j=b->index(aa,e,d);if(j>=0){fmpz_mul(u.x,b->c[j].x,g.x);if(sign>0)fmpz_add(t.x,t.x,u.x);else fmpz_sub(t.x,t.x,u.x);}}ck(fmpz_equal(t.x,z->c[k].x),"backadd_coefficient_identity");}return z;}
static void convolution(P&z,const P&a,const P&b,int sign){Z t; a.each([&](int aa,int ae,int ad,const fmpz_t av){b.each([&](int ba,int be,int bd,const fmpz_t bv){if((pairs&8191)==0)guard();int k=z.index(aa+ba,ae+be,ad+bd);ck(k>=0,"convolution_support");fmpz_mul(t.x,av,bv);bits(t.x);if(sign>0)fmpz_add(z.c[k].x,z.c[k].x,t.x);else fmpz_sub(z.c[k].x,z.c[k].x,t.x);bits(z.c[k].x);pairs++;});});}
static unique_ptr<P> mul(const Spec&s,const P&a,const P&b){auto z=make_unique<P>(s);Z D;fmpz_mul(D.x,a.den.x,b.den.x);ck(fmpz_equal(D.x,z->den.x),"multiplication_denominator_exact");convolution(*z,a,b,1);
 auto verify=make_unique<P>(s);for(size_t i=0;i<z->c.size();i++)fmpz_set(verify->c[i].x,z->c[i].x);convolution(*verify,b,a,-1);for(auto&v:verify->c)ck(fmpz_is_zero(v.x),"reverse_product_fullQQ_identity");return z;}
static Spec readspec(istream&in){Spec s;in>>s.a0>>s.a1>>s.e0>>s.e1>>s.d0>>s.d1>>s.lo>>s.hi>>s.den;return s;}
static void writepoly(const P&p,const string&path){ofstream o(path);o<<"{\"denominator\":\""<<zs(p.den.x)<<"\",\"integer_coefficients\":[";bool first=true;p.each([&](int a,int e,int d,const fmpz_t v){if(!first)o<<',';first=false;o<<"[["<<a<<','<<e<<','<<d<<"],\""<<zs(v)<<"\"]";});o<<"]}\n";ck(o.good(),"output_stream");}
int main(int argc,char**argv){if(argc==2&&string(argv[1])=="--help"){cout<<"CRAMER17D001 three fixed Cramer numerators,12 source polys,54 nodes\n";return 0;}if(argc!=3)return 2;started=chrono::steady_clock::now();string out=argv[2];
 try{ifstream in(argv[1]);string magic;int n;in>>magic>>n>>deadline;ck(magic=="CRAMER17D001"&&n==54&&deadline>0&&deadline<=48,"header");defs.resize(n);val.resize(n);users.assign(n,0);owned+=n*20;guard();
 for(int i=0;i<n;i++){auto&d=defs[i];in>>d.id>>d.op>>d.left>>d.right;d.spec=readspec(in);if(d.op=='i')d.atom=readsparse(in);else{ck(d.op=='a'||d.op=='s'||d.op=='m',"opcode");ck(d.left<i&&d.right<i&&d.left>=0&&d.right>=0,"source_order");users[d.left]++;users[d.right]++;}}
 array<int,3> roots;for(int&i:roots){in>>i;ck(i>=0&&i<n,"output_root");users[i]++;}string trailing;ck(!(in>>trailing),"trailing_input");
 for(int i=0;i<n;i++){current=i;phase="Cramer_fixed_DAG";guard();auto&d=defs[i];if(d.op=='i'){val[i]=dense(d.spec,*d.atom);d.atom.reset();}else if(d.op=='m')val[i]=mul(d.spec,*val[d.left],*val[d.right]);else val[i]=add(d.spec,val[d.left].get(),val[d.right].get(),d.op=='s'?-1:1);
 if(d.op!='i')for(int j:{d.left,d.right})if(--users[j]==0)val[j].reset();}
 array<long,3> count{};phase="stream_publication";for(int j=0;j<3;j++){val[roots[j]]->each([&](int,int,int,const fmpz_t v){ck(fmpz_bits(v)<=217,"final_coefficient_height");count[j]++;});writepoly(*val[roots[j]],out+"/Cramer"+to_string(j)+".json");}
 ofstream r(out+"/native.json");r<<"{\"all_passed\":true,\"all12_literal_inputs_and54_nodes\":true,\"all_product_reverseQQ_checks\":true,\"all_backadd_QQ_checks\":true,\"terms\":["<<count[0]<<','<<count[1]<<','<<count[2]<<"],\"integer_products\":"<<pairs<<",\"peak_owned_numeric_metadata\":"<<peak<<",\"observed_integer_bits\":"<<maxbits<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"}\n";cout<<"three numerators complete\n";return 0;
 }catch(const exception&e){ofstream r(out+"/native-failure.json");r<<"{\"all_passed\":false,\"reason\":\""<<e.what()<<"\",\"phase\":\""<<phase<<"\",\"cursor\":"<<current<<",\"no_result_claim\":true}";cerr<<e.what()<<'\n';return 1;}}
