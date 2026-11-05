#include <flint/fmpz.h>
#include <gmp.h>
#include <cstdint>
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
static long owned=0;static int DB=108,CB=1096,RB=2048,OWN=33275;
struct Z { fmpz_t x; Z(){fmpz_init(x);}~Z(){fmpz_clear(x);} Z(const Z&)=delete;Z&operator=(const Z&)=delete;Z(Z&&z)noexcept{fmpz_init(x);fmpz_swap(x,z.x);} };
struct P {int cap,d=-1;vector<Z>c;explicit P(int n):cap(n),c(n){owned+=n;}~P(){owned-=cap;}P(const P&)=delete;P&operator=(const P&)=delete;
 void zero(){for(auto &z:c)fmpz_zero(z.x);d=-1;}void trim(){d=cap-1;while(d>=0&&fmpz_is_zero(c[d].x))--d;}};
static vector<unique_ptr<P>> M,B;
static array<int,12> rowperm,colperm;
static Z H,coef;
static chrono::steady_clock::time_point beginTime;
static double secLimit=45;
static long long multiplications=0,scalarDivisions=0,updatedCells=0;
static int stage=-1,ci=-1,cj=-1,xi=-1,xj=-1,divshift=-1,divindex=-1,maxObserved=0;
static string phase="initial",pathPrefix;static int signDet=1;
static P&at(int i,int j){return *M[i*12+j];}
static double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-beginTime).count();}
static string zz(const fmpz_t x){char*p=fmpz_get_str(nullptr,10,x);string s(p);flint_free(p);return s;}
static void require(bool c,const string&s){if(!c)throw runtime_error(s);}
static void bitcheck(const fmpz_t x,int b){int v=fmpz_bits(x);maxObserved=max(maxObserved,v);require(v<=b,"integer_bit_guard");}
static void guard(){require(owned+15000+1494<=OWN,"owned_scalar_guard");if(elapsed()>secLimit)throw runtime_error("soft_time_guard");}
static void writepoly(ostream&o,const P&p){o<<'[';bool first=true;for(int g=0;g<p.cap;g++)if(!fmpz_is_zero(p.c[g].x)){if(!first)o<<',';first=false;o<<'['<<g<<",\""<<zz(p.c[g].x)<<"\"]";}o<<']';}
static Z product_tmp,sum_tmp;
static void addmul_checked(const fmpz_t xx,const fmpz_t yy,fmpz_t dst,int sign){
 require(fmpz_bits(xx)+fmpz_bits(yy)<=RB,"before_product_bit_guard");fmpz_mul(product_tmp.x,xx,yy);bitcheck(product_tmp.x,RB);
 require(max(fmpz_bits(dst),fmpz_bits(product_tmp.x))+1<=RB,"before_sum_bit_guard");
 if(sign>0)fmpz_add(sum_tmp.x,dst,product_tmp.x);else fmpz_sub(sum_tmp.x,dst,product_tmp.x);bitcheck(sum_tmp.x,RB);fmpz_set(dst,sum_tmp.x);}
static void copy(P&d,const P&s){require(d.cap>s.d,"copy_capacity");d.zero();for(int i=0;i<=s.d;i++)fmpz_set(d.c[i].x,s.c[i].x);d.d=s.d;}
static void mul(P&z,const P&a,const P&b,const string&label){z.zero();phase=label;xi=xj=-1;if(a.d<0||b.d<0)return;require(a.d+b.d<z.cap,"multiply_degree_guard");
 for(int i=0;i<=a.d;i++){if(fmpz_is_zero(a.c[i].x))continue;for(int j=0;j<=b.d;j++){if(fmpz_is_zero(b.c[j].x))continue;xi=i;xj=j;if((multiplications&4095)==0)guard();addmul_checked(a.c[i].x,b.c[j].x,z.c[i+j].x,1);bitcheck(z.c[i+j].x,RB);multiplications++;}}z.trim();xi=xj=-1;}
static void subtract(P&z,const P&a,const P&b){phase="subtract_products";z.zero();int n=max(a.d,b.d);require(n<z.cap,"subtract_capacity");for(int i=0;i<=n;i++){xi=i;guard();require(max(fmpz_bits(a.c[i].x),fmpz_bits(b.c[i].x))+1<=RB,"before_subtract_bit_guard");fmpz_sub(z.c[i].x,a.c[i].x,b.c[i].x);bitcheck(z.c[i].x,RB);}z.trim();xi=-1;}
static void exactdivide(P&q,P&r,const P&den){q.zero();phase="exact_division_setup";require(den.d>=0,"zero_previous_pivot");
 while(r.d>=den.d){divshift=r.d-den.d;divindex=-1;phase="exact_division_quotient";guard();require(divshift<=DB&&divshift<q.cap,"quotient_degree_guard");require(fmpz_divisible(r.c[r.d].x,den.c[den.d].x),"nonexact_integer_leading_division");
  fmpz_divexact(coef.x,r.c[r.d].x,den.c[den.d].x);scalarDivisions++;require(fmpz_cmpabs(coef.x,H.x)<=0,"quotient_minor_height_guard");require(fmpz_is_zero(q.c[divshift].x),"duplicate_quotient_degree");fmpz_set(q.c[divshift].x,coef.x);
  for(int j=0;j<=den.d;j++){divindex=j;phase="exact_division_subtract";if((multiplications&4095)==0)guard();addmul_checked(coef.x,den.c[j].x,r.c[j+divshift].x,-1);bitcheck(r.c[j+divshift].x,RB);multiplications++;}r.trim();
 }require(r.d<0,"nonzero_polynomial_remainder");q.trim();require(q.d<=DB,"minor_degree_guard");for(int i=0;i<=q.d;i++)require(fmpz_cmpabs(q.c[i].x,H.x)<=0,"minor_coefficient_height_guard");divshift=divindex=-1;}
static string input_sha;static long long published=0;static int serial_degree=-1;
static void u16(ostream&o,uint16_t v){o.write(reinterpret_cast<char*>(&v),2);}static void u32(ostream&o,uint32_t v){o.write(reinterpret_cast<char*>(&v),4);}static void u64(ostream&o,uint64_t v){o.write(reinterpret_cast<char*>(&v),8);}
static void putz(ostream&o,const fmpz_t z){mpz_t m;mpz_init(m);fmpz_get_mpz(m,z);uint8_t sign=mpz_sgn(m)<0;mpz_abs(m,m);unsigned char b[256];size_t n=0;require(mpz_sizeinbase(m,2)<=2048,"binary_bit_guard");mpz_export(b,&n,1,1,1,0,m);require(n<=256,"binary_length_guard");o.write(reinterpret_cast<char*>(&sign),1);u16(o,n);o.write(reinterpret_cast<char*>(b),n);mpz_clear(m);}
static void putpoly(ostream&o,const P&p){uint16_t n=0;for(int k=0;k<=p.d;k++)if(!fmpz_is_zero(p.c[k].x))n++;u16(o,n);for(int k=0;k<=p.d;k++)if(!fmpz_is_zero(p.c[k].x)){serial_degree=k;u16(o,k);putz(o,p.c[k].x);}serial_degree=-1;}
static void checkpoint(const string&why){ofstream o(pathPrefix+"/guard-M12AS001.bin",ios::binary);o.write("M12AS001",8);u16(o,input_sha.size());o<<input_sha;u16(o,why.size());o<<why;u16(o,phase.size());o<<phase;
 int v[]={stage,ci,cj,xi,xj,divshift,divindex,serial_degree,maxObserved};for(int x:v)u32(o,(uint32_t)x);u64(o,multiplications);u64(o,scalarDivisions);u64(o,updatedCells);putz(o,H.x);putz(o,coef.x);putz(o,product_tmp.x);putz(o,sum_tmp.x);u16(o,M.size());for(auto&p:M)putpoly(o,*p);u16(o,B.size());for(auto&p:B)putpoly(o,*p);o.close();}
int main(int argc,char**argv){
 if(argc==3&&string(argv[1])=="--binary-fixture"){beginTime=chrono::steady_clock::now();ofstream o(argv[2],ios::binary);o.write("M12FX001",8);Z z;fmpz_set_si(z.x,-1);putz(o,z.x);fmpz_zero(z.x);putz(o,z.x);fmpz_one(z.x);fmpz_mul_2exp(z.x,z.x,2047);fmpz_sub_ui(z.x,z.x,1);putz(o,z.x);P p(3);fmpz_set_si(p.c[0].x,-1);fmpz_set(p.c[2].x,z.x);p.d=2;putpoly(o,p);return 0;}
 if(argc==2&&string(argv[1])=="--help"){cout<<"M12AXIS001 input outputdir SHA256: fixed12 integer Bareiss, binary M12BP001/M12AS001, hard2048/33275\n";return 0;}
 if(argc!=4)return 2;pathPrefix=argv[2];input_sha=argv[3];beginTime=chrono::steady_clock::now();
 try{ifstream in(argv[1]);string magic,hs;int n;in>>magic>>n>>DB>>CB>>RB>>OWN>>secLimit;require(magic=="M12AXIS001"&&n==12&&DB==108&&CB==1096&&RB==2048&&OWN==33275&&secLimit>0&&secLimit<=45&&input_sha.size()==64,"header");in>>hs;require(fmpz_set_str(H.x,hs.c_str(),10)==0&&fmpz_bits(H.x)==1096,"source_H");
 M.reserve(144);for(int i=0;i<144;i++)M.emplace_back(make_unique<P>(109));B.reserve(5);for(int i=0;i<5;i++)B.emplace_back(make_unique<P>(217));
 for(int k=0;k<144;k++){int i,j,nt;in>>i>>j>>nt;require(i==k/12&&j==k%12&&nt>=0&&nt<=10,"matrix_input_index");int old=-1;for(int q=0;q<nt;q++){int d;string z;in>>d>>z;require(d>old&&d<=9&&d>=0,"matrix_degree");old=d;require(fmpz_set_str(at(i,j).c[d].x,z.c_str(),10)==0&&!fmpz_is_zero(at(i,j).c[d].x),"matrix_integer");require(fmpz_bits(at(i,j).c[d].x)<=118,"input_source_height");}at(i,j).trim();}
 string trailing;require(!(in>>trailing),"input_trailing");for(int i=0;i<12;i++)rowperm[i]=colperm[i]=i;fmpz_one(B[4]->c[0].x);B[4]->d=0;
 for(int k=0;k<12;k++){stage=k;ci=cj=-1;phase="fixed_source_pivot";guard();if(k<11)require(at(k,k).d>=0,"unexpected_zero_old11_fixed_pivot");
 ofstream patch(pathPrefix+"/stage-"+to_string(k)+"-M12BP001.bin",ios::binary);patch.write("M12BP001",8);u16(patch,k);u16(patch,12);u16(patch,input_sha.size());patch<<input_sha;putpoly(patch,at(k,k));u32(patch,(11-k)*(11-k));
 if(k<11){for(int i=k+1;i<12;i++)for(int j=k+1;j<12;j++){ci=i;cj=j;mul(*B[0],at(k,k),at(i,j),"pivot_product");mul(*B[1],at(i,k),at(k,j),"cross_product");subtract(*B[2],*B[0],*B[1]);exactdivide(*B[3],*B[2],*B[4]);phase="commit_exact_quotient";guard();copy(at(i,j),*B[3]);updatedCells++;u16(patch,i);u16(patch,j);putpoly(patch,*B[3]);}
 copy(*B[4],at(k,k));for(int i=k+1;i<12;i++)at(i,k).zero();}patch.close();}
 phase="complete";bool nonzero=at(11,11).d>=0;ofstream o(pathPrefix+"/native.json");o<<"{\"all_passed\":true,\"fixed_old11_pivots_only\":true,\"all_exact_ZZ_divisions_remainder_ZERO\":true,\"determinant_C12\":";writepoly(o,at(11,11));o<<",\"axis_rank\":"<<(nonzero?12:11)<<",\"axis_determinant_NONZERO\":"<<(nonzero?"true":"false")<<",\"all12_binary_witnesses_complete\":true,\"integer_products\":"<<multiplications<<",\"exact_integer_divisions\":"<<scalarDivisions<<",\"observed_integer_bits\":"<<maxObserved<<",\"owned_with_ALL_metadata_and_input\":"<<owned+15000+1494<<",\"seconds\":"<<elapsed()<<"}\n";cout<<"complete axis_rank="<<(nonzero?12:11)<<" products="<<multiplications<<'\n';return 0;
 }catch(const exception&e){if(M.size()==144&&B.size()==5)checkpoint(e.what());ofstream o(pathPrefix+"/native-failure.json");o<<"{\"all_passed\":false,\"reason\":\""<<e.what()<<"\",\"binary_state_saved\":true,\"no_determinant_result_claim\":true}";cerr<<e.what()<<'\n';return 1;}}
