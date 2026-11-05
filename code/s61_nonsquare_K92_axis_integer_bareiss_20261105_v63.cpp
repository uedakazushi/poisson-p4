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
static long owned=0;static int DB=88,CB=864,RB=1738,OWN=27812;
struct Z { fmpz_t x; Z(){fmpz_init(x);}~Z(){fmpz_clear(x);} Z(const Z&)=delete;Z&operator=(const Z&)=delete;Z(Z&&z)noexcept{fmpz_init(x);fmpz_swap(x,z.x);} };
struct P {int cap,d=-1;vector<Z>c;explicit P(int n):cap(n),c(n){owned+=n;}~P(){owned-=cap;}P(const P&)=delete;P&operator=(const P&)=delete;
 void zero(){for(auto &z:c)fmpz_zero(z.x);d=-1;}void trim(){d=cap-1;while(d>=0&&fmpz_is_zero(c[d].x))--d;}};
static vector<unique_ptr<P>> M,B;
static array<int,11> rowperm,colperm;
static Z H,coef;
static chrono::steady_clock::time_point beginTime;
static double secLimit=45;
static long long multiplications=0,scalarDivisions=0,updatedCells=0;
static int stage=-1,ci=-1,cj=-1,xi=-1,xj=-1,divshift=-1,divindex=-1,maxObserved=0;
static string phase="initial",pathPrefix;static int signDet=1;
static P&at(int i,int j){return *M[i*11+j];}
static double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-beginTime).count();}
static string zz(const fmpz_t x){char*p=fmpz_get_str(nullptr,10,x);string s(p);flint_free(p);return s;}
static void require(bool c,const string&s){if(!c)throw runtime_error(s);}
static void bitcheck(const fmpz_t x,int b){int v=fmpz_bits(x);maxObserved=max(maxObserved,v);require(v<=b,"integer_bit_guard");}
static void guard(){require(owned+15000+1158<=OWN,"owned_scalar_guard");if(elapsed()>secLimit)throw runtime_error("soft_time_guard");}
static void writepoly(ostream&o,const P&p){o<<'[';bool first=true;for(int g=0;g<p.cap;g++)if(!fmpz_is_zero(p.c[g].x)){if(!first)o<<',';first=false;o<<'['<<g<<",\""<<zz(p.c[g].x)<<"\"]";}o<<']';}
static void copy(P&d,const P&s){require(d.cap>s.d,"copy_capacity");d.zero();for(int i=0;i<=s.d;i++)fmpz_set(d.c[i].x,s.c[i].x);d.d=s.d;}
static void mul(P&z,const P&a,const P&b,const string&label){z.zero();phase=label;xi=xj=-1;if(a.d<0||b.d<0)return;require(a.d+b.d<z.cap,"multiply_degree_guard");
 for(int i=0;i<=a.d;i++){if(fmpz_is_zero(a.c[i].x))continue;for(int j=0;j<=b.d;j++){if(fmpz_is_zero(b.c[j].x))continue;xi=i;xj=j;if((multiplications&4095)==0)guard();fmpz_addmul(z.c[i+j].x,a.c[i].x,b.c[j].x);bitcheck(z.c[i+j].x,RB);multiplications++;}}z.trim();xi=xj=-1;}
static void subtract(P&z,const P&a,const P&b){phase="subtract_products";z.zero();int n=max(a.d,b.d);require(n<z.cap,"subtract_capacity");for(int i=0;i<=n;i++){xi=i;guard();fmpz_sub(z.c[i].x,a.c[i].x,b.c[i].x);bitcheck(z.c[i].x,RB);}z.trim();xi=-1;}
static void exactdivide(P&q,P&r,const P&den){q.zero();phase="exact_division_setup";require(den.d>=0,"zero_previous_pivot");
 while(r.d>=den.d){divshift=r.d-den.d;divindex=-1;phase="exact_division_quotient";guard();require(divshift<=DB&&divshift<q.cap,"quotient_degree_guard");require(fmpz_divisible(r.c[r.d].x,den.c[den.d].x),"nonexact_integer_leading_division");
  fmpz_divexact(coef.x,r.c[r.d].x,den.c[den.d].x);scalarDivisions++;require(fmpz_cmpabs(coef.x,H.x)<=0,"quotient_minor_height_guard");require(fmpz_is_zero(q.c[divshift].x),"duplicate_quotient_degree");fmpz_set(q.c[divshift].x,coef.x);
  for(int j=0;j<=den.d;j++){divindex=j;phase="exact_division_subtract";if((multiplications&4095)==0)guard();fmpz_submul(r.c[j+divshift].x,coef.x,den.c[j].x);bitcheck(r.c[j+divshift].x,RB);multiplications++;}r.trim();
 }require(r.d<0,"nonzero_polynomial_remainder");q.trim();require(q.d<=DB,"minor_degree_guard");for(int i=0;i<=q.d;i++)require(fmpz_cmpabs(q.c[i].x,H.x)<=0,"minor_coefficient_height_guard");divshift=divindex=-1;}
static void checkpoint(const string&why){ofstream o(pathPrefix+"/guard-state-K92AB001.json");o<<"{\"version\":\"K92AB001\",\"reason\":\""<<why<<"\",\"stage\":"<<stage<<",\"row\":"<<ci<<",\"column\":"<<cj<<",\"phase\":\""<<phase<<"\",\"multiply_next_i\":"<<xi<<",\"multiply_next_j\":"<<xj<<",\"division_shift\":"<<divshift<<",\"division_next_j\":"<<divindex<<",\"next_operation_NOT_committed\":true,\"current_quotient_coefficient\":\""<<zz(coef.x)<<"\",\"matrix\":[";
 for(int k=0;k<121;k++){if(k)o<<',';writepoly(o,*M[k]);}o<<"],\"five_buffers\":[";for(int k=0;k<5;k++){if(k)o<<',';writepoly(o,*B[k]);}o<<"],\"row_permutation\":[";for(int k=0;k<11;k++){if(k)o<<',';o<<rowperm[k];}o<<"],\"column_permutation\":[";for(int k=0;k<11;k++){if(k)o<<',';o<<colperm[k];}o<<"],\"sign\":"<<signDet<<",\"integer_products\":"<<multiplications<<",\"exact_integer_divisions\":"<<scalarDivisions<<",\"owned_cells\":"<<owned<<",\"seconds\":"<<elapsed()<<",\"no_rank_or_determinant_result_claim\":true}";}
int main(int argc,char**argv){if(argc==2&&string(argv[1])=="--help"){cout<<"K92AXIS001 inputfile outputdirectory; finite ZZ[a] Bareiss only\n";return 0;}if(argc!=3)return 2;pathPrefix=argv[2];beginTime=chrono::steady_clock::now();
 try{ifstream in(argv[1]);string magic,hs;int n;in>>magic>>n>>DB>>CB>>RB>>OWN>>secLimit;require(magic=="K92AXIS001"&&n==11&&DB==88&&CB==864&&RB==1738&&OWN==27812&&secLimit>0&&secLimit<=45,"input_header");in>>hs;require(fmpz_set_str(H.x,hs.c_str(),10)==0&&fmpz_bits(H.x)==CB,"input_minor_bound");
  M.reserve(121);for(int k=0;k<121;k++)M.emplace_back(make_unique<P>(DB+1));B.reserve(5);for(int k=0;k<5;k++)B.emplace_back(make_unique<P>(2*DB+1));
  for(int k=0;k<121;k++){int i,j,nt;in>>i>>j>>nt;require(i==k/11&&j==k%11&&nt>=0,"input_polynomial_index");int old=-1;for(int h=0;h<nt;h++){int d;string s;in>>d>>s;require(d>old&&d<=8&&d>=0,"input_degree");old=d;require(fmpz_set_str(at(i,j).c[d].x,s.c_str(),10)==0&&!fmpz_is_zero(at(i,j).c[d].x),"input_integer");}at(i,j).trim();}
  string tail;require(!(in>>tail),"input_trailing_data");for(int k=0;k<11;k++)rowperm[k]=colperm[k]=k;fmpz_one(B[4]->c[0].x);B[4]->d=0;int rank=0;
  for(int k=0;k<11;k++){stage=k;ci=cj=-1;phase="pivot_select";guard();int ri=-1,rj=-1;for(int i=k;i<11&&ri<0;i++)for(int j=k;j<11;j++)if(at(i,j).d>=0){ri=i;rj=j;break;}
   if(ri<0){rank=k;break;}if(ri!=k){for(int j=0;j<11;j++)swap(M[ri*11+j],M[k*11+j]);swap(rowperm[ri],rowperm[k]);signDet=-signDet;}if(rj!=k){for(int i=0;i<11;i++)swap(M[i*11+rj],M[i*11+k]);swap(colperm[rj],colperm[k]);signDet=-signDet;}
   require(at(k,k).d>=0,"selected_zero_pivot");rank=k+1;ofstream patch(pathPrefix+"/stage-"+to_string(k)+"-patch.json");patch<<"{\"stage\":"<<k<<",\"row_swap\":["<<k<<','<<ri<<"],\"column_swap\":["<<k<<','<<rj<<"],\"pivot_polynomial\":";writepoly(patch,at(k,k));patch<<",\"previous_pivot_reference\":\""<<(k?"previous_completed_stage_pivot":"constant1")<<"\",\"updates\":[";bool first=true;
   if(k<10){for(int i=k+1;i<11;i++)for(int j=k+1;j<11;j++){ci=i;cj=j;mul(*B[0],at(k,k),at(i,j),"multiply_pivot_oldij");mul(*B[1],at(i,k),at(k,j),"multiply_oldik_oldkj");subtract(*B[2],*B[0],*B[1]);exactdivide(*B[3],*B[2],*B[4]);phase="commit_quotient";guard();copy(at(i,j),*B[3]);updatedCells++;if(!first)patch<<',';first=false;patch<<"{\"row\":"<<i<<",\"column\":"<<j<<",\"quotient\":";writepoly(patch,at(i,j));patch<<",\"remainder_ZERO\":true}";}
    copy(*B[4],at(k,k));for(int i=k+1;i<11;i++)at(i,k).zero();}
   patch<<"],\"column_below_pivot_ZEROed\":true,\"complete\":true}";patch.close();ci=cj=-1;
  }
  require(rank>=6&&rank<=11,"rank_below_known6minor");phase="final";ofstream o(pathPrefix+"/native.json");o<<"{\"all_passed\":true,\"rank_over_QQ_a\":"<<rank<<",\"all_exact_ZZ_polynomial_divisions_remainder_ZERO\":true,\"last_nonzero_pivot_minor\":";writepoly(o,at(rank-1,rank-1));o<<",\"determinant\":";
  if(rank==11){if(signDet<0)for(int h=0;h<=at(10,10).d;h++)fmpz_neg(at(10,10).c[h].x,at(10,10).c[h].x);writepoly(o,at(10,10));}else o<<"[]";
  o<<",\"row_permutation\":[";for(int k=0;k<11;k++){if(k)o<<',';o<<rowperm[k];}o<<"],\"column_permutation\":[";for(int k=0;k<11;k++){if(k)o<<',';o<<colperm[k];}o<<"],\"permutation_determinant_sign\":"<<signDet<<",\"remaining_block_ZERO\":"<<(rank<11?"true":"false")<<",\"updated_cells\":"<<updatedCells<<",\"integer_products\":"<<multiplications<<",\"exact_integer_leading_divisions\":"<<scalarDivisions<<",\"observed_raw_integer_bits\":"<<maxObserved<<",\"owned_polynomial_coefficient_capacity\":"<<owned<<",\"whole_owned_with_input_and_metadata\":"<<owned+1158+15000<<",\"seconds\":"<<elapsed()<<"}";cout<<"rank="<<rank<<" products="<<multiplications<<" cells="<<updatedCells<<" seconds="<<elapsed()<<'\n';return 0;
 }catch(const exception&e){if(M.size()==121&&B.size()==5)checkpoint(e.what());ofstream o(pathPrefix+"/native-failure.json");o<<"{\"all_passed\":false,\"reason\":\""<<e.what()<<"\",\"exact_state_preserved\":true,\"no_rank_or_determinant_claim\":true}";cerr<<e.what()<<'\n';return 1;}
}
