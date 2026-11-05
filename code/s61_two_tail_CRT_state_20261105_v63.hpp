#pragma once
#include <flint/fmpz.h>
#include <gmp.h>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
namespace tt {
struct Z{fmpz_t v;Z(){fmpz_init(v);}Z(const Z&o){fmpz_init(v);fmpz_set(v,o.v);}Z(Z&&o)noexcept{fmpz_init(v);fmpz_swap(v,o.v);}~Z(){fmpz_clear(v);}};
struct Rec{uint32_t key;Z c;};using Poly=std::map<uint64_t,Z>;
struct Stop:std::runtime_error{using std::runtime_error::runtime_error;};
inline constexpr uint32_t P[40]={32003,32009,32027,32029,32051,32057,32059,32063,32069,32077,32083,32089,32099,32117,32119,32141,32143,32159,32173,32183,32189,32191,32203,32213,32233,32237,32251,32257,32261,32297,32299,32303,32309,32321,32323,32327,32341,32353,32359,32363};
inline constexpr uint64_t WRAD=uint64_t(1)<<26,MA=uint64_t(1)<<20,ME=uint64_t(1)<<10;
inline constexpr size_t OWN=60000,META=16000;
enum F{MODE,PI,EUST,EUPASS,SHIFT,JNEXT,IL,FACTOR,NORMI,CRTI,CRSUB,CRRES,CRU,CRINV,ACCEPTED,NORM_STAGE,IDX,WHICH,INT_STAGE,ACTIVE_READY,ACTIVE_KEY,HKEY,QKEY,HINDEX,MODOPS,PRIMEOPS,BIGOPS,PEAK,SIGN,QB0,RB0,QB2,RB2,REXP,VERIFIED0,VERIFIED2,ELAPSED_MS,HAS_CRT,HAS_FACTOR,INDEG0,INDEG2,PREV0,PREV2};
enum ZZ{GAMMA,MODULUS,NEWMOD,CONTENT,LC,NORM,PENDING,ACTIVE,TMP,L10,L12,LC0,LC2,D,ZCOUNT};
enum Modes{LOAD=1,SEED=2,FIELD=3,MONIC=4,CRT=5,BALANCE=6,CONTENT_GCD=7,PRIMITIVE=8,CANDIDATE=9,HEIGHT=10,INT_LOAD=11,INT_DIV=12,VERIFY_LOAD=13,VERIFY_MUL=14,DONE=15,FIXTURE=99};
struct State{
 std::array<uint64_t,64> u{};std::array<Z,ZCOUNT> z;
 std::array<std::vector<Rec>,2> in;std::vector<uint32_t> f,g;std::vector<Z> crt;
 Poly H,R,Q;std::array<std::array<uint64_t,8>,40> ledger{};
 std::array<std::string,3> paths;std::chrono::steady_clock::time_point started;
 size_t owned()const{return META+2*(in[0].capacity()+in[1].capacity()+H.size()+R.size()+Q.size())+f.capacity()+g.capacity()+crt.capacity();}
 void own(size_t extra=0){size_t x=owned()+extra;if(x>OWN)throw Stop("60000_owned_before_commit");u[PEAK]=std::max<uint64_t>(u[PEAK],x);}
 void clock(){u[ELAPSED_MS]=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-started).count();if(u[ELAPSED_MS]>45000)throw Stop("45s_soft_prefix_stop");}
 void bit(const Z&x,unsigned b){if(fmpz_bits(x.v)>b)throw Stop("integer_bit_before_commit");}
 void modreserve(uint64_t n){if(u[MODOPS]+n>6000000000ULL||u[PRIMEOPS]+n>150000000ULL)throw Stop("modular_operation_before_commit");}
 uint32_t mul(uint32_t a,uint32_t b,uint32_t p){modreserve(1);++u[MODOPS];++u[PRIMEOPS];return uint32_t((uint64_t(a)*b)%p);}
 uint32_t inverse(uint32_t a,uint32_t p){if(!a)throw Stop("zero_modular_inverse");modreserve(64);uint32_t r=1;for(uint32_t n=p-2;n;n>>=1,a=mul(a,a,p))if(n&1)r=mul(r,a,p);return r;}
 void bigreserve(uint64_t n=1){if(u[BIGOPS]+n>2100000ULL)throw Stop("bigint_operation_before_commit");}
};
inline std::string zs(const Z&z){char*p=fmpz_get_str(nullptr,10,z.v);std::string s(p);flint_free(p);return s;}
inline void xyz(uint64_t legacy,unsigned&a,unsigned&e,unsigned&d){a=legacy/MA;e=(legacy%MA)/ME;d=legacy%ME;}
inline uint64_t weighted(uint64_t legacy){unsigned a,e,d;xyz(legacy,a,e,d);return uint64_t(2*a+e)*WRAD+legacy;}
inline uint32_t kron(uint64_t legacy){unsigned a,e,d;xyz(legacy,a,e,d);return 646*a+19*e+d;}
inline void release(std::vector<uint32_t>&x){std::vector<uint32_t>().swap(x);}
inline void trim(std::vector<uint32_t>&x){while(!x.empty()&&!x.back())x.pop_back();}
struct Writer{
 std::ofstream o;uint64_t n=0,limit;
 Writer(const std::string&p,uint64_t lim):o(p,std::ios::binary),limit(lim){if(!o)throw Stop("artifact_open_IO_failure");}
 void bytes(const void*p,size_t m){if(n+m>limit)throw Stop("binary_output_before_write_bound");o.write(static_cast<const char*>(p),m);n+=m;if(!o)throw Stop("artifact_device_IO_failure");}
 void b(uint8_t x){bytes(&x,1);}void w(uint32_t x){uint8_t z[4];for(int i=0;i<4;++i)z[i]=x>>(8*i);bytes(z,4);}void q(uint64_t x){uint8_t z[8];for(int i=0;i<8;++i)z[i]=x>>(8*i);bytes(z,8);}
 void text(const std::string&s){if(s.size()>4096)throw Stop("metadata_string_bound");w(s.size());bytes(s.data(),s.size());}
 void z(const Z&x,unsigned cap){if(fmpz_bits(x.v)>cap)throw Stop("binary_scalar_type_bound");unsigned char buf[2048];size_t m=0;mpz_t v;mpz_init(v);fmpz_get_mpz(v,x.v);bool neg=mpz_sgn(v)<0;mpz_abs(v,v);mpz_export(buf,&m,-1,1,0,0,v);mpz_clear(v);if(m>sizeof(buf))throw Stop("binary_scalar_buffer_bound");b(neg);w(m);bytes(buf,m);}
 void map(const Poly&p,unsigned cap){w(p.size());for(auto&[key,c]:p){q(key);z(c,cap);}}
 void close(){o.flush();if(!o)throw Stop("artifact_flush_IO_failure");o.close();}
};
inline void savePrefix(State&s,const std::string&path,const std::string&why){
 s.u[ELAPSED_MS]=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-s.started).count(); // Never rethrow the time guard while saving it.
 Writer o(path,16*1024*1024);o.bytes("TTCRT001",8);o.text(why);for(auto x:s.u)o.q(x);for(auto&x:s.z)o.z(x,8193);
 for(auto&row:s.ledger)for(auto x:row)o.q(x);for(auto&p:s.paths)o.text(p);
 for(auto&v:s.in){o.w(v.capacity());o.w(v.size());for(auto&r:v){o.w(r.key);o.z(r.c,301);}}
 for(auto*v:{&s.f,&s.g}){o.w(v->capacity());o.w(v->size());for(auto x:*v)o.w(x);}
 o.w(s.crt.capacity());o.w(s.crt.size());for(auto&x:s.crt)o.z(x,600);
 o.map(s.H,600);o.map(s.R,8192);o.map(s.Q,8192);o.close();
}
inline void saveMap(const Poly&p,const std::string&path,const char*magic,uint32_t which,uint64_t limit,unsigned cap){Writer o(path,limit);o.bytes(magic,8);o.w(which);o.map(p,cap);o.close();}
inline void saveCRT(State&s,const std::string&path){Writer o(path,1024*1024);o.bytes("TTCRTF01",8);o.z(s.z[MODULUS],600);o.z(s.z[GAMMA],600);o.z(s.z[CONTENT],600);o.w(s.crt.size());for(auto&x:s.crt)o.z(x,600);o.close();}
}
