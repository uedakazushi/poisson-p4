// Fixed weight/K only. Integer truncated subset determinant, no parameter point.
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
constexpr int K=32,N=11;
struct Z{fmpz_t z;Z(){fmpz_init(z);}~Z(){fmpz_clear(z);}};
long owned=0,peak=0,inputSlots=0;uint64_t updates=0;int rowCur=-1,maskCur=-1,colCur=-1,degreeCur=-1,gapCur=-1,bitsBound=0;bool complete=false;string folder;
chrono::steady_clock::time_point begun;
struct Stop:runtime_error{using runtime_error::runtime_error;};
void check(bool p,const string&s){if(!p)throw runtime_error(s);}
void clockcheck(){if(chrono::duration<double>(chrono::steady_clock::now()-begun).count()>45)throw Stop("45s_guard");}
struct U{array<fmpz_t,K+1>c;U(){if(owned+K+1+inputSlots>135000)throw Stop("owned_coefficient_guard");for(auto&z:c)fmpz_init(z);owned+=K+1;peak=max(peak,owned+inputSlots);}~U(){for(auto&z:c)fmpz_clear(z);owned-=K+1;}};
struct Term{int gap;unique_ptr<Z>value;};
vector<vector<vector<Term>>>M(N,vector<vector<Term>>(N));
map<int,unique_ptr<U>>previous,nextLayer;
string str(const fmpz_t z){char*s=fmpz_get_str(nullptr,10,z);string r=s;flint_free(s);return r;}
void savePrefix(){ofstream f(folder+"/partial-state.json");f<<"{\"version\":\"K92T0321\",\"complete\":false,\"row\":"<<rowCur<<",\"mask\":"<<maskCur<<",\"column\":"<<colCur<<",\"degree\":"<<degreeCur<<",\"gap\":"<<gapCur<<",\"pending_update_uncommitted\":true,\"previous\":[";bool comma=false;
 auto writeLayer=[&](const map<int,unique_ptr<U>>&m){for(auto&[mask,p]:m){if(comma)f<<",";comma=true;f<<"["<<mask<<",[";for(int k=0;k<=K;k++){if(k)f<<",";f<<"\""<<str(p->c[k])<<"\"";}f<<"]]";check(f.tellp()<16*1024*1024,"prefix_RAW_guard");}};
 writeLayer(previous);f<<"],\"next\":[";comma=false;writeLayer(nextLayer);f<<"]}\n";check(f.good(),"prefix_write");}
int main(int argc,char**argv){check(argc==3,"args");begun=chrono::steady_clock::now();folder=argv[2];int code=0;Z Dprod;fmpz_one(Dprod.z);
 try{ifstream f(argv[1]);int n,k;f>>n>>k>>bitsBound;check(n==N&&k==K&&bitsBound>0&&bitsBound<1024,"fixed_input");for(int i=0;i<N;i++){string s;f>>s;Z d;check(!fmpz_set_str(d.z,s.c_str(),10)&&fmpz_sgn(d.z)>0,"row_LCM");fmpz_mul(Dprod.z,Dprod.z,d.z);}
  for(int i=0;i<N;i++)for(int j=0;j<N;j++){int nt;f>>nt;check(nt>=0&&nt<=K+1,"entry_count");int prev=-1;for(int t=0;t<nt;t++){int gap;string s;f>>gap>>s;check(gap>prev&&gap<=K,"gap");prev=gap;auto z=make_unique<Z>();check(!fmpz_set_str(z->z,s.c_str(),10)&&!fmpz_is_zero(z->z)&&fmpz_bits(z->z)<=bitsBound,"entry_integer");M[i][j].push_back({gap,move(z)});inputSlots++;}}
  string extra;check(!(f>>extra),"input_tail");check(inputSlots==78,"selected78");previous[0]=make_unique<U>();fmpz_one(previous[0]->c[0]);Z tmp;
  for(int row=0;row<N;row++){rowCur=row;for(auto it=previous.begin();it!=previous.end();){int mask=it->first;maskCur=mask;U*p=it->second.get();for(int j=0;j<N;j++)if(!(mask>>j&1)&&!M[row][j].empty()){colCur=j;int nm=mask|1<<j;if(!nextLayer.count(nm))nextLayer[nm]=make_unique<U>();U*q=nextLayer[nm].get();bool neg=__builtin_popcount(unsigned(mask)>>(j+1))%2;
    for(int k=0;k<=K;k++)if(!fmpz_is_zero(p->c[k]))for(auto&t:M[row][j])if(k+t.gap<=K){degreeCur=k;gapCur=t.gap;if(!(updates&1023))clockcheck();fmpz_mul(tmp.z,p->c[k],t.value->z);if(neg)fmpz_neg(tmp.z,tmp.z);fmpz_add(tmp.z,tmp.z,q->c[k+t.gap]);if(fmpz_bits(tmp.z)>bitsBound)throw Stop("BEFORE_COMMIT_integer_height");fmpz_set(q->c[k+t.gap],tmp.z);updates++;}
   }it=previous.erase(it);}previous=move(nextLayer);nextLayer.clear();clockcheck();}
  check(previous.count((1<<N)-1)==1,"fullmask");complete=true;ofstream out(folder+"/native.json");out<<"{\"all_passed\":true,\"K\":32,\"row_denominator_product\":\""<<str(Dprod.z)<<"\",\"cleared_integer_coefficients\":[";U*p=previous.at((1<<N)-1).get();for(int k=0;k<=K;k++){if(k)out<<",";out<<"\""<<str(p->c[k])<<"\"";}out<<"],\"coefficient_updates\":"<<updates<<",\"peak_initialized_integer_slots\":"<<peak<<",\"peak_plus15000scratch\":"<<peak+15000<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-begun).count()<<"}\n";check(out.good(),"output");
 }catch(const exception&e){code=2;savePrefix();ofstream o(folder+"/native-failure.json");o<<"{\"all_passed\":false,\"reason\":\""<<e.what()<<"\",\"coefficient_updates\":"<<updates<<",\"peak_initialized_integer_slots\":"<<peak<<"}\n";}
 previous.clear();nextLayer.clear();check(owned==0,"released");flint_cleanup();return code;
}
