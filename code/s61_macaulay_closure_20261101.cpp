// Exact bounded degree-three closure over a prime field.  Newly exposed
// quadratic/linear rows are multiplied by every variable and fed back.
// The binary state is scratch data: it preserves expensive elimination for
// a later bounded continuation, and is not itself a research certificate.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
using Row=std::vector<std::pair<int,int>>;
template<class T> void put(std::ostream&o,const T&v){o.write(reinterpret_cast<const char*>(&v),sizeof(T));}
template<class T> void get(std::istream&i,T&v){i.read(reinterpret_cast<char*>(&v),sizeof(T));}
void putrow(std::ostream&o,const Row&r){uint32_t n=r.size();put(o,n);for(auto e:r){uint16_t k=e.first,v=e.second;put(o,k);put(o,v);}}
void getrow(std::istream&i,Row&r){uint32_t n;get(i,n);r.resize(n);for(auto&e:r){uint16_t k,v;get(i,k);get(i,v);e={k,v};}}
void putrows(std::ostream&o,const std::vector<Row>&r){uint32_t n=r.size();put(o,n);for(auto&a:r)putrow(o,a);}
void getrows(std::istream&i,std::vector<Row>&r){uint32_t n;get(i,n);r.resize(n);for(auto&a:r)getrow(i,a);}
int power(int a,int n,int p){int r=1;for(;n;n>>=1,a=int(int64_t(a)*a%p))if(n&1)r=int(int64_t(r)*a%p);return r;}
int main(int argc,char**argv){
 if(argc!=6)return 2;
 std::ifstream in(argv[1]);int prime,ncols,nvars,nseeds;in>>prime>>ncols>>nvars>>nseeds;
 if(prime>=65536||ncols>=65536||nvars!=45)return 3;
 std::vector<int>degree(ncols),lookup(46*46*46,-1);std::vector<std::vector<int>>mon(ncols);
 for(int k=0;k<ncols;k++){in>>degree[k];mon[k].resize(degree[k]);int code=0,mul=1;for(auto&v:mon[k]){in>>v;code+=(v+1)*mul;mul*=46;}lookup[code]=k;}
 std::vector<Row> seeds(nseeds);for(auto&r:seeds){int n;in>>n;r.resize(n);for(auto&e:r)in>>e.first>>e.second;}
 if(!in)return 4;
 std::vector<std::vector<int>>shift(nvars,std::vector<int>(ncols,-1));
 for(int j=0;j<nvars;j++)for(int k=0;k<ncols;k++)if(degree[k]<=2){auto m=mon[k];m.push_back(j);std::sort(m.begin(),m.end());int code=0,mul=1;for(int v:m){code+=(v+1)*mul;mul*=46;}shift[j][k]=lookup[code];}
 std::vector<Row> basis(ncols),pending=seeds,low_new=seeds,discovered;
 uint64_t cursor=0,processed=0,entries=0;int rank=0,round=0,counts[4]={0};
 bool seeded=false;std::string resume=argv[5];
 if(resume!="-"){
  std::ifstream state(resume,std::ios::binary);uint64_t magic;get(state,magic);if(magic!=0x533631434c4f5331ULL)return 5;
  int saved_prime,saved_cols;get(state,saved_prime);get(state,saved_cols);if(saved_prime!=prime||saved_cols!=ncols)return 6;
  get(state,cursor);get(state,processed);get(state,entries);get(state,rank);get(state,round);get(state,seeded);for(auto&v:counts)get(state,v);
  getrows(state,basis);getrows(state,pending);getrows(state,low_new);getrows(state,discovered);if(!state)return 7;
 }
 double seconds=std::stod(argv[2]);uint64_t entry_limit=std::stoull(argv[3]);auto start=std::chrono::steady_clock::now();std::string status="complete";
 while(true){
  if(cursor==pending.size()){
   if(seeded&&low_new.empty())break;
   pending.clear();cursor=0;
   for(int j=0;j<nvars;j++)for(const Row&r:low_new){Row product;product.reserve(r.size());for(auto e:r)product.emplace_back(shift[j][e.first],e.second);std::sort(product.begin(),product.end());pending.push_back(std::move(product));}
   low_new.clear();seeded=true;round++;
   std::cerr<<"round="<<round<<" pending="<<pending.size()<<" rank="<<rank<<" low="<<counts[1]<<","<<counts[2]<<std::endl;
  }
  Row row=pending[cursor];
  while(!row.empty()){
   if(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()>seconds){status="time_limit";goto finish;}
   if(entries>entry_limit){status="entry_limit";goto finish;}
   int pivot=row.front().first;
   if(basis[pivot].empty()){
    int inv=power(row.front().second,prime-2,prime);for(auto&e:row)e.second=int(int64_t(e.second)*inv%prime);
    if(seeded&&degree[pivot]<=2){low_new.push_back(row);discovered.push_back(row);std::cerr<<"new_low degree="<<degree[pivot]<<" terms="<<row.size()<<" processed="<<processed<<std::endl;}
    entries+=row.size();basis[pivot]=std::move(row);rank++;counts[degree[pivot]]++;break;
   }
   const Row&old=basis[pivot];int scale=row.front().second;Row next;next.reserve(row.size()+old.size());size_t i=0,j=0;
   while(i<row.size()||j<old.size()){
    int key,value;
    if(j==old.size()||(i<row.size()&&row[i].first<old[j].first)){key=row[i].first;value=row[i++].second;}
    else if(i==row.size()||old[j].first<row[i].first){key=old[j].first;value=int((prime-int64_t(scale)*old[j++].second%prime)%prime);}
    else{key=row[i].first;value=int((row[i++].second-int64_t(scale)*old[j++].second%prime+prime)%prime);}
    if(value)next.emplace_back(key,value);
   }row.swap(next);
  }
  cursor++;processed++;
 }
 finish:
 // A partially reduced current row is discarded and retried from its original
 // pending row upon continuation.  The stored basis only contains completed rows.
 std::ofstream state(argv[4],std::ios::binary);uint64_t magic=0x533631434c4f5331ULL;put(state,magic);put(state,prime);put(state,ncols);
 put(state,cursor);put(state,processed);put(state,entries);put(state,rank);put(state,round);put(state,seeded);for(auto v:counts)put(state,v);
 putrows(state,basis);putrows(state,pending);putrows(state,low_new);putrows(state,discovered);state.close();
 std::cout<<"{\"status\":\""<<status<<"\",\"prime\":"<<prime<<",\"round\":"<<round<<",\"processed\":"<<processed<<",\"pending_cursor\":"<<cursor<<",\"pending_size\":"<<pending.size()<<",\"rank\":"<<rank<<",\"pivot_counts\":["<<counts[0]<<","<<counts[1]<<","<<counts[2]<<","<<counts[3]<<"],\"entries\":"<<entries<<",\"elapsed_seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<",\"new_lower_degree_rows\":[";
 for(size_t i=0;i<discovered.size();i++){if(i)std::cout<<',';std::cout<<'[';for(size_t j=0;j<discovered[i].size();j++){if(j)std::cout<<',';std::cout<<'['<<discovered[i][j].first<<','<<discovered[i][j].second<<']';}std::cout<<']';}std::cout<<"]}"<<std::endl;
}
