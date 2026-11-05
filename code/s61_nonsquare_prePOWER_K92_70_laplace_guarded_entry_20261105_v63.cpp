// ONE changed association, reuses completed source values; never calls sourceEvaluate.
#define main frozen_old_entry_NOT_called
#include "s61_nonsquare_prePOWER_K92_det11_guarded_entry_20261104_v63.cpp"
#undef main

constexpr uint64_t TOTAL_RAW_CAP=32ULL*1024*1024;
int group70=-1,term70=-1,groupsDone70=0,top70=-1,bottom70=-1,weighted70=-1,acc70=-1;
string stage70="selective_restore",oldStateHash,layoutHash;
vector<vector<Poly*>> matrix70;
vector<Poly*> weight70;
uint64_t skipped70=0,loaded70=0,written70=0,writeLimit70=0;

template<class T>T get70(istream&f){using U=make_unsigned_t<T>;U v=0;for(size_t i=0;i<sizeof(T);i++){int c=f.get();assert(c!=EOF);v|=U(uint8_t(c))<<(8*i);}return static_cast<T>(v);}
void readZ70(istream&f,Big&z,int limit,bool keep){uint16_t n=get70<uint16_t>(f);int neg=f.get();assert((neg==0||neg==1)&&n<=(limit+7)/8);if(!keep){f.ignore(n);assert(f.good());return;}fmpz_zero(z.z);for(int i=0;i<n;i++){int c=f.get();assert(c!=EOF&&(i||c));fmpz_mul_2exp(z.z,z.z,8);fmpz_add_ui(z.z,z.z,uint8_t(c));}if(neg)fmpz_neg(z.z,z.z);assert(bits(z.z)<=limit);}
void raw70(ostream&o,const void*p,size_t n){assert(written70+n<=writeLimit70);o.write(static_cast<const char*>(p),n);assert(o.good());written70+=n;}
template<class T>void put70(ostream&o,T x){using U=make_unsigned_t<T>;U v=static_cast<U>(x);unsigned char b[sizeof(T)];for(size_t i=0;i<sizeof(T);i++)b[i]=(v>>(8*i))&255;raw70(o,b,sizeof(T));}
void putZ70(ostream&o,const fmpz_t z,int limit){assert(bits(z)<=limit);char*p=fmpz_get_str(nullptr,16,z);string h=p;flint_free(p);bool neg=!h.empty()&&h[0]=='-';if(neg)h.erase(h.begin());if(h.size()%2)h='0'+h;uint16_t n=fmpz_is_zero(z)?0:h.size()/2;put70<uint16_t>(o,n);put70<uint8_t>(o,neg);for(int i=0;i<n;i++)put70<uint8_t>(o,stoul(h.substr(2*i,2),nullptr,16));}
string ids70(){ostringstream h;h<<"[";for(size_t i=0;i<matrix70.size();i++){if(i)h<<",";h<<"[";for(size_t j=0;j<matrix70[i].size();j++){if(j)h<<",";h<<matrix70[i][j]->id;}h<<"]";}h<<"]";return h.str();}
void checkpoint70(bool complete,const string&why){
 assert(alive.size()<=512&&cells<=CELL_CAP);uint64_t existing=0;for(auto&e:filesystem::directory_iterator(folder))if(e.is_regular_file())existing+=e.file_size();assert(existing<6*1024*1024);writeLimit70=TOTAL_RAW_CAP-existing;written70=0;
 ostringstream h;h<<"{\"version\":\"K9270P01\",\"old_state_sha256\":\""<<oldStateHash<<"\",\"layout_sha256\":\""<<layoutHash<<"\",\"complete\":"<<(complete?"true":"false")<<",\"stage\":\""<<stage70<<"\",\"reason\":\""<<jsonEscape(why)<<"\",\"group\":"<<group70<<",\"term\":"<<term70<<",\"groups_completed\":"<<groupsDone70<<",\"terms_completed\":"<<completedTerms<<",\"top_id\":"<<top70<<",\"bottom_id\":"<<bottom70<<",\"weighted_id\":"<<weighted70<<",\"acc_id\":"<<acc70<<",\"matrix_ids\":"<<ids70()<<",\"weight_ids\":[";
 for(size_t i=0;i<weight70.size();i++){if(i)h<<",";h<<weight70[i]->id;}h<<"],\"cursor_phase\":\""<<cur.phase<<"\",\"row\":"<<cur.row<<",\"mask\":"<<cur.mask<<",\"column\":"<<cur.column<<",\"updates\":"<<updates<<",\"source_DAG_evaluations\":0,\"reused_input_terms\":"<<loaded70<<",\"skipped_old_terms\":"<<skipped70<<",\"known_output_a_shift\":"<<outputAShift<<",\"peak_cells\":"<<peakcells<<",\"peak_polys\":"<<peakpolys<<"}";
 string header=h.str();assert(header.size()<65536);ofstream o(folder+"/state.bin",ios::binary);raw70(o,"K9270P01",8);put70<uint32_t>(o,header.size());raw70(o,header.data(),header.size());put70<uint32_t>(o,alive.size());
 for(auto&[id,p]:alive){put70<uint32_t>(o,id);put70<uint8_t>(o,p->role);put70<int32_t>(o,p->source);put70<int32_t>(o,p->term);put70<int32_t>(o,p->row);put70<int32_t>(o,p->mask);put70<uint8_t>(o,p->complete);putZ70(o,fmpq_numref(p->c.q),CONTENT_BITS);putZ70(o,fmpq_denref(p->c.q),CONTENT_BITS);put70<uint32_t>(o,p->a.size());for(auto&[m,c]:p->a){for(int e:m)put70<int32_t>(o,e);putZ70(o,c->z,COEFF_BITS);}}
 put70<uint8_t>(o,cur.pending);if(cur.pending){put70<int32_t>(o,cur.dst);put70<int32_t>(o,cur.lhs);put70<int32_t>(o,cur.rhs);put70<uint64_t>(o,cur.left);put70<uint64_t>(o,cur.right);for(auto x:cur.key)put70<int64_t>(o,x);putZ70(o,cur.scale.z,PENDING_BITS);put70<uint8_t>(o,cur.valueKnown);putZ70(o,cur.value.z,PENDING_BITS);}o.flush();assert(o.good());assert(existing+written70<=TOTAL_RAW_CAP);
}
void status70(bool complete,const string&why){ofstream o(folder+"/native.json");o<<"{\"all_passed\":"<<(complete?"true":"false")<<",\"full_scalar_complete\":"<<(complete?"true":"false")<<",\"output_nonzero\":"<<(complete?(alive.at(acc70)->a.empty()?"false":"true"):"null")<<",\"output_terms\":"<<(acc70>=0&&alive.count(acc70)?alive.at(acc70)->a.size():0)<<",\"groups_completed\":"<<groupsDone70<<",\"H3_terms_completed\":"<<completedTerms<<",\"source_DAG_evaluations\":0,\"reused_source_terms\":"<<loaded70<<",\"old_terms_skipped_without_allocation\":"<<skipped70<<",\"coefficient_updates\":"<<updates<<",\"peak_owned_cells\":"<<peakcells<<",\"peak_owned_plus_work\":"<<peakcells+15000<<",\"peak_active_polys\":"<<peakpolys<<",\"reason\":\""<<jsonEscape(why)<<"\",\"all_active_state_saved\":true,\"state_bytes\":"<<written70<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";}

void restoreInputs70(const string&path,const vector<vector<int>>&roots){
 set<int>wanted;for(auto&r:roots)for(int i:r)wanted.insert(i);for(int i=593;i<=598;i++)wanted.insert(i);
 ifstream f(path,ios::binary);char magic[8];f.read(magic,8);assert(string(magic,8)=="K92SP001");uint32_t n=get70<uint32_t>(f);assert(n<65536);string h(n,' ');f.read(h.data(),n);assert(h.find("\"source_completed\":572")!=string::npos);uint32_t np=get70<uint32_t>(f);assert(np<=1024);map<int,Poly*>src,lit;
 for(uint32_t i=0;i<np;i++){get70<uint32_t>(f);int role=get70<uint8_t>(f),s=get70<int32_t>(f),t=get70<int32_t>(f);get70<int32_t>(f);get70<int32_t>(f);bool done=get70<uint8_t>(f);bool keep=wanted.count(s)||(role==2&&t>=0&&t<4);Big num,den;readZ70(f,num,CONTENT_BITS,keep);readZ70(f,den,CONTENT_BITS,keep);Poly*p=nullptr;
  if(keep){assert(done&&fmpz_sgn(num.z)>0&&fmpz_sgn(den.z)>0);Rat c;fmpq_set_fmpz_frac(c.q,num.z,den.z);checkRat(c);p=makePoly(c,role,s,t);p->input=true;}
  uint32_t nt=get70<uint32_t>(f);for(uint32_t j=0;j<nt;j++){Mon m;for(int&x:m)x=get70<int32_t>(f);Big v;readZ70(f,v,COEFF_BITS,keep);if(keep){assert(!fmpz_is_zero(v.z));insertInput(p,m,v.z);}}
  if(keep){p->complete=true;loaded70+=nt;if(s>=0)src[s]=p;else lit[t]=p;}else skipped70+=nt;
 }
 bool pending=get70<uint8_t>(f);if(pending){for(int i=0;i<3;i++)get70<int32_t>(f);for(int i=0;i<2;i++)get70<uint64_t>(f);for(int i=0;i<3;i++)get70<int64_t>(f);Big z;readZ70(f,z,PENDING_BITS,false);get70<uint8_t>(f);readZ70(f,z,PENDING_BITS,false);}assert(f.peek()==EOF&&loaded70==25266&&skipped70==109734);
 matrix70.assign(11,vector<Poly*>(8));for(int i=0;i<11;i++)for(int j=0;j<8;j++){assert(src.count(roots[i][j]));matrix70[i][j]=src.at(roots[i][j]);}
 assert(lit.size()==4);weight70={lit.at(0),lit.at(1),src.at(593),lit.at(2),lit.at(3),src.at(594),src.at(597),src.at(596),src.at(595),src.at(598)};assert(alive.size()==98&&cells==25266);cur.pending=false;
}

void fixture70(){stage70="fixture";Rat ca=parseRat("1/6"),cb=parseRat("1/5");Poly*a=makePoly(ca,5),*b=makePoly(cb,5);Big z;fmpz_set_si(z.z,3);insertInput(a,{1,0,0},z.z);fmpz_set_si(z.z,-2);insertInput(a,{0,1,0},z.z);fmpz_set_si(z.z,5);insertInput(b,{0,0,1},z.z);fmpz_one(z.z);insertInput(b,{0,0,0},z.z);Rat c=qmul(ca,cb);Poly*p=makePoly(c,5);addMul(p,a,b,z.z);assert(p->a.size()==4&&fmpz_equal_si(p->a.at({1,0,1})->z,15));a->complete=b->complete=p->complete=true;Poly*t=makePoly(Rat(),5);long saved=cells;bool stop=false;cells=CELL_CAP;try{insertInput(t,{0,0,0},z.z);}catch(const Stop&){stop=true;}cells=saved;assert(stop&&t->a.empty());fmpz_set_si(z.z,7);insertInput(t,{0,0,0},z.z);cur.dst=t->id;cur.key={0,0,0};cur.pending=true;cur.valueKnown=true;fmpz_one(cur.value.z);fmpz_mul_2exp(cur.value.z,cur.value.z,1024);stop=false;try{commit(t);}catch(const Stop&){stop=true;}assert(stop&&fmpz_equal_si(t->a.at({0,0,0})->z,7));t->complete=true;cur.pending=false;checkpoint70(true,"fixture");}

int main(int argc,char**argv){began=chrono::steady_clock::now();assert(argc==3||argc==6);folder=argv[2];int code=0;
 try{
  if(string(argv[1])=="--fixture"){fixture70();alive.clear();flint_cleanup();return 0;}
  assert(argc==6);oldStateHash=argv[3];layoutHash=argv[5];ifstream layout(argv[4]);int nr,nc;layout>>nr>>nc;assert(nr==11&&nc==8);vector<vector<int>>roots(11,vector<int>(8));for(auto&r:roots)for(int&i:r)layout>>i;layout>>expectedGlobalC>>outputAShift;assert(outputAShift==5);string extra;assert(!(layout>>extra));restoreInputs70(argv[1],roots);
  Rat gc=parseRat(expectedGlobalC);Poly*acc=makePoly(gc,4);acc70=acc->id;vector<array<int,3>>selected={{0,2,4},{0,2,5},{0,4,5},{1,2,4},{1,2,5},{1,4,5},{2,3,4},{2,3,5},{2,4,5},{3,4,5}};vector<int>sign={1,-1,-1,-1,1,1,-1,-1,-1,1};
  for(int a=0;a<5;a++)for(int b=a+1;b<6;b++)for(int c=b+1;c<7;c++)for(int d=c+1;d<8;d++){
   group70++;term70=-1;vector<int>cols={a,b,c,d},other;for(int j=0;j<8;j++)if(find(cols.begin(),cols.end(),j)==cols.end())other.push_back(j);
   stage70="top4";vector<vector<Poly*>>T(4,vector<Poly*>(4));vector<vector<bool>>ts(4,vector<bool>(4));for(int i=0;i<4;i++)for(int j=0;j<4;j++){T[i][j]=matrix70[i][cols[j]];ts[i][j]=!T[i][j]->a.empty();}Poly*top=determinant(T,ts,4,1000+group70);top70=top->id;
   for(int ti=0;ti<10;ti++){
    term70=ti;stage70="bottom4";vector<int>rr={4};for(int i=0;i<6;i++)if(find(selected[ti].begin(),selected[ti].end(),i)==selected[ti].end())rr.push_back(5+i);vector<vector<Poly*>>B(4,vector<Poly*>(4));vector<vector<bool>>bs(4,vector<bool>(4));for(int i=0;i<4;i++)for(int j=0;j<4;j++){B[i][j]=matrix70[rr[i]][other[j]];bs[i][j]=!B[i][j]->a.empty();}Poly*bottom=determinant(B,bs,4,group70*10+ti);bottom70=bottom->id;
    stage70="weighted_bottom";Rat wc=qmul(weight70[ti]->c,bottom->c);Poly*wb=makePoly(wc,6,-1,ti);weighted70=wb->id;Big one;fmpz_one(one.z);addMul(wb,weight70[ti],bottom,one.z);wb->complete=true;destroy(bottom);bottom70=-1;
    stage70="global_accumulate";Big scale;integerRatio(scale,top->c,wb->c,acc->c);if(sign[ti]*((a+b+c+d)%2?-1:1)<0)fmpz_neg(scale.z,scale.z);addMul(acc,top,wb,scale.z,-outputAShift);destroy(wb);weighted70=-1;completedTerms++;clockcheck();
   }
   destroy(top);top70=-1;groupsDone70++;stage70="group_complete";clockcheck();
  }
  assert(group70==69&&groupsDone70==70&&completedTerms==700);acc->complete=true;stage70="complete";cur.pending=false;for(auto&[m,z]:acc->a)assert(m[0]>=0&&m[1]>=0&&m[2]>=0);checkpoint70(true,"");status70(true,"");
 }catch(const Stop&e){code=2;checkpoint70(false,e.what());status70(false,e.what());}
 catch(const exception&e){code=3;try{checkpoint70(false,string("error:")+e.what());status70(false,string("error:")+e.what());}catch(...){ofstream f(folder+"/fatal.txt");f<<e.what()<<'\n';}}
 alive.clear();assert(cells==0);flint_cleanup();return code;
}
