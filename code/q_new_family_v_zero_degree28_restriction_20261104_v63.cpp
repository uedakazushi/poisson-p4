// One exactQQ source28 restriction, reusing the reviewed ambient16 gradients.
#define main unused_frozen_v0_26_main
#include "q_new_family_v_zero_degree26_restriction_20261104_v63.cpp"
#undef main
P readqpoly(istream&in){P p;int n;in>>n;assert(n>=0);fmpq_t c;fmpq_init(c);ulong e[4];for(int i=0;i<n;i++){string q;in>>q;for(auto&z:e)in>>z;assert(fmpq_set_str(c,q.c_str(),10)==0);fmpq_canonicalise(c);fmpq_mpoly_set_coeff_fmpq_ui(p.v,c,e,ctx);}fmpq_clear(c);return p;}
P one28(){P p;fmpq_mpoly_one(p.v,ctx);return p;}
P power28(const P&a,int n){P p=one28(),b=a;while(n){if(n&1)p=mul(p,b);n>>=1;if(n)b=mul(b,b);}return p;}
int main(int argc,char**argv){
 assert(argc==3);auto started=chrono::steady_clock::now();fmpq_mpoly_ctx_init(ctx,4,ORD_DEGLEX);int status=0;{
 try{
  string dir=argv[2];ifstream in(argv[1]);int dim,width;in>>dim>>width;assert(dim==16&&width==136);vector<P>beta;for(int i=0;i<16;i++)beta.push_back(readqpoly(in));
  vector<pair<int,int>>mon(136);for(auto&x:mon)in>>x.first>>x.second;vector<vector<long>>cc(32,vector<long>(136)),hc(512,vector<long>(136));for(auto&r:cc)for(auto&x:r)in>>x;for(auto&r:hc)for(auto&x:r)in>>x;
  vector<string>Oi(256);for(auto&x:Oi)in>>x;int ignored;for(int i=0;i<45;i++)in>>ignored;
  int n;in>>n;assert(n==19);vector<P>low;for(int i=0;i<n;i++)low.push_back(readqpoly(in));in>>n;assert(n==128);vector<P>grads;for(int i=0;i<n;i++)grads.push_back(readqpoly(in));in>>n;assert(n==8);vector<P>traces;for(int i=0;i<n;i++)traces.push_back(readqpoly(in));
  vector<int>selected24(67),selected28(108);for(auto&x:selected24)in>>x;for(auto&x:selected28)in>>x;string extra;assert(in&&!(in>>extra));
  vector<P>quad;for(auto[i,j]:mon)quad.push_back(mul(beta[i],beta[j]));M C(16),D(16),Theta(256),Hess(256);
  for(int row=0;row<32;row++)for(int k=0;k<136;k++)if(cc[row][k]){P t=scale(quad[k],cc[row][k],4);(row<16?C[row]:D[row-16])=add(row<16?C[row]:D[row-16],t);}
  for(int row=0;row<512;row++)for(int k=0;k<136;k++)if(hc[row][k]){P t=scale(quad[k],hc[row][k],2);(row<256?Theta[row]:Hess[row-256])=add(row<256?Theta[row]:Hess[row-256],t);}
  for(int i=0;i<8;i++){P e;for(int j=0;j<16;j++)e=add(e,mul(beta[j],grads[16*i+j]));assert(equal(e,scale(traces[i],i<5?8:12)));}
  auto poisson=[&](int f,int g){P p;for(int i=0;i<16;i++)for(int j=0;j<16;j++)if(Oi[16*i+j]!="0"&&!zero(grads[16*f+i])&&!zero(grads[16*g+j]))p=add(p,scaleq(mul(grads[16*f+i],grads[16*g+j]),Oi[16*i+j]));return p;};
  assert(equal(poisson(0,1),low[0])&&equal(poisson(0,2),low[1])&&equal(poisson(1,2),low[2]));
  M mu(16);for(int i=0;i<16;i++)mu[i]=add(scale(C[i],-1),scale(D[i],2));P F=trpair(mu,mu,4);assert(equal(F,low[9]));
  vector<P>pencil(5);vector<int>perm={0,1,2,3};
  do{int sign=1;for(int i=0;i<4;i++)for(int j=i+1;j<4;j++)if(perm[i]>perm[j])sign=-sign;vector<P>cur={scale(one28(),sign)};
   for(int i=0;i<4;i++){vector<P>next(cur.size()+1);for(int k=0;k<int(cur.size());k++){next[k]=add(next[k],mul(cur[k],C[4*i+perm[i]]));next[k+1]=add(next[k+1],mul(cur[k],D[4*i+perm[i]]));}cur=move(next);}for(int k=0;k<5;k++)pencil[k]=add(pencil[k],cur[k]);
  }while(next_permutation(perm.begin(),perm.end()));
  assert(equal(mul(F,F),add(scale(pencil[0],400),scale(pencil[4],1600))));
  P H6,H8;for(int column=0;column<16;column++){
   vector<P>level(16);for(int i=0;i<16;i++)level[i]=Hess[16*i+column];
   for(int degree=1;degree<=7;degree++){vector<P>next(16);for(int i=0;i<16;i++)for(int j=0;j<16;j++)if(!zero(Theta[16*i+j])&&!zero(level[j]))next[i]=add(next[i],mul(Theta[16*i+j],level[j]));level=move(next);if(degree==5)H6=add(H6,level[column]);if(degree==7)H8=add(H8,level[column]);}
  }assert(equal(H6,low[18]));
  vector<P>prim12={low[15],low[16],low[17],H6},q14={low[0],low[1],low[2]};
  vector<P>gens={F,pencil[0],pencil[1],pencil[2],pencil[3],pencil[4]};gens.insert(gens.end(),prim12.begin(),prim12.end());gens.push_back(H8);gens.insert(gens.end(),q14.begin(),q14.end());writepolys(dir+"/source-generators14.json",gens);
  vector<P>all24;
  for(int r=0;r<=3;r++){vector<int>indices;auto rec=[&](auto&&self,int next,int left)->void{if(!left){P p=power28(F,6-2*r);for(int j:indices)p=mul(p,pencil[j+1]);all24.push_back(p);return;}for(int j=next;j<4;j++){indices.push_back(j);self(self,j,left-1);indices.pop_back();}};rec(rec,0,r);}assert(all24.size()==35);
  vector<P>spectral12={power28(F,3)};for(int j=1;j<=4;j++)spectral12.push_back(mul(F,pencil[j]));
  for(auto&p:prim12)for(auto&q:spectral12)all24.push_back(mul(p,q));for(int i=0;i<4;i++)for(int j=i;j<4;j++)all24.push_back(mul(prim12[i],prim12[j]));all24.push_back(mul(H8,mul(F,F)));for(int j=1;j<=4;j++)all24.push_back(mul(H8,pencil[j]));assert(all24.size()==70);
  vector<pair<int,int>>pairs4,pairs3;for(int i=0;i<4;i++)for(int j=i;j<4;j++)pairs4.push_back({i,j});for(int i=0;i<3;i++)for(int j=i;j<3;j++)pairs3.push_back({i,j});
  map<E,vector<ulong>>coeffs;vector<long>lengths;fmpz_t zz;fmpz_init(zz);ulong ex[4];int column=0;
  for(int idx:selected28){P p;
   if(idx<67)p=mul(F,all24[selected24[idx]]);
   else if(idx<107){int k=idx-67;auto[i,j]=pairs4[k%10];p=mul(prim12[k/10],mul(pencil[i+1],pencil[j+1]));}
   else if(idx<111)p=mul(H8,prim12[idx-107]);
   else{auto[i,j]=pairs3[idx-111];p=mul(q14[i],q14[j]);}
   long len=fmpq_mpoly_length(p.v,ctx);assert(len<=3696&&len<=100000);lengths.push_back(len);ulong scalar=qmod(p.v->content);
   for(long k=0;k<len;k++){fmpq_mpoly_get_term_exp_ui(ex,p.v,k,ctx);assert(ex[1]%2==0&&(2*ex[0]+3*ex[1])%4==0&&ex[1]<=8&&ex[0]+ex[1]<=28&&ex[2]+ex[3]<=(28-3*ex[1])/2);E e={(2*ex[0]+3*ex[1])/4,ex[1]/2,ex[2],ex[3]};auto&v=coeffs[e];if(v.empty())v.resize(108);fmpz_mpoly_get_term_coeff_fmpz(zz,p.v->zpoly,k,ctx->zctx);v[column]=scalar*fmpz_fdiv_ui(zz,PRIME)%PRIME;}
   string path=dir+"/column"+to_string(column)+".json";ofstream f(path);f<<"{\"universal_candidate_index\":"<<idx<<",\"xi0_power\":14,\"HLsq_polynomial\":";writepoly(f,p);f<<"}\n";f.close();saved(path);cerr<<"COLUMN28 "<<column<<" CANDIDATE "<<idx<<" TERMS "<<len<<endl;column++;
  }
  map<int,vector<ulong>>basis;vector<E>rows;vector<vector<ulong>>rawrows;
  for(auto&[e,raw]:coeffs){vector<ulong>v=raw;for(auto&[pivot,b]:basis){ulong z=v[pivot];if(z)for(int j=pivot;j<108;j++)v[j]=(v[j]+PRIME-z*b[j]%PRIME)%PRIME;}int pivot=0;while(pivot<108&&!v[pivot])pivot++;if(pivot==108)continue;ulong inv=power(v[pivot],PRIME-2);for(int j=pivot;j<108;j++)v[j]=v[j]*inv%PRIME;basis[pivot]=move(v);rows.push_back(e);rawrows.push_back(raw);if(basis.size()==108)break;}
  auto dense=rawrows;int rank=0;ulong determinant=1;
  for(int col=0;col<108&&rank<int(dense.size());col++){int p=rank;while(p<int(dense.size())&&!dense[p][col])p++;if(p==int(dense.size()))continue;if(p!=rank){swap(dense[p],dense[rank]);determinant=(PRIME-determinant)%PRIME;}ulong z=dense[rank][col];determinant=determinant*z%PRIME;ulong inv=power(z,PRIME-2);for(int j=col;j<108;j++)dense[rank][j]=dense[rank][j]*inv%PRIME;for(int k=rank+1;k<int(dense.size());k++){ulong z=dense[k][col];if(z)for(int j=col;j<108;j++)dense[k][j]=(dense[k][j]+PRIME-z*dense[rank][j]%PRIME)%PRIME;}rank++;}assert(rank==int(basis.size()));
  ofstream f(dir+"/native-certificate.json");f<<"{\"all_passed\":true,\"prime\":32003,\"rank\":"<<rank<<",\"minor_mod32003\":"<<(rank==108?determinant:0)<<",\"selected_universal_candidates\":";ints(f,selected28);f<<",\"coefficient_monomials_r0_B0_s0_q0\":[";for(size_t i=0;i<rows.size();i++){if(i)f<<',';f<<'[';for(int j=0;j<4;j++){if(j)f<<',';f<<rows[i][j];}f<<']';}f<<"],\"selected_coefficient_matrix\":[";for(size_t i=0;i<rawrows.size();i++){if(i)f<<',';f<<'[';for(int j=0;j<108;j++){if(j)f<<',';f<<rawrows[i][j];}f<<']';}f<<"],\"column_term_counts\":[";for(size_t i=0;i<lengths.size();i++){if(i)f<<',';f<<lengths[i];}f<<"],\"ambient128_gradient_Euler_and14_brackets_verified\":true,\"physical_theta5H_matches_saved\":true,\"physical_theta7H_new\":true,\"all108_parity_and_support\":true,\"exactQQ_columns_saved\":108,\"max_terms\":"<<maxterms<<",\"products\":"<<products<<",\"bytes_saved\":"<<savedbytes<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"}\n";f.close();saved(dir+"/native-certificate.json");fmpz_clear(zz);cerr<<"SOURCE28_FINAL_RANK "<<rank<<" MINOR "<<determinant<<endl;
 }catch(const exception&e){cerr<<"BOUNDED_SOURCE28_STOP "<<e.what()<<endl;status=2;}
 }fmpq_mpoly_ctx_clear(ctx);flint_cleanup();return status;
}
