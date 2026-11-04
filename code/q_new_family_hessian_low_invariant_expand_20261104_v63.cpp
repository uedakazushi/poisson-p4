// ONE capped four-output source expansion; T^2=B torus preprocessing is exact.
#define main unused_original_low_expansion_main
#include "q_new_family_low_invariant_expand_20261104_v63.cpp"
#undef main
using Const=vector<long>;
M matmul(const M&a,int n,int k,const M&b,int m){M c(n*m);for(int i=0;i<n;i++)for(int j=0;j<m;j++)for(int l=0;l<k;l++)if(!fmpz_mpoly_is_zero(a[i*k+l].v,ctx)&&!fmpz_mpoly_is_zero(b[l*m+j].v,ctx))c[i*m+j]=add(c[i*m+j],mul(a[i*k+l],b[l*m+j]));return c;}
M cleft(const Const&a,int n,int k,const M&b,int m){M c(n*m);for(int i=0;i<n;i++)for(int j=0;j<m;j++)for(int l=0;l<k;l++)if(a[i*k+l])c[i*m+j]=add(c[i*m+j],scale(a[i*k+l],b[l*m+j]));return c;}
P powp(P a,int n){P b;fmpz_mpoly_one(b.v,ctx);while(n){if(n&1)b=mul(b,a);n>>=1;if(n)a=mul(a,a);}return b;}
P divide2(const P&p){P q;fmpz_t two;fmpz_init_set_ui(two,2);assert(fmpz_mpoly_scalar_divides_fmpz(q.v,p.v,two,ctx));fmpz_clear(two);return q;}
P determinant4(const M&a){P p;vector<int>perm={0,1,2,3};do{int sign=1;for(int i=0;i<4;i++)for(int j=i+1;j<4;j++)if(perm[i]>perm[j])sign=-sign;P term;fmpz_mpoly_set_si(term.v,sign,ctx);for(int i=0;i<4;i++)term=mul(term,a[4*i+perm[i]]);p=add(p,term);}while(next_permutation(perm.begin(),perm.end()));return p;}
P contract(const M&a,const M&b){assert(a.size()==160&&b.size()==160);P p;for(int i=0;i<160;i++)if(!fmpz_mpoly_is_zero(a[i].v,ctx)&&!fmpz_mpoly_is_zero(b[i].v,ctx))p=add(p,mul(a[i],b[i]));return p;}
string save_torus(const string&dir,const string&name,const P&p,int degree,long scalar_factor,int cpower){
 long n=fmpz_mpoly_length(p.v,ctx);assert(n);if(n>OUTPUT_TERMS)throw runtime_error("100k output term cap");
 fmpz_t den,g,z;fmpz_init(den);fmpz_init(g);fmpz_init(z);fmpz_set_ui(den,1080);fmpz_pow_ui(den,den,degree);fmpz_mul_ui(den,den,scalar_factor);fmpz_set(g,den);ulong e[4];long cancelB=degree/2,bits=0;
 for(long k=0;k<n;k++){fmpz_mpoly_get_term_coeff_fmpz(z,p.v,k,ctx);fmpz_gcd(g,g,z);fmpz_mpoly_get_term_exp_ui(e,p.v,k,ctx);assert(e[3]%2==0&&(e[0]+e[1])%2==0);cancelB=min(cancelB,long(e[3]/2));}
 fmpz_divexact(den,den,g);ofstream out(dir+"/"+name+".json");ostringstream head;
 head<<"{\"name\":\""<<name<<"\",\"variables\":[\"a=u/X\",\"Y=X^2\",\"W\",\"B=T^2\"],\"source_c_power\":"<<cpower
 <<",\"denominator_integer\":\""<<zstr(den)<<"\",\"denominator_B_power\":"<<degree/2-cancelB<<",\"denominator_D_power\":"<<degree
 <<",\"D\":\"147456*Y-(3072*W-1)^2\",\"T_squared_equals_B_preconditioned\":true,\"terms\":[";writebytes(out,head.str());
 for(long k=0;k<n;k++){fmpz_mpoly_get_term_coeff_fmpz(z,p.v,k,ctx);fmpz_divexact(z,z,g);bits=max(bits,long(fmpz_bits(z)));fmpz_mpoly_get_term_exp_ui(e,p.v,k,ctx);ostringstream row;if(k)row<<',';row<<"[\""<<zstr(z)<<"\","<<e[0]<<','<<(e[0]+e[1])/2<<','<<e[2]<<','<<e[3]/2-cancelB<<']';writebytes(out,row.str());}
 writebytes(out,"]}\n");out.close();ostringstream meta;meta<<"{\"name\":\""<<name<<"\",\"terms\":"<<n<<",\"max_coefficient_bits\":"<<bits<<",\"additional_B_cancelled\":"<<cancelB<<",\"denominator_B_power\":"<<degree/2-cancelB<<",\"denominator_integer\":\""<<zstr(den)<<"\"}";
 cerr<<"SAVED "<<name<<" TERMS "<<n<<" BITS "<<bits<<" TOTAL_RAW_BYTES "<<output_bytes<<endl;fmpz_clear(den);fmpz_clear(g);fmpz_clear(z);return meta.str();
}
int main(int argc,char**argv){assert(argc==3);fmpz_mpoly_ctx_init(ctx,4,ORD_DEGLEX);int result=0;auto start=chrono::steady_clock::now();{
 try{
  ifstream in(argv[1]);int vars,count,width;long L;in>>vars>>count>>width>>L;assert(vars==4&&count==16&&width==136&&L==1080);vector<P>b(16);for(auto&p:b)readpoly(in,p);vector<pair<int,int>>mon(width);for(auto&ij:mon)in>>ij.first>>ij.second;
  vector<vector<long>>coeff(32,vector<long>(width));for(auto&r:coeff)for(auto&z:r)in>>z;
  Const O(256),GHat(100);vector<Const>rho(10,Const(256));for(auto&z:O)in>>z;for(auto&z:GHat)in>>z;for(auto&r:rho)for(auto&z:r)in>>z;string extra;assert(in&&!(in>>extra));
  vector<P>bp;for(auto ij:mon)bp.push_back(mul(b[ij.first],b[ij.second]));M C(16),D(16);for(int r=0;r<32;r++){P sum;for(int j=0;j<width;j++)if(coeff[r][j])sum=add(sum,scale(coeff[r][j],bp[j]));(r<16?C[r]:D[r-16])=sum;}bp.clear();
  vector<string>stats;
  {M C2=mm(C,C),C3=mm(C2,C),D2=mm(D,D),D3=mm(D2,D),DC=mm(D,C),DCD=mm(DC,D);stats.push_back(save_torus(argv[2],"CD000101",trpair(C3,DCD),12,4096,6));stats.push_back(save_torus(argv[2],"CD000111",trpair(C3,D3),12,4096,6));}
  // Mhat=2mu(b) is integral, mu^T=-C_original+2D_original, C,D here are4times.
  M Mhat(16);for(int i=0;i<4;i++)for(int j=0;j<4;j++)Mhat[4*i+j]=divide2(add(scale(-1,C[4*j+i]),scale(2,D[4*j+i])));
  M S(16);for(int j=0;j<4;j++){S[j]=scale(-1,Mhat[4+j]);S[4+j]=Mhat[j];S[8+j]=scale(-1,Mhat[12+j]);S[12+j]=Mhat[8+j];}for(int i=0;i<4;i++)for(int j=0;j<4;j++)assert(fmpz_mpoly_equal(S[4*i+j].v,S[4*j+i].v,ctx));
  M Theta(256);int index=0;for(int i=0;i<4;i++)for(int j=i;j<4;j++){for(int k=0;k<256;k++)if(rho[index][k])Theta[k]=add(Theta[k],scale(rho[index][k],S[4*i+j]));index++;}assert(index==10);
  M Z(160);for(int j=0;j<10;j++){M v=cleft(rho[j],16,16,b,1);for(int i=0;i<16;i++)Z[10*i+j]=v[i];}M V=cleft(O,16,16,Z,10),VG(160);for(int i=0;i<16;i++)for(int a=0;a<10;a++)for(int j=0;j<10;j++)if(GHat[10*a+j])VG[10*i+a]=add(VG[10*i+a],scale(GHat[10*a+j],V[10*i+j]));
  P s1=divide2(tr(mm(Mhat,Mhat))),s2=determinant4(Mhat);P s12=mul(s1,s1),s13=mul(s12,s1);
  P pure6=add(scale(264,s13),scale(408,mul(s1,s2)));
  P pure8=add(add(scale(1032,mul(s12,s12)),scale(3488,mul(s12,s2))),scale(-4208,mul(s2,s2)));
  C.clear();D.clear();Mhat.clear();S.clear();b.clear();V.clear();
  M level=move(Z);
  for(int k=1;k<=7;k++){
   M next=matmul(Theta,16,16,level,10);level=move(next);cerr<<"VECTOR_LEVEL "<<k<<" MAX_TERMS "<<max_terms<<endl;
   if(k==5)stats.push_back(save_torus(argv[2],"theta5H",add(pure6,contract(level,VG)),12,64,6));
   if(k==7)stats.push_back(save_torus(argv[2],"theta7H",add(pure8,contract(level,VG)),16,256,8));
  }
  ofstream report(string(argv[2])+"/native-report.json");report<<"{\"all_passed\":true,\"outputs\":[";for(int i=0;i<int(stats.size());i++){if(i)report<<',';report<<stats[i];}report<<"],\"max_intermediate_terms\":"<<max_terms<<",\"polynomial_multiplications\":"<<multiplications<<",\"coefficient_output_bytes\":"<<output_bytes<<",\"all_T_parity_and_u_X_parity_verified\":true,\"integer_scale_theta\":2,\"vector_recursion_levels\":7,\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";cerr<<"ALL_FOUR_HESSIAN_SOURCE_OUTPUTS_PASS\n";
 }catch(const exception&e){cerr<<"BOUNDED_SOURCE_STOP "<<e.what()<<endl;result=2;}
 }fmpz_mpoly_ctx_clear(ctx);flint_cleanup();return result;
}
