// Structural repair on the SAME formal plane; old H is a 16x16 Hessian.
// Reuse immutable polynomial/rank helpers, not the old main or its computations.
#define main unused_cd_only_main
#include "q_new_family_full_degree24_universal_basis_20261104_v63.cpp"
#undef main
using Const=vector<ulong>;
M matmul(const M&a,int n,int k,const M&b,int m){M c(n*m);for(int i=0;i<n;i++)for(int j=0;j<m;j++)for(int l=0;l<k;l++)c[i*m+j]=add(c[i*m+j],mul(a[i*k+l],b[l*m+j]));return c;}
M constleft(const Const&a,int n,int k,const M&b,int m){M c(n*m);for(int i=0;i<n;i++)for(int j=0;j<m;j++)for(int l=0;l<k;l++)if(a[i*k+l])c[i*m+j]=add(c[i*m+j],scale(a[i*k+l],b[l*m+j]));return c;}
M transp(const M&a,int n,int m){M b(n*m);for(int i=0;i<n;i++)for(int j=0;j<m;j++)b[j*n+i]=a[i*m+j];return b;}
P tracepair16(const M&a,const M&b){P p;for(int i=0;i<16;i++)for(int j=0;j<16;j++)p=add(p,mul(a[16*i+j],b[16*j+i]));return p;}
P readcoef(istream&in,int degree){P p;ulong e[2];for(int i=0;i<=degree;i++)for(int j=0;j<=degree-i;j++){ulong z;in>>z;assert(z<PRIME);e[0]=i;e[1]=j;if(z)nmod_mpoly_set_coeff_ui_ui(p.v,z,e,ctx);}assert(in);return p;}
Const readconst(istream&in,int n){Const a(n);for(auto&z:a){in>>z;assert(z<PRIME);}assert(in);return a;}
int main(int argc,char**argv){assert(argc==3);nmod_mpoly_ctx_init(ctx,2,ORD_DEGLEX,PRIME);int result=0;auto start=chrono::steady_clock::now();{
 try{
  ifstream in(argv[1]);int dim;ulong prime;in>>dim>>prime;assert(dim==16&&prime==PRIME);vector<P>beta(16);ulong e[2];for(int j=0;j<3;j++)for(int i=0;i<16;i++){ulong z;in>>z;e[0]=j==1;e[1]=j==2;nmod_mpoly_set_coeff_ui_ui(beta[i].v,z,e,ctx);}
  Const O=readconst(in,256),Oi=readconst(in,256),Gi=readconst(in,100);vector<Const>rep;for(int j=0;j<10;j++)rep.push_back(readconst(in,256));
  P F=readcoef(in,4);vector<P>E;for(int i=0;i<5;i++)E.push_back(readcoef(in,8));vector<P>cd12,cd8;for(int i=0;i<20;i++)cd12.push_back(readcoef(in,12));for(int i=0;i<36;i++)cd8.push_back(readcoef(in,16));string extra;assert(!(in>>extra));
  M V(160);for(int j=0;j<10;j++){M rb=constleft(rep[j],16,16,beta,1),v=constleft(O,16,16,rb,1);for(int i=0;i<16;i++)V[10*i+j]=v[i];}
  M mu(10);for(int j=0;j<10;j++)for(int i=0;i<16;i++)mu[j]=add(mu[j],mul(beta[i],V[10*i+j]));for(auto&p:mu)p=scale(modpow(2,PRIME-2),p);
  M co=constleft(Gi,10,10,mu,1),theta(256);for(int k=0;k<256;k++)for(int j=0;j<10;j++)if(rep[j][k])theta[k]=add(theta[k],scale(rep[j][k],co[j]));
  M left=constleft(Oi,16,16,V,10),right=constleft(Gi,10,10,transp(V,16,10),16),H=matmul(left,16,10,right,16);for(int k=0;k<256;k++)H[k]=add(H[k],theta[k]);
  for(const auto&a:{theta,H}){M oa=constleft(O,16,16,a,16);for(int i=0;i<16;i++)for(int j=0;j<16;j++)assert(nmod_mpoly_equal(oa[16*i+j].v,oa[16*j+i].v,ctx));}
  map<string,M>table;table["0"]=theta;table["1"]=H;
  auto prefix=[&](int len){for(int k=0;k<(1<<len);k++){string w;for(int j=len-1;j>=0;j--)w+=char('0'+((k>>j)&1));table[w]=matmul(table.at(w.substr(0,len-1)),16,16,w.back()=='0'?theta:H,16);}};
  prefix(2);prefix(3);auto w6=necklaces(6),w8=necklaces(8);vector<P>old6,old8;for(auto&w:w6)old6.push_back(tracepair16(table.at(w.substr(0,3)),table.at(w.substr(3))));
  vector<P>deg12=cd12;for(auto&p:old6)deg12.push_back(p);Rank r12;for(int i=0;i<(int)deg12.size();i++)r12.push(coeffs(deg12[i],12),i);
  bool gate12=r12.rank()==9&&r12.selected.size()==9&&vector<int>(r12.selected.begin(),r12.selected.begin()+8)==vector<int>({0,1,2,3,4,5,9,10});
  vector<P>basis12,prim12,prod16,deg16,all24;Rank r16prod,r16,r24;bool gate16prod=false,gate16=false,gate24=false,old8attempted=false;int chosen_cd8=-1,chosen_old8=-1;
  if(gate12){for(int i:r12.selected)basis12.push_back(deg12[i]);for(int i=5;i<9;i++)prim12.push_back(basis12[i]);for(auto&p:basis12)prod16.push_back(mul(F,p));for(int i=0;i<5;i++)for(int j=i;j<5;j++)prod16.push_back(mul(E[i],E[j]));
   assert(prod16.size()==24);for(int i=0;i<24;i++)r16prod.push(coeffs(prod16[i],16),i);gate16prod=r16prod.rank()==19;
   if(gate16prod){deg16=prod16;for(auto&p:cd8)deg16.push_back(p);for(int i=0;i<60;i++)r16.push(coeffs(deg16[i],16),i);
    if(r16.rank()==19){old8attempted=true;prefix(4);for(auto&w:w8)old8.push_back(tracepair16(table.at(w.substr(0,4)),table.at(w.substr(4))));for(auto&p:old8)deg16.push_back(p);for(int i=60;i<96;i++)r16.push(coeffs(deg16[i],16),i);}
    gate16=r16.rank()==20;
    if(gate16){int selected=r16.selected.back();P prim16;if(selected<60){chosen_cd8=selected-24;assert(chosen_cd8>=0);prim16=cd8[chosen_cd8];}else{chosen_old8=selected-60;assert(chosen_old8>=0);prim16=old8[chosen_old8];}
     for(int r=0;r<=3;r++){vector<int>indices;auto rec=[&](auto&&self,int next,int left)->void{if(!left){P q=power(F,6-2*r);for(int k:indices)q=mul(q,E[k+1]);all24.push_back(q);return;}for(int k=next;k<4;k++){indices.push_back(k);self(self,k,left-1);indices.pop_back();}};rec(rec,0,r);}assert(all24.size()==35);
     vector<P>spectral12={power(F,3)};for(int i=1;i<5;i++)spectral12.push_back(mul(F,E[i]));for(auto&q:prim12)for(auto&s:spectral12)all24.push_back(mul(q,s));for(int i=0;i<4;i++)for(int j=i;j<4;j++)all24.push_back(mul(prim12[i],prim12[j]));for(auto&q:E)all24.push_back(mul(prim16,q));assert(all24.size()==70);
     for(int i=0;i<70;i++)r24.push(coeffs(all24[i],24),i);gate24=r24.rank()==67;
    }
   }
  }
  string status=!gate12?"degree12_repair_gate_failed":!gate16prod?"degree16_product_gate_failed":!gate16?"degree16_repair_gate_failed":!gate24?"degree24_rank_gate_failed":"complete_Inv24_basis_certified";
  ofstream out(argv[2]);out<<"{\"all_passed\":true,\"status\":\""<<status<<"\",\"prime\":32003,\"coefficient_order\":\"i=0..d,j=0..d-i: s^i*t^j\",\"H_is16x16_Hessian_not4x4_D\":true,\"old_H_symplectic_matrix_identity_verified\":true,\"old_trace8_attempted\":"<<(old8attempted?"true":"false")<<",\"old_trace6_columns\":";polys(out,old6,12);out<<",\"old_trace8_columns\":";polys(out,old8,16);
  out<<",\"degree12_columns\":";polys(out,deg12,12);out<<",\"degree12_rank\":"<<r12.rank()<<",\"degree12_selected_columns\":";ints(out,r12.selected);out<<",\"degree12_selected_coefficient_rows\":";ints(out,r12.pivots);
  out<<",\"degree16_product_rank\":"<<r16prod.rank()<<",\"degree16_columns\":";polys(out,deg16,16);out<<",\"degree16_rank\":"<<r16.rank()<<",\"degree16_selected_columns\":";ints(out,r16.selected);out<<",\"degree16_selected_coefficient_rows\":";ints(out,r16.pivots);out<<",\"primitive16_CD_trace8_index\":"<<chosen_cd8<<",\"primitive16_old_trace8_index\":"<<chosen_old8;
  out<<",\"degree24_columns\":";polys(out,all24,24);out<<",\"degree24_rank\":"<<r24.rank()<<",\"degree24_selected_columns\":";ints(out,r24.selected);out<<",\"degree24_selected_coefficient_rows\":";ints(out,r24.pivots);out<<",\"max_intermediate_terms\":"<<maxterms<<",\"polynomial_multiplications\":"<<multiplications<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";out.close();cerr<<status<<" RANKS "<<r12.rank()<<' '<<r16prod.rank()<<' '<<r16.rank()<<' '<<r24.rank()<<" OLD8 "<<old8attempted<<endl;
 }catch(const exception&e){cerr<<"BOUNDED_REPAIR_STOP "<<e.what()<<endl;result=2;}
 }nmod_mpoly_ctx_clear(ctx);flint_cleanup();return result;
}
