// ONE universal degree26 same-plane basis job. ALL gradients are AMBIENT16.
#define main unused_immutable_CD_only_main
#include "q_new_family_full_degree24_universal_basis_20261104_v63.cpp"
#undef main
#include <set>
P mul26(const P&a,const P&b){P c;nmod_mpoly_mul(c.v,a.v,b.v,ctx);++multiplications;long n=nmod_mpoly_length(c.v,ctx);maxterms=max(maxterms,n);if(n>378)throw runtime_error("formal degree26 support exceeded378");return c;}
P readco(istream&in,int degree){P p;ulong e[2];for(int i=0;i<=degree;i++)for(int j=0;j<=degree-i;j++){ulong z;in>>z;assert(z<PRIME);e[0]=i;e[1]=j;if(z)nmod_mpoly_set_coeff_ui_ui(p.v,z,e,ctx);}return p;}
using Grad=vector<P>;
int main(int argc,char**argv){assert(argc==3);nmod_mpoly_ctx_init(ctx,2,ORD_DEGLEX,PRIME);int result=0;auto started=chrono::steady_clock::now();{
 try{
  ifstream in(argv[1]);int dim,width;ulong prime;in>>dim>>width>>prime;assert(dim==16&&width==136&&prime==PRIME);vector<P>b(16);ulong exponent[2];for(int j=0;j<3;j++)for(int i=0;i<16;i++){ulong z;in>>z;assert(z<PRIME);exponent[0]=j==1;exponent[1]=j==2;nmod_mpoly_set_coeff_ui_ui(b[i].v,z,exponent,ctx);}
  vector<pair<int,int>>mon(width);for(auto&ij:mon)in>>ij.first>>ij.second;vector<vector<long>>cc(32,vector<long>(width));for(auto&r:cc)for(auto&z:r)in>>z;vector<ulong>Oi(256);for(auto&z:Oi)in>>z;for(int i=0;i<16;i++)for(int j=0;j<16;j++)assert((Oi[16*i+j]+Oi[16*j+i])%PRIME==0);
  P F4=readco(in,4);vector<P>Inv12;for(int i=0;i<9;i++)Inv12.push_back(readco(in,12));string extra;assert(in&&!(in>>extra));
  M C(16),D(16);vector<vector<M>>derivatives(2,vector<M>(16,M(16)));ulong inv4=modpow(4,PRIME-2);
  // First differentiate the136 ambient beta_i*beta_j monomials in EACH beta_k.
  // The b_i polynomial substitution is performed only when evaluating that derivative.
  for(int r=0;r<32;r++)for(int k=0;k<136;k++)if(cc[r][k]){auto[i,j]=mon[k];P term=scale(cc[r][k],mul(b[i],b[j]));int letter=r/16,pos=r%16;(letter==0?C[pos]:D[pos])=add(letter==0?C[pos]:D[pos],term);
   derivatives[letter][i][pos]=add(derivatives[letter][i][pos],scale(cc[r][k],b[j]));
   derivatives[letter][j][pos]=add(derivatives[letter][j][pos],scale(cc[r][k],b[i]));
  }
  for(auto&p:C)p=scale(inv4,p);for(auto&p:D)p=scale(inv4,p);for(auto&a:derivatives)for(auto&k:a)for(auto&p:k)p=scale(inv4,p);
  for(int letter=0;letter<2;letter++)for(int pos=0;pos<16;pos++){P euler;for(int k=0;k<16;k++)euler=add(euler,mul(b[k],derivatives[letter][k][pos]));assert(nmod_mpoly_equal(euler.v,scale(2,letter==0?C[pos]:D[pos]).v,ctx));}
  vector<string>w4={"0000","0001","0011","0101","0111"},w6={"000101","000111","011011"},allwords=w4;allwords.insert(allwords.end(),w6.begin(),w6.end());set<string>prefixes;
  for(auto&w:allwords){for(int r=0;r<int(w.size());r++){string rest=w.substr(r+1)+w.substr(0,r);for(int k=1;k<=int(rest.size());k++)prefixes.insert(rest.substr(0,k));}for(int half=0;half<2;half++){string h=w.substr(half*w.size()/2,w.size()/2);for(int k=1;k<=int(h.size());k++)prefixes.insert(h.substr(0,k));}}
  vector<string>order(prefixes.begin(),prefixes.end());sort(order.begin(),order.end(),[](auto&a,auto&b){return a.size()==b.size()?a<b:a.size()<b.size();});map<string,M>table;table["0"]=C;table["1"]=D;for(auto&w:order)if(w.size()>1)table[w]=mm(table.at(w.substr(0,w.size()-1)),w.back()=='0'?C:D);
  vector<P>traces;vector<Grad>grads;
  for(auto&w:allwords){int n=w.size();traces.push_back(trpair(table.at(w.substr(0,n/2)),table.at(w.substr(n/2))));vector<M>partial(2,M(16));for(int r=0;r<n;r++){string rest=w.substr(r+1)+w.substr(0,r);for(int p=0;p<16;p++)partial[w[r]-'0'][p]=add(partial[w[r]-'0'][p],table.at(rest)[p]);}
   Grad grad(16);for(int k=0;k<16;k++)for(int letter=0;letter<2;letter++)for(int i=0;i<4;i++)for(int j=0;j<4;j++)grad[k]=add(grad[k],mul(partial[letter][4*i+j],derivatives[letter][k][4*j+i]));
   P euler;for(int k=0;k<16;k++)euler=add(euler,mul(b[k],grad[k]));assert(nmod_mpoly_equal(euler.v,scale(2*n,traces.back()).v,ctx));grads.push_back(move(grad));
  }
  auto poisson=[&](int f,int g){P p;for(int i=0;i<16;i++)for(int j=0;j<16;j++)if(Oi[16*i+j])p=add(p,scale(Oi[16*i+j],mul(grads[f][i],grads[g][j])));return p;};
  vector<P>Inv8(traces.begin(),traces.begin()+5),degree14;vector<pair<int,int>>pairs14;for(int i=0;i<5;i++)for(int j=i+1;j<5;j++){pairs14.push_back({i,j});degree14.push_back(poisson(i,j));}
  Rank r14;for(int i=0;i<10;i++)r14.push(coeffs(degree14[i],14),i);bool gate14=r14.rank()==3;vector<P>basis14,degree18,degree26;Rank r18lower,r18,r26;bool gate18lower=false,gate18=false,gate26=false;vector<int>primitive18;
  if(gate14){for(int i:r14.selected)basis14.push_back(degree14[i]);for(auto&p:basis14)degree18.push_back(mul(F4,p));for(int i=0;i<3;i++)r18lower.push(coeffs(degree18[i],18),i);gate18lower=r18lower.rank()==3;
   if(gate18lower){for(int i=0;i<5;i++)for(int j=5;j<8;j++)degree18.push_back(poisson(i,j));for(int i=0;i<18;i++)r18.push(coeffs(degree18[i],18),i);gate18=r18.rank()==9&&r18.selected.size()==9&&vector<int>(r18.selected.begin(),r18.selected.begin()+3)==vector<int>({0,1,2});
    if(gate18){primitive18.assign(r18.selected.begin()+3,r18.selected.end());for(auto&p:basis14)for(auto&q:Inv12)degree26.push_back(mul26(p,q));for(int i:primitive18)for(auto&q:Inv8)degree26.push_back(mul26(degree18[i],q));assert(degree26.size()==57);for(int i=0;i<57;i++)r26.push(coeffs(degree26[i],26),i);gate26=r26.rank()==45;}
   }
  }
  string status=!gate14?"Inv14_rank_gate_failed":!gate18lower?"F4_Inv14_rank_gate_failed":!gate18?"Inv18_rank_gate_failed":!gate26?"Inv26_rank_gate_failed":"complete_QQ_Inv26_basis_certified";
  ofstream out(argv[2]);out<<"{\"all_passed\":true,\"status\":\""<<status<<"\",\"prime\":32003,\"ambient16_derivatives_before_substitution\":true,\"Poisson_convention\":\"gradient(f)^T OmegaInverse gradient(g)\",\"all32_quadratic_and_all8_trace_Euler_identities\":true,\"coefficient_order\":\"i=0..d,j=0..d-i: s^i*t^j\",\"Inv8_trace4_columns\":";polys(out,Inv8,8);out<<",\"trace6_columns\":";polys(out,vector<P>(traces.begin()+5,traces.end()),12);
  out<<",\"ambient16_gradient_columns\":[";for(int i=0;i<8;i++){if(i)out<<',';polys(out,grads[i],i<5?7:11);}out<<']';
  out<<",\"degree14_columns\":";polys(out,degree14,14);out<<",\"degree14_rank\":"<<r14.rank()<<",\"degree14_selected_columns\":";ints(out,r14.selected);out<<",\"degree14_selected_coefficient_rows\":";ints(out,r14.pivots);
  out<<",\"degree18_F4_product_rank\":"<<r18lower.rank()<<",\"degree18_columns\":";polys(out,degree18,18);out<<",\"degree18_rank\":"<<r18.rank()<<",\"degree18_selected_columns\":";ints(out,r18.selected);out<<",\"degree18_selected_coefficient_rows\":";ints(out,r18.pivots);out<<",\"primitive18_selected_candidate_columns\":";ints(out,primitive18);
  out<<",\"degree26_columns\":";polys(out,degree26,26);out<<",\"degree26_rank\":"<<r26.rank()<<",\"degree26_selected_columns\":";ints(out,r26.selected);out<<",\"degree26_selected_coefficient_rows\":";ints(out,r26.pivots);out<<",\"max_intermediate_terms\":"<<maxterms<<",\"polynomial_multiplications\":"<<multiplications<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"}\n";cerr<<status<<" RANKS "<<r14.rank()<<' '<<r18.rank()<<' '<<r26.rank()<<endl;
 }catch(const exception&e){cerr<<"BOUNDED_UNIVERSAL26_STOP "<<e.what()<<endl;result=2;}
 }nmod_mpoly_ctx_clear(ctx);flint_cleanup();return result;
}
