// Fixed56 row replacements and exact univariate gcd only. No rank search.
#define main frozen_q5_row_main_NOT_CALLED
#include "s61_nonsquare_P100_q5_pure98_row_20261104_v63.cpp"
#undef main
#include <sstream>
bool same(const P&a,const P&b){return fmpq_poly_equal(a.v,b.v);}
bool zz(const P&p){return fmpq_poly_is_zero(p.v);}
P subp(const P&a,const P&b){return plusp(a,negp(b));}
P readpoly(istream&f){long n;f>>n;P p;fmpq_t q;fmpq_init(q);for(long i=0;i<n;i++){long e;string c;f>>e>>c;assert(!fmpq_set_str(q,c.c_str(),10));fmpq_canonicalise(q);fmpq_poly_set_coeff_fmpq(p.v,e,q);}fmpq_clear(q);guard(p);return p;}
P divscalar(const P&p,const string&c){P q;fmpq_t s;fmpq_init(s);assert(!fmpq_set_str(s,c.c_str(),10));fmpq_canonicalise(s);assert(!fmpq_is_zero(s));fmpq_poly_scalar_div_fmpq(q.v,p.v,s);fmpq_clear(s);guard(q);return q;}
struct Split{P full,Q,U;string c="1";vector<int>e;};vector<P>F;
Split strip(const P&p){Split s;s.full=p;s.Q=p;s.U=one();s.e=vector<int>(5);if(zz(p))return s;
 for(int k=0;k<5;k++)for(int n=0;n<=1000;n++){assert(n<1000);P q,r;fmpq_poly_divrem(q.v,r.v,s.Q.v,F[k].v);guard(q);guard(r);if(!zz(r))break;assert(same(mulp(q,F[k]),s.Q));s.Q=q;s.e[k]++;s.U=mulp(s.U,F[k]);}
 fmpq_t c;fmpq_init(c);fmpq_poly_get_coeff_fmpq(c,s.Q.v,fmpq_poly_degree(s.Q.v));char*t=fmpq_get_str(nullptr,10,c);s.c=t;flint_free(t);fmpq_clear(c);s.Q=divscalar(s.Q,s.c);assert(same(scalar(mulp(s.U,s.Q),s.c),p));return s;}
struct CapBuf:streambuf{string data;streamsize xsputn(const char*s,streamsize n)override{if(data.size()+n>27*1024*1024)throw runtime_error("27MiB_native_JSON_guard");data.append(s,n);return n;}int overflow(int c)override{if(c!=EOF){char q=c;xsputn(&q,1);}return c;}};
void vecints(ostream&o,const vector<int>&v){o<<'[';for(int i=0;i<int(v.size());i++){if(i)o<<',';o<<v[i];}o<<']';}
void emit_split(ostream&o,const Split&s){o<<"{\"full\":";jp(o,s.full);o<<",\"primitive\":";jp(o,s.Q);o<<",\"unit_polynomial\":";jp(o,s.U);o<<",\"QQ_unit\":\""<<s.c<<"\",\"five_unit_multiplicities\":";vecints(o,s.e);o<<'}';}
int main(int argc,char**argv){assert(argc==3);began=chrono::steady_clock::now();string folder=argv[2];try{
 ifstream f(argv[1]);f>>input_reserve;P z0=shiftp(one(),1);F={z0,plusp(z0,pc("20")),plusp(z0,pc("60")),plusp(scalar(z0,"203"),pc("2700")),plusp(z0,pc("12"))};P D=readpoly(f);vector<vector<P>>S(8,vector<P>(8)),H=S;
 for(auto&row:S)for(auto&p:row)p=readpoly(f);for(auto&row:H)for(auto&p:row)p=readpoly(f);
 stage="old_selected_adjugate_binding";for(int i=0;i<8;i++)for(int j=0;j<8;j++){P l,r;for(int k=0;k<8;k++){l=plusp(l,mulp(S[i][k],H[k][j]));r=plusp(r,mulp(H[i][k],S[k][j]));}assert(same(l,i==j?D:P())&&same(r,i==j?D:P()));}
 vector<vector<int>>ex(7,vector<int>(5));for(auto&row:ex)for(auto&v:row)f>>v;vector<vector<P>>V(7,vector<P>(8));vector<P>LCD(7);
 stage="all56_new_LCD_source_clearings";for(int i=0;i<7;i++){P d=one();for(int k=0;k<5;k++)for(int j=0;j<ex[i][k];j++)d=mulp(d,F[k]);LCD[i]=d;
  for(int j=0;j<8;j++){R v=inputr(f);V[i][j]=mulp(v.n,quotient(d,v.d));assert(same(mulp(V[i][j],v.d),mulp(v.n,d)));}}
 string extra;assert(!(f>>extra));stage="oldD_knownunit_strip";vector<Split>splits;vector<P>weights;splits.push_back(strip(D));weights.push_back(one());P G=splits[0].Q;
 cout<<"OLD_D_DEGREE "<<fmpq_poly_degree(D.v)<<" UNKNOWN_PRIMITIVE_DEGREE "<<fmpq_poly_degree(G.v)<<'\n'<<flush;
 vector<pair<int,int>>labels;vector<int>gdegrees;int computed=0,reductions=0;
 for(int light=0;light<7&&fmpq_poly_degree(G.v)>0;light++)for(int replace=0;replace<8&&fmpq_poly_degree(G.v)>0;replace++){
  stage="fixed_minor_"+to_string(light)+"_"+to_string(replace);P p,q;for(int k=0;k<8;k++)p=plusp(p,mulp(V[light][k],H[k][replace]));for(int k=7;k>=0;k--)q=plusp(q,mulp(V[light][k],H[k][replace]));assert(same(p,q));
  splits.push_back(strip(p));labels.push_back({light,replace});computed++;weights.emplace_back();
  if(!zz(p)){
   P g;fmpq_poly_gcd(g.v,G.v,splits.back().Q.v);guard(g);
   if(!same(g,G)){P x,y,h;fmpq_poly_xgcd(h.v,x.v,y.v,G.v,splits.back().Q.v);guard(h);guard(x);guard(y);assert(same(h,g)&&same(plusp(mulp(x,G),mulp(y,splits.back().Q)),g));
    for(int k=0;k<int(weights.size())-1;k++)weights[k]=mulp(weights[k],x);weights.back()=y;G=g;reductions++;
   }
  }
  P check;for(int k=0;k<int(weights.size());k++)check=plusp(check,mulp(weights[k],splits[k].Q));assert(same(check,G));gdegrees.push_back(fmpq_poly_degree(G.v));
 }
 stage="literal_final_gcd_and_actual_minor_Bezout";P check;for(int k=0;k<int(weights.size());k++){check=plusp(check,mulp(weights[k],splits[k].Q));if(!zz(splits[k].Q))quotient(splits[k].Q,G);}assert(same(check,G));
 vector<int>Epow(5);vector<int>active;for(int i=0;i<int(weights.size());i++)if(!zz(weights[i])){active.push_back(i);for(int j=0;j<5;j++)Epow[j]=max(Epow[j],splits[i].e[j]);}
 P E=one();for(int j=0;j<5;j++)for(int k=0;k<Epow[j];k++)E=mulp(E,F[j]);vector<P>actualcof(weights.size());P actualcheck;
 for(int i:active){actualcof[i]=divscalar(mulp(weights[i],quotient(E,splits[i].U)),splits[i].c);actualcheck=plusp(actualcheck,mulp(actualcof[i],splits[i].full));}assert(same(actualcheck,mulp(E,G)));
 CapBuf buffer;ostream o(&buffer);o.exceptions(ios::badbit|ios::failbit);o<<"{\"all_passed\":true,\"old_minor_strip\":";emit_split(o,splits[0]);o<<",\"computed_fixed_replacements\":"<<computed<<",\"stopped_at_knownunit_gcd\":"<<(fmpq_poly_degree(G.v)==0?"true":"false")<<",\"candidates\":[";
 for(int i=0;i<computed;i++){if(i)o<<',';o<<"{\"K2_light_position\":"<<labels[i].first<<",\"replaced_old_selected_row_position\":"<<labels[i].second<<",\"exact_forward_reverse_8dot_verified\":true,\"split\":";emit_split(o,splits[i+1]);o<<'}';}
 o<<"],\"gcd_degree_after_each_candidate\":";vecints(o,gdegrees);o<<",\"final_gcd\":";jp(o,G);o<<",\"literal_primitive_Bezout_contributors\":[";
 for(int i=0;i<int(active.size());i++){if(i)o<<',';int k=active[i];o<<"{\"basis_index\":"<<k<<",\"coefficient\":";jp(o,weights[k]);o<<'}';}
 o<<"],\"common_actual_minor_unit_powers5\":";vecints(o,Epow);o<<",\"common_actual_minor_unit_polynomial\":";jp(o,E);o<<",\"actual_minor_polynomial_Bezout_contributors\":[";
 for(int i=0;i<int(active.size());i++){if(i)o<<',';int k=active[i];o<<"{\"basis_index\":"<<k<<",\"coefficient\":";jp(o,actualcof[k]);o<<'}';}
 o<<"],\"literal_sum_Ai_actualMinor_i_equals_E_times_gcd_verified\":true,\"all56_selected_LCD_clearings_verified\":true,\"all128_old_adjugate_products_verified\":true,\"unknown_gcd_factors_NOT_inverted\":true,\"gcd_reduction_count\":"<<reductions<<",\"peak_live_polynomial_terms_with_literal_reserve\":"<<peak_live+input_reserve<<",\"literal_input_reserve\":"<<input_reserve<<",\"serialized_terms\":"<<stored<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";
 ofstream of(folder+"/native.json");of<<buffer.data;
 }catch(const exception&e){ofstream o(folder+"/frontier.json");o<<"{\"all_passed\":false,\"stage\":\""<<stage<<"\",\"reason\":\""<<e.what()<<"\"}\n";return 2;}flint_cleanup();}
