// One cached literal K7 pure cubic, physicalW product only; no source query or rank.
// The renamed historical main is never called.
#define main frozen_q5_row_main_NOT_CALLED
#include "s61_nonsquare_P100_q5_pure98_row_20261104_v63.cpp"
#undef main
bool knowncurveunit(const R&r){
 if(zero(r))return false;P p=r.n;
 for(int step=0;step<100;step++){P g;fmpq_poly_gcd(g.v,p.v,unitproduct.v);guard(g);if(fmpq_poly_degree(g.v)==0)break;p=quotient(p,g);}
 return fmpq_poly_degree(p.v)==0;
}
int main(int argc,char**argv){assert(argc==3);began=chrono::steady_clock::now();string folder=argv[2];try{
 ifstream f(argv[1]);f>>input_reserve;assert(input_reserve<=8746);
 z=shiftp(one(),1);P zp20=plusp(z,pc("20"));ap={one(),scalar(zp20,"-10/9")};ep={one(),scalar(plusp(z,pc("100")),"1/9")};dp={one(),scalar(mulp(z,zp20),"-3/100")};
 R aa(ap[1],shiftp(one(),2)),ee(ep[1],z),dd(dp[1]);R CC=sub(sub(add(mul(rc("100"),aa),mul(rc("18"),powr(ee,2))),mul(rc("9"),ee)),rc("2"));
 unitproduct=mulp(mulp(mulp(mulp(z,zp20),plusp(z,pc("60"))),plusp(scalar(z,"203"),pc("2700"))),plusp(z,pc("12")));
 vector<R>oldpi(84);for(auto&p:oldpi)p=inputr(f);vector<Mon>oldmon(84);for(auto&m:oldmon)for(auto&i:m)f>>i;
 vector<vector<R>>Bphys(7,vector<R>(7));for(auto&row:Bphys)for(auto&v:row)v=inputr(f);long k;f>>k;map<long,AP>na;for(int j=0;j<k;j++){long id;f>>id;AP a;a.n=readc(f);a.d=readc(f);na[id]=a;}
 f>>k;vector<ON>sn;for(int j=0;j<k;j++)sn.push_back(readnode(f));f>>k;vector<ON>lin;for(int j=0;j<k;j++)lin.push_back(readnode(f));
 vector<long>outroots(98),outlabels(98),inputs(239);for(auto&i:outroots)f>>i;for(auto&i:outlabels)f>>i;for(auto&i:inputs)f>>i;
 vector<array<int,8>>jointmon(239);for(auto&m:jointmon)for(auto&i:m)f>>i;string extra;assert(!(f>>extra));
 vector<Mon>m3=oldmon,m4;sort(m3.begin(),m3.end());map<Mon,int>i3,i4;for(int i=0;i<84;i++)i3[m3[i]]=i;
 for(auto&m:jointmon)if(m[7]==0){Mon q;copy_n(m.begin(),7,q.begin());i4[q]=m4.size();m4.push_back(q);}assert(m4.size()==210);
 stage="cached_actual_K7_pure84_values";for(auto&v:oldpi)assert(knowncurveunit(R(v.d)));
 {ofstream o(folder+"/physical-cubic84.json");o<<"{\"kernel_genuine_index\":7,\"actual_nonpivot135\":68,\"values84\":[";for(int i=0;i<84;i++){if(i)o<<',';jr(o,oldpi[i]);}o<<"]}\n";}
 stage="accepted_cached_physical_to_joint_frame";for(auto&row:Bphys)for(auto&v:row)assert(knowncurveunit(R(v.d)));for(int k=0;k<7;k++)assert(eq(Bphys[3][k],k==6?rc("1"):R()));
 stage="one_joint_cubic84";vector<R>pi(84);
 for(int m=0;m<84;m++){map<Mon,R>image;Mon zm{};image.emplace(zm,rc("1"));for(int j=0;j<7;j++)for(int v=0;v<oldmon[m][j];v++){
  map<Mon,R>next;for(auto&p:image)for(int t=0;t<7;t++)if(!zero(Bphys[j][t])){Mon q=p.first;q[t]++;auto it=next.find(q);R term=mul(p.second,Bphys[j][t]);if(it==next.end())next.emplace(q,move(term));else it->second=add(it->second,term);}image=move(next);
 }for(auto&p:image)if(!zero(oldpi[m]))pi[i3.at(p.first)]=add(pi[i3.at(p.first)],mul(oldpi[m],p.second));}
 oldpi.clear();
 {ofstream o(folder+"/joint-cubic-frame.json");o<<"{\"joint_cubic_monomials84\":[";for(int i=0;i<84;i++){if(i)o<<',';o<<'[';for(int j=0;j<7;j++){if(j)o<<',';o<<m3[i][j];}o<<']';}o<<"],\"joint_cubic_values84\":[";for(int i=0;i<84;i++){if(i)o<<',';jr(o,pi[i]);}o<<"],\"physical_light_inverse7x7\":[";for(int i=0;i<7;i++){if(i)o<<',';o<<'[';for(int j=0;j<7;j++){if(j)o<<',';jr(o,Bphys[i][j]);}o<<']';}o<<"]}\n";}
 stage="NF117_scalar_closure_once";map<long,R>nav,sv;for(auto&p:na){R n=eval(p.second.n),d=eval(p.second.d);assert(knowncurveunit(d));nav.emplace(p.first,divunit(n,d,"accepted_NF117_atomic_named_unit_denominator"));}na.clear();
 for(auto&n:sn)sv.emplace(n.id,valnode(n,sv,nav));nav.clear();map<long,int>iv;for(int j=0;j<239;j++)iv[inputs[j]]=j;
 vector<int>nonzeros,units;long largestrow=0;
 for(int light=3;light<4;light++){
  stage="K7W_product_"+to_string(light);vector<R>psi(210);
  for(int m=0;m<84;m++)if(!zero(pi[m]))for(int j=0;j<7;j++)if(!zero(Bphys[light][j])){Mon q=m3[m];q[j]++;psi[i4.at(q)]=add(psi[i4.at(q)],mul(pi[m],Bphys[light][j]));}
  for(int m=0;m<210;m++){R check;for(int j=0;j<7;j++)if(m4[m][j]){Mon q=m4[m];q[j]--;check=add(check,mul(pi[i3.at(q)],Bphys[light][j]));}assert(eq(check,psi[m]));}
  map<long,int>uses;for(auto&n:lin){uses[n.id]+=0;for(long c:children(n))uses[c]++;}for(long r:outroots)uses[r]++;map<long,R>lv;
  for(auto&n:lin){R v;if(n.op=='v'){int c=iv.at(n.x);assert(jointmon[c][7]==0);Mon q;copy_n(jointmon[c].begin(),7,q.begin());v=psi[i4.at(q)];}else if(n.op=='z')v=R();else{assert(n.op=='l');for(auto&p:n.terms)v=add(v,mul(sv.at(p.first),lv.at(p.second)));}lv.emplace(n.id,move(v));for(long c:children(n))if(--uses[c]==0)lv.erase(c);}
  psi.clear();vector<R>out;for(long r:outroots)out.push_back(lv.at(r));lv.clear();int nz=0,un=0;long rowterms=0;
  for(auto&r:out){nz+=!zero(r);assert(knowncurveunit(R(r.d)));un+=knowncurveunit(r);rowterms+=rterms(r);}largestrow=max(largestrow,rowterms);nonzeros.push_back(nz);units.push_back(un);
  {ofstream o(folder+"/pure98-row"+to_string(light)+".json");o<<"{\"pure_kernel_genuine_index\":7,\"physical_light\":"<<vector<int>({0,1,2,3,4,6,7})[light]<<",\"target122_columns\":[";for(int j=0;j<98;j++){if(j)o<<',';o<<outlabels[j];}o<<"],\"values98\":[";for(int j=0;j<98;j++){if(j)o<<',';jr(o,out[j]);}o<<"],\"all210_combination_identities_checked\":true,\"all98_denominators_known_curve_units\":true,\"nonzero_count\":"<<nz<<",\"individual_known_unit_count\":"<<un<<"}\n";}
  // out, every linear value, and the210 input are released before the next row.
 }
 ofstream o(folder+"/audit.json");o<<"{\"all_passed\":true,\"cached_actual_K7_coefficients\":84,\"source_rows\":1,\"all210_joint_combination_identities_checked\":true,\"all98_NF_outputs_evaluated\":true,\"NF_scalar_closure_evaluated_once\":true,\"one_row_nonzero_counts\":[";
 for(int j=0;j<1;j++){if(j)o<<',';o<<nonzeros[j];}o<<"],\"one_row_individual_unit_counts\":[";for(int j=0;j<1;j++){if(j)o<<',';o<<units[j];}o<<"],\"max_row_RF_terms\":"<<largestrow<<",\"peak_live_RF_terms_with_inputreserve\":"<<peak_live+input_reserve<<",\"input_term_reserve\":"<<input_reserve<<",\"serialized_terms\":"<<stored<<",\"peak_polynomial_length\":"<<peak<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";
 }catch(const exception&e){ofstream o(folder+"/frontier.json");o<<"{\"all_passed\":false,\"stage\":\""<<stage<<"\",\"reason\":\""<<e.what()<<"\",\"peak_live_RF_terms_with_inputreserve\":"<<peak_live+input_reserve<<"}\n";return 2;}flint_cleanup();}
