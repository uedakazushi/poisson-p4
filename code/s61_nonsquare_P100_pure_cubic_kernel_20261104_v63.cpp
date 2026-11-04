// Stream existing core polynomials; no oracle, quartic NF or rank job.
#define main frozen_q5_row_main_NOT_CALLED
#include "s61_nonsquare_P100_q5_pure98_row_20261104_v63.cpp"
#undef main
R nextcore(istream&f){C c=readc(f);return eval(c);}
bool curveunit(const R&r){if(zero(r))return false;P p=r.n;for(int i=0;i<100;i++){P g;fmpq_poly_gcd(g.v,p.v,unitproduct.v);guard(g);if(fmpq_poly_degree(g.v)==0)break;p=quotient(p,g);}return fmpq_poly_degree(p.v)==0;}
void emitvec(ostream&o,const vector<R>&v){o<<'[';for(int j=0;j<int(v.size());j++){if(j)o<<',';jr(o,v[j]);}o<<']';}
int main(int argc,char**argv){assert(argc==3);began=chrono::steady_clock::now();string folder=argv[2];try{
 ifstream f(argv[1]);f>>input_reserve;z=shiftp(one(),1);P zp20=plusp(z,pc("20"));ap={one(),scalar(zp20,"-10/9")};ep={one(),scalar(plusp(z,pc("100")),"1/9")};dp={one(),scalar(mulp(z,zp20),"-3/100")};
 unitproduct=mulp(mulp(mulp(mulp(z,zp20),plusp(z,pc("60"))),plusp(scalar(z,"203"),pc("2700"))),plusp(z,pc("12")));
 vector<int>Pids={1,6,8,9,11,14};vector<vector<R>>H(16,vector<R>(7)),sig(7,vector<R>(16)),piv(6,vector<R>(84));vector<R>eta(7);
 stage="raw_H_and_cached_sigma";for(auto&row:H)for(auto&v:row)v=nextcore(f);for(auto&row:sig)for(auto&v:row)v=inputr(f);for(auto&v:eta)v=inputr(f);for(auto&row:piv)for(auto&v:row)v=inputr(f);
 assert(eq(eta[3],rc("1")));for(auto&row:H){R s;for(int j=0;j<7;j++)s=add(s,mul(row[j],eta[j]));assert(zero(s));}
 for(int j=0;j<7;j++)for(int k=0;k<7;k++){R s;for(int p=0;p<16;p++)s=add(s,mul(sig[j][p],H[p][k]));R want=j==k?rc("1"):R();if(k==3)want=sub(want,eta[j]);assert(eq(s,want));}
 for(int p=0;p<16;p++)if(find(Pids.begin(),Pids.end(),p)==Pids.end())for(int j=0;j<7;j++)assert(zero(sig[j][p]));
 auto alpha=[&](int i){vector<R>a(6);for(int p=0;p<6;p++)for(int j=0;j<7;j++)a[p]=add(a[p],mul(H[i][j],sig[j][Pids[p]]));return a;};
 for(int r=0;r<6;r++){auto a=alpha(Pids[r]);for(int p=0;p<6;p++)assert(eq(a[p],rc(p==r?"1":"0")));}
 R Cbar(pc("-25/9"));Cbar=divunit(mul(Cbar,R(plusp(z,pc("60")))),R(z),"known_z_unit");R delta=mul(rc("-64000"),powr(Cbar,2));assert(curveunit(delta));
 long count;f>>count;assert(count==10);vector<int>ids;vector<long>rowterms;vector<int>nzs;long pure_checks=0;
 for(int t=0;t<10;t++){
  int i,mode;f>>i>>mode;ids.push_back(i);assert(find(Pids.begin(),Pids.end(),i)==Pids.end());stage="kernel_nonpivot_"+to_string(i);vector<R>pi(84);
  for(auto&v:pi)v=mode==1?inputr(f):nextcore(f);
  auto a=alpha(i);vector<R>old(16),now(16);old[i]=rc("1");for(int p=0;p<6;p++)old[Pids[p]]=neg(a[p]);now=old;now[6]=divunit(old[6],delta,"actual_delta7_source_unit");now[1]=add(old[1],old[6]);
  for(int j=0;j<7;j++){R l,r;for(int s=0;s<16;s++){l=add(l,mul(old[s],H[s][j]));R h=s==6?mul(delta,sub(H[6][j],H[1][j])):H[s][j];r=add(r,mul(now[s],h));}assert(zero(l)&&zero(r));}
  for(int j=0;j<16;j++)if(find(Pids.begin(),Pids.end(),j)==Pids.end())assert(eq(old[j],rc(j==i?"1":"0")));
  vector<R>out(84);int nz=0;long terms=0;
  for(int m=0;m<84;m++){
   R v=pi[m];for(int p=0;p<6;p++)v=sub(v,mul(a[p],piv[p][m]));R reverse=pi[m];for(int p=5;p>=0;p--)reverse=add(reverse,mul(old[Pids[p]],piv[p][m]));assert(eq(v,reverse));pure_checks++;
   R newer=pi[m];for(int p=0;p<6;p++){R c=Pids[p]==6?mul(delta,sub(piv[1][m],piv[0][m])):piv[p][m];newer=add(newer,mul(now[Pids[p]],c));}assert(eq(v,newer));pure_checks++;
   assert(curveunit(R(v.d)));out[m]=v;nz+=!zero(v);terms+=rterms(v);
  }
  for(auto&v:old)assert(curveunit(R(v.d)));for(auto&v:now)assert(curveunit(R(v.d)));rowterms.push_back(terms);nzs.push_back(nz);
  {ofstream o(folder+"/pure-cubic84-"+to_string(i)+".json");o<<"{\"nonpivot_genuine_index\":"<<i<<",\"old16_source_cofactors\":";emitvec(o,old);o<<",\"actual_delta7_source_cofactors\":";emitvec(o,now);o<<",\"pure84_values\":";emitvec(o,out);o<<",\"both_full_source_frame_identities_verified\":true,\"all7_Hlinear_coefficients_zero\":true,\"whole_curve_denominators_only\":true}\n";}
 }
 string extra;assert(!(f>>extra));ofstream o(folder+"/audit.json");o<<"{\"all_passed\":true,\"all16_rawH_eta_identities\":true,\"all49_sigma_lift_identities\":true,\"all36_pivot_projection_identities\":true,\"all140_kernel_H_identities_both_frames\":true,\"all100_nonpivot_identity_block_entries\":true,\"pure_coefficient_identity_checks\":"<<pure_checks<<",\"nonpivot_indices\":[";
 for(int i=0;i<10;i++){if(i)o<<',';o<<ids[i];}o<<"],\"pure84_nonzero_counts\":[";for(int i=0;i<10;i++){if(i)o<<',';o<<nzs[i];}o<<"],\"ten_RF_termcounts\":[";for(int i=0;i<10;i++){if(i)o<<',';o<<rowterms[i];}o<<"],\"peak_live_RF_terms_with_stream_input\":"<<peak_live+input_reserve<<",\"input_stream_reserve\":"<<input_reserve<<",\"serialized_terms\":"<<stored<<",\"native_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-began).count()<<"}\n";
 }catch(const exception&e){ofstream o(folder+"/frontier.json");o<<"{\"all_passed\":false,\"stage\":\""<<stage<<"\",\"reason\":\""<<e.what()<<"\"}\n";return 2;}flint_cleanup();}
