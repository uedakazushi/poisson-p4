// Independently accepted peer arithmetic only; historical main is never called.
#define main previous_independent_pure98_main_NOT_CALLED
#include "s61_nonsquare_P100_pure98_review_20261104_v63.cpp"
#undef main
int main(){try{auto start=std::chrono::steady_clock::now();
 for(int k=0;k<5;k++){factors[k]=readP();powers[k]={constant("1")};}
 P z;fmpq_poly_set_coeff_si(z.v,1,1);z.sync();P zp20=pa(z,constant("20"));
 R a(pm(constant("-10/9"),zp20),{2,0,0,0,0}),e(pm(constant("1/9"),pa(z,constant("100"))),{1,0,0,0,0});
 R d(pm(constant("-3/100"),pm(z,zp20)));ap={rc("1"),a};ep={rc("1"),e};dp={rc("1"),d};
 std::vector<V>H(16,V(7)),sigma(7,V(16)),piv(6,V(84));V eta(7);
 for(auto&row:H)for(auto&v:row)v=eval(readS());monomials.clear();
 for(auto&row:sigma)for(auto&v:row)v=readR();for(auto&v:eta)v=readR();for(auto&row:piv)for(auto&v:row)v=readR();
 std::vector<int>Pids={1,6,8,9,11,14},ids={0,2,3,4,5,7,10,12,13,15};
 check(eq(eta[3],rc("1")),"eta physical coordinate3 unit1");
 for(int p=0;p<16;p++)check(zero(sigma[3][p]),"sigma row3 literallyzero");
 for(auto&row:H){R v;for(int j=0;j<7;j++)v=ra(v,rm(row[j],eta[j]));check(zero(v),"all16 H eta identities");}
 for(int j=0;j<7;j++)for(int k=0;k<7;k++){R v;for(int p=0;p<16;p++)v=ra(v,rm(sigma[j][p],H[p][k]));R want=j==k?rc("1"):R();if(k==3)want=rs(want,eta[j]);check(eq(v,want),"all49 source lift projection identities");}
 for(int p=0;p<16;p++)if(std::find(Pids.begin(),Pids.end(),p)==Pids.end())for(int j=0;j<7;j++)check(zero(sigma[j][p]),"sigma support only six pivots");
 auto alpha=[&](int i){V v(6);for(int p=0;p<6;p++)for(int j=0;j<7;j++)v[p]=ra(v[p],rm(H[i][j],sigma[j][Pids[p]]));return v;};
 for(int i=0;i<6;i++){auto v=alpha(Pids[i]);for(int j=0;j<6;j++)check(eq(v[j],rc(i==j?"1":"0")),"all36 pivot section identities");}
 R C=ra(ra(rm(rc("100"),a),rm(rc("18"),rm(e,e))),ra(rm(rc("-9"),e),rc("-2")));
 R expectedC(pm(constant("-25/9"),pa(z,constant("60"))),{1,0,0,0,0});check(eq(C,expectedC),"actual P100 C formula");
 R delta=rm(rc("-64000"),rm(C,C)),idelta=inv_unit(delta);
 long source_checks=0,Hchecks=0,identityblock=0;
 for(int i:ids){int actual,mode;std::cin>>actual>>mode;check(actual==i&&mode==int(i==0),"fixed nonpivot/source stream order");
  V ci(84);for(auto&v:ci)v=mode?readR():eval(readS());monomials.clear();
  auto al=alpha(i);V old(16),now(16);old[i]=rc("1");for(int p=0;p<6;p++)old[Pids[p]]=rn(al[p]);now=old;now[6]=rm(old[6],idelta);now[1]=ra(old[1],old[6]);
  for(int p=0;p<16;p++)check(eq(old[p],readR()),"all160 original source coefficients");
  for(int p=0;p<16;p++)check(eq(now[p],readR()),"all160 actualdelta source coefficients");
  for(int j=0;j<7;j++){R v,w;for(int p=0;p<16;p++){v=ra(v,rm(old[p],H[p][j]));R hp=p==6?rm(delta,rs(H[6][j],H[1][j])):H[p][j];w=ra(w,rm(now[p],hp));}check(zero(v)&&zero(w),"all140 source Hkernel identities");Hchecks+=2;}
  for(int k:ids){check(eq(old[k],rc(i==k?"1":"0")),"nonpivot I10 block");identityblock++;}
  for(int m=0;m<84;m++){
   R v=ci[m];for(int p=0;p<6;p++)v=rs(v,rm(al[p],piv[p][m]));
   R reverse=ci[m];for(int p=5;p>=0;p--)reverse=ra(reverse,rm(old[Pids[p]],piv[p][m]));
   R converted=ci[m];for(int p=0;p<6;p++){R cp=Pids[p]==6?rm(delta,rs(piv[1][m],piv[0][m])):piv[p][m];converted=ra(converted,rm(now[Pids[p]],cp));}
   check(eq(v,reverse)&&eq(v,converted)&&eq(v,readR()),"all1680 full pure source frame equalities and840 savedvalues");source_checks+=2;
  }
 }
 std::string extra;check(!(std::cin>>extra),"no trailing fixed source data");
 std::cout<<"{\"all_passed\":true,\"all16_Heta_49sigma_36pivot\":true,\"Hkernel_identity_checks\":"<<Hchecks<<",\"nonpivot_identity_block_checks\":"<<identityblock
 <<",\"pure_source_frame_identity_checks\":"<<source_checks<<",\"saved_values_checked\":840,\"whole_curve_named_denominators\":true,\"peak_polynomial_slots\":"<<PEAK
 <<",\"no_image_rank_or_quartic_job\":true,\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
