// One pinned cached C68 cubic times physical v3=joint W. No other source row.
#define main previous_independent_pure98_main_NOT_CALLED
#include "s61_nonsquare_P100_pure98_review_20261104_v63.cpp"
#undef main
int main(){try{
 auto start=std::chrono::steady_clock::now();
 for(int k=0;k<5;k++){factors[k]=readP();powers[k]={constant("1")};}
 P z;fmpq_poly_set_coeff_si(z.v,1,1);z.sync();P zp20=pa(z,constant("20"));
 R a(pm(constant("-10/9"),zp20),{2,0,0,0,0}),e(pm(constant("1/9"),pa(z,constant("100"))),{1,0,0,0,0});
 R d(pm(constant("-3/100"),pm(z,zp20)));ap={rc("1"),a};ep={rc("1"),e};dp={rc("1"),d};
 V parent(84);for(auto&v:parent)v=readR();auto oldmon=readMon(84);
 for(auto&v:parent)check(eq(v,readR()),"all84 saved physical values equal accepted literal C68 cache");
 std::vector<V>B(7,V(7));for(auto&row:B)for(auto&v:row)v=readR();
 for(auto&row:B)for(auto&v:row)check(eq(v,readR()),"all49 joint inverse entries match accepted frame");
 for(int j=0;j<7;j++)check(eq(B[3][j],rc(j==6?"1":"0")),"physical v3 is exactly joint W");
 auto mon3=oldmon;std::sort(mon3.begin(),mon3.end());check(mon3==readMon(84),"same joint84 monomial order");
 std::map<Mon,int>i3;for(int j=0;j<84;j++)i3[mon3[j]]=j;V joint(84);
 for(int j=0;j<84;j++){
  std::map<Mon,R>image;Mon zz{};image.emplace(zz,rc("1"));
  for(int k=0;k<7;k++)for(int c=0;c<oldmon[j][k];c++){
   std::map<Mon,R>next;
   for(auto&term:image)for(int l=0;l<7;l++)if(!zero(B[k][l])){
    Mon m=term.first;m[l]++;R val=rm(term.second,B[k][l]);auto it=next.find(m);
    if(it==next.end())next.emplace(m,std::move(val));else it->second=ra(it->second,val);
   }image=std::move(next);
  }
  for(auto&term:image)joint[i3.at(term.first)]=ra(joint[i3.at(term.first)],rm(parent[j],term.second));
 }
 parent.clear();for(auto&v:joint)check(eq(v,readR()),"all84 actual joint C68 coefficients");
 auto atom=readAtoms();auto sn=readNodes();check(sn.size()==244,"exact NF scalar244 closure");
 std::map<long,R>scalars;for(auto&n:sn)scalars.emplace(n.id,value(n,scalars,atom));atom.clear();sn.clear();
 auto ln=readNodes();check(ln.size()==547,"exact NF linear547 closure");
 std::vector<long>roots(98),labels(98),inputs(239);for(auto&i:roots)std::cin>>i;for(auto&i:labels)std::cin>>i;for(auto&i:inputs)std::cin>>i;
 std::vector<std::array<int,8>>allmon(239);for(auto&m:allmon)for(int&i:m)std::cin>>i;
 std::vector<Mon>mon4;std::map<Mon,int>i4;
 for(auto&m:allmon)if(!m[7]){Mon q;std::copy_n(m.begin(),7,q.begin());i4[q]=mon4.size();mon4.push_back(q);}
 check(mon4.size()==210,"complete pure210 target universe");
 V psi(210);
 for(int j=0;j<84;j++)for(int l=0;l<7;l++)if(!zero(B[3][l])){Mon m=mon3[j];m[l]++;psi[i4.at(m)]=ra(psi[i4.at(m)],rm(joint[j],B[3][l]));}
 for(int j=0;j<210;j++){
  R reverse;for(int l=0;l<7;l++)if(mon4[j][l]){Mon m=mon4[j];m[l]--;reverse=ra(reverse,rm(joint[i3.at(m)],B[3][l]));}
  check(eq(reverse,psi[j]),"all210 full product coefficient identities");
 }
 std::map<long,int>inputcol;for(int j=0;j<239;j++)inputcol[inputs[j]]=j;
 std::map<long,long>uses;for(auto&n:ln){uses[n.id]+=0;for(long c:children(n))uses[c]++;}for(long i:roots)uses[i]++;
 std::map<long,R>v;
 for(auto&n:ln){
  R q;if(n.op=='v'){auto&m=allmon[inputcol.at(n.a)];check(!m[7],"NF pure closure has no H input");Mon a;std::copy_n(m.begin(),7,a.begin());q=psi[i4.at(a)];}
  else q=value(n,v,scalars);
  v.emplace(n.id,std::move(q));for(long c:children(n))if(--uses[c]==0)v.erase(c);
 }
 int nz=0,units=0;
 for(int j=0;j<98;j++){R q=v.at(roots[j]);check(eq(q,readR()),"all98 literal C68*W NF outputs");nz+=!zero(q);units+=named_unit(q.n);if(--uses[roots[j]]==0)v.erase(roots[j]);}
 check(nz==40&&units==0,"row support metadata40 and no individual named unit");
 std::string extra;check(!(std::cin>>extra),"no trailing source row");
 std::cout<<"{\"all_passed\":true,\"source_rows\":1,\"accepted_parent84_physical_and_joint_verified\":true,\"all49_frame_entries_verified\":true,\"physical_v3_equals_joint_W\":true,\"all210_product_coefficients_checked\":true,\"all98_NF_outputs_checked\":true,\"nonzero_count\":"<<nz<<",\"individual_named_unit_count\":"<<units<<",\"peak_polynomial_slots\":"<<PEAK<<",\"no_other_row_source_rank_or_factor_query\":true,\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
 return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
