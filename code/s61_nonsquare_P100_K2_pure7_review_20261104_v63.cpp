// Generic cached-parent84 -> seven actual products/NF98 independent checker.
#define main previous_independent_pure98_main_NOT_CALLED
#include "s61_nonsquare_P100_pure98_review_20261104_v63.cpp"
#undef main
int main(){try{auto start=std::chrono::steady_clock::now();
 for(int k=0;k<5;k++){factors[k]=readP();powers[k]={constant("1")};}
 P z;fmpq_poly_set_coeff_si(z.v,1,1);z.sync();P zp20=pa(z,constant("20"));
 R a(pm(constant("-10/9"),zp20),{2,0,0,0,0}),e(pm(constant("1/9"),pa(z,constant("100"))),{1,0,0,0,0});
 R d(pm(constant("-3/100"),pm(z,zp20)));ap={rc("1"),a};ep={rc("1"),e};dp={rc("1"),d};
 V parent84(84);for(auto&v:parent84)v=readR();auto oldmon=readMon(84);for(auto&v:parent84)check(eq(v,readR()),"all84 physical values equal accepted parent");
 std::vector<V>back(7,V(7));for(auto&row:back)for(auto&v:row)v=readR();for(auto&row:back)for(auto&v:row)check(eq(v,readR()),"same49 accepted wholecurve inverse values");
 auto mon3=oldmon;std::sort(mon3.begin(),mon3.end());std::map<Mon,int>i3;for(int i=0;i<84;i++)i3[mon3[i]]=i;
 check(mon3==readMon(84),"joint84 monomial order");V joint(84);
 for(int m=0;m<84;m++){
  std::map<Mon,R>image;Mon zz{};image.emplace(zz,rc("1"));
  for(int j=0;j<7;j++)for(int c=0;c<oldmon[m][j];c++){
   std::map<Mon,R>next;for(auto&v:image)for(int k=0;k<7;k++)if(!zero(back[j][k])){
    Mon q=v.first;q[k]++;R term=rm(v.second,back[j][k]);auto it=next.find(q);if(it==next.end())next.emplace(q,std::move(term));else it->second=ra(it->second,term);
   }image=std::move(next);
  }
  for(auto&v:image)joint[i3.at(v.first)]=ra(joint[i3.at(v.first)],rm(parent84[m],v.second));
 }
 parent84.clear();for(auto&v:joint)check(eq(v,readR()),"all84 actual joint parent values");
 auto atoms=readAtoms();auto sn=readNodes();check(sn.size()==244,"same scalar244 closure");std::map<long,R>scalar;
 for(auto&n:sn)scalar.emplace(n.id,value(n,scalar,atoms));atoms.clear();sn.clear();
 auto ln=readNodes();check(ln.size()==547,"same linear547 closure");std::vector<long>roots(98),labels(98),inputs(239);
 for(auto&i:roots)std::cin>>i;for(auto&i:labels)std::cin>>i;for(auto&i:inputs)std::cin>>i;
 std::vector<std::array<int,8>>allmon(239);for(auto&m:allmon)for(int&i:m)std::cin>>i;
 std::vector<Mon>mon4;std::map<Mon,int>i4;for(auto&m:allmon)if(!m[7]){Mon q;std::copy_n(m.begin(),7,q.begin());i4[q]=mon4.size();mon4.push_back(q);}check(mon4.size()==210,"pure210 targets");
 std::map<long,int>inputcol;for(int i=0;i<239;i++)inputcol[inputs[i]]=i;
 std::vector<int>nzs;
 for(int light=0;light<7;light++){
  V psi(210);
  for(int c=0;c<84;c++)for(int j=0;j<7;j++)if(!zero(joint[c])&&!zero(back[light][j])){Mon m=mon3[c];m[j]++;psi[i4.at(m)]=ra(psi[i4.at(m)],rm(joint[c],back[light][j]));}
  for(int m=0;m<210;m++){R v;for(int j=0;j<7;j++)if(mon4[m][j]){Mon q=mon4[m];q[j]--;v=ra(v,rm(joint[i3.at(q)],back[light][j]));}check(eq(v,psi[m]),"all1470 full product identities");}
  std::map<long,long>uses;for(auto&n:ln){uses[n.id]+=0;for(long c:children(n))uses[c]++;}for(long i:roots)uses[i]++;
  std::map<long,R>values;
  for(auto&n:ln){R v;if(n.op=='v'){auto&m=allmon[inputcol.at(n.a)];check(!m[7],"NF inputs preserve pure Hdegree0");Mon q;std::copy_n(m.begin(),7,q.begin());v=psi[i4.at(q)];}else v=value(n,values,scalar);values.emplace(n.id,std::move(v));for(long c:children(n))if(--uses[c]==0)values.erase(c);}
  int nz=0;for(int j=0;j<98;j++){R v=values.at(roots[j]);check(eq(v,readR()),"all686 exact pure98 source outputs");nz+=!zero(v);if(--uses[roots[j]]==0)values.erase(roots[j]);}nzs.push_back(nz);
 }
 std::string extra;check(!(std::cin>>extra),"no trailing source values");
 std::cout<<"{\"all_passed\":true,\"cached_parent84_physical_and_joint_verified\":true,\"same49_accepted_inverse\":true,\"all1470_actual_products_and686_NF_values\":true,\"nonzero_counts\":[";
 for(int i=0;i<7;i++){if(i)std::cout<<',';std::cout<<nzs[i];}
 std::cout<<"],\"peak_polynomial_slots\":"<<PEAK<<",\"no_rank_or_new_source_query\":true,\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
