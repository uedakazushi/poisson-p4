// Reuse only the independently accepted peer's QQ[z] arithmetic/node reader.
#define main previous_independent_pure_row_main_NOT_CALLED
#include "s61_nonsquare_P100_pure98_review_20261104_v63.cpp"
#undef main
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_mat.h>
int main(){try{
 auto start=std::chrono::steady_clock::now();
 for(int k=0;k<5;k++){factors[k]=readP();powers[k]={constant("1")};}
 P z;fmpq_poly_set_coeff_si(z.v,1,1);z.sync();P zp20=pa(z,constant("20"));
 std::vector<P> rawf={z,zp20,pa(z,constant("60")),pa(pm(constant("203"),z),constant("2700")),pa(z,constant("12"))};
 std::vector<int> selectedRows(8),selectedCols(8);for(auto&i:selectedRows)std::cin>>i;for(auto&i:selectedCols)std::cin>>i;
 std::vector<std::array<int,5>> ex(8);std::vector<P> lcd(8);
 for(int i=0;i<8;i++){lcd[i]=constant("1");for(int k=0;k<5;k++){std::cin>>ex[i][k];for(int j=0;j<ex[i][k];j++)lcd[i]=pm(lcd[i],rawf[k]);}}
 R aa(pm(constant("-10/9"),zp20),{2,0,0,0,0}),ee(pm(constant("1/9"),pa(z,constant("100"))),{1,0,0,0,0});
 R dd(pm(constant("-3/100"),pm(z,zp20)));ap={rc("1"),aa};ep={rc("1"),ee};dp={rc("1"),dd};
 std::vector<R> oldpi(84);for(auto&v:oldpi)v=eval(readS());auto oldmon=readMon(84);
 for(auto&v:oldpi)check(eq(v,readR()),"all84 full saved G18 values");monomials.clear();
 R CC=ra(ra(rm(rc("100"),aa),rm(rc("18"),rm(ee,ee))),ra(rm(rc("-9"),ee),rc("-2")));
 R jac=rm(rc("3200/81"),rm(rm(aa,CC),rm(dd,dd))),ijac=inv_unit(jac);
 std::vector<V>NN(7,V(7)),backward(7,V(7)),ff(7,V(7)),forward(7,V(7));
 for(auto&row:NN)for(auto&v:row)v=rm(eval(readS()),ijac);
 int order[7]={0,1,2,6,3,5,4};
 for(int i=0;i<7;i++)for(int j=0;j<7;j++)backward[i][j]=i==5?rs(ra(NN[5][j],NN[1][j]),rm(rc("3"),NN[0][j])):NN[order[i]][j];NN.clear();
 for(auto&row:ff)for(auto&v:row)v=eval(readS());
 int Hi[7][7];for(auto&row:Hi)for(int&v:row)std::cin>>v;
 for(int i=0;i<7;i++)for(int j=0;j<7;j++)for(int k=0;k<7;k++)forward[i][j]=ra(forward[i][j],rm(ff[i][k],rc(std::to_string(Hi[j][k]))));ff.clear();
 for(int i=0;i<7;i++)for(int j=0;j<7;j++){
  R x,y;for(int k=0;k<7;k++){x=ra(x,rm(forward[i][k],backward[k][j]));y=ra(y,rm(backward[i][k],forward[k][j]));}
  check(eq(x,rc(i==j?"1":"0"))&&eq(y,rc(i==j?"1":"0")),"both coordinate inverse products");check(eq(backward[i][j],readR()),"all49 literal inverse values");
 }forward.clear();monomials.clear();
 auto mon3=oldmon;std::sort(mon3.begin(),mon3.end());std::map<Mon,int>i3;
 for(int i=0;i<84;i++)i3[mon3[i]]=i;auto savedmon3=readMon(84);check(savedmon3==mon3,"joint monomial84 ordering");
 V joint(84);
 for(int m=0;m<84;m++){
  std::map<Mon,R> image;Mon zz{};image.emplace(zz,rc("1"));
  for(int j=0;j<7;j++)for(int n=0;n<oldmon[m][j];n++){
   std::map<Mon,R> next;for(auto&pair:image)for(int k=0;k<7;k++)if(!zero(backward[j][k])){
    Mon q=pair.first;q[k]++;R v=rm(pair.second,backward[j][k]);auto it=next.find(q);if(it==next.end())next.emplace(q,std::move(v));else it->second=ra(it->second,v);
   }image=std::move(next);
  }
  for(auto&pair:image)joint[i3.at(pair.first)]=ra(joint[i3.at(pair.first)],rm(oldpi[m],pair.second));
 }
 oldpi.clear();for(auto&v:joint)check(eq(v,readR()),"all84 joint G18 values");
 auto nfAtoms=readAtoms();auto scalarNodes=readNodes();check(scalarNodes.size()==244,"fixed scalar244 closure");std::map<long,R>scalar;
 for(auto&n:scalarNodes)scalar.emplace(n.id,value(n,scalar,nfAtoms));nfAtoms.clear();scalarNodes.clear();
 auto linearNodes=readNodes();check(linearNodes.size()==547,"fixed linear547 closure");std::vector<long>roots(98),labels(98),inputs(239);
 for(auto&i:roots)std::cin>>i;for(auto&i:labels)std::cin>>i;for(auto&i:inputs)std::cin>>i;
 std::vector<std::array<int,8>>allmon(239);for(auto&m:allmon)for(int&i:m)std::cin>>i;
 std::map<Mon,int>i4;std::vector<Mon>mon4;for(auto&m:allmon)if(!m[7]){Mon q;std::copy_n(m.begin(),7,q.begin());i4[q]=mon4.size();mon4.push_back(q);}check(mon4.size()==210,"pure210 inputs");
 std::map<long,int>column;for(int i=0;i<239;i++)column[inputs[i]]=i;
 std::vector<std::vector<P>>selectedBySource(8,std::vector<P>(8));long clearings=0;
 auto clearrow=[&](const V&row,int source){for(int j=0;j<98;j++){
   R cleared=rm(row[j],R(lcd[source]));for(int k=0;k<5;k++)check(!cleared.e[k],"rowLCD wholecurve unit clears denominator");clearings++;
   for(int c=0;c<8;c++)if(selectedCols[c]==j)selectedBySource[source][c]=cleared.n;
  }};
 V q5(98);for(auto&v:q5)v=readR();clearrow(q5,0);q5.clear();
 for(int light=0;light<7;light++){
  V psi(210);
  for(int c=0;c<84;c++)for(int j=0;j<7;j++)if(!zero(joint[c])&&!zero(backward[light][j])){
   Mon m=mon3[c];m[j]++;psi[i4.at(m)]=ra(psi[i4.at(m)],rm(joint[c],backward[light][j]));
  }
  for(int m=0;m<210;m++){R v;for(int j=0;j<7;j++)if(mon4[m][j]){Mon n=mon4[m];n[j]--;v=ra(v,rm(joint[i3.at(n)],backward[light][j]));}check(eq(v,psi[m]),"all1470 product coefficient identities");}
  std::map<long,long>uses;for(auto&n:linearNodes){uses[n.id]+=0;for(long c:children(n))uses[c]++;}for(long i:roots)uses[i]++;
  std::map<long,R>linear;
  for(auto&n:linearNodes){R v;if(n.op=='v'){auto&m=allmon[column.at(n.a)];check(!m[7],"NF pure dependencies have H degree0");Mon q;std::copy_n(m.begin(),7,q.begin());v=psi[i4.at(q)];}else v=value(n,linear,scalar);linear.emplace(n.id,std::move(v));for(long c:children(n))if(--uses[c]==0)linear.erase(c);}
  V out(98);for(int i=0;i<98;i++){out[i]=linear.at(roots[i]);check(eq(out[i],readR()),"all686 actual NF outputs");if(--uses[roots[i]]==0)linear.erase(roots[i]);}clearrow(out,light+1);
 }
 backward.clear();joint.clear();scalar.clear();linearNodes.clear();check(clearings==784,"all784 rowLCD identities");
 std::vector<std::vector<P>>S(8,std::vector<P>(8)),adj(8,std::vector<P>(8));
 for(int i=0;i<8;i++)for(int j=0;j<8;j++){S[i][j]=readP();check(fmpq_poly_equal(S[i][j].v,selectedBySource[selectedRows[i]][j].v),"all64 actual selected source coefficients");}selectedBySource.clear();
 P delta=readP();check(!fmpq_poly_is_zero(delta.v),"nonzero200degree polynomial");
 for(auto&row:adj)for(auto&v:row)v=readP();
 // A different fixed determinant algorithm, not the producer's bordered search.
 fmpz_poly_mat_t Z;fmpz_poly_mat_init(Z,8,8);fmpz_t den,di,dj;fmpz_init(den);fmpz_init(di);fmpz_init(dj);fmpz_one(den);fmpq_poly_t scaled;fmpq_poly_init(scaled);
 for(int i=0;i<8;i++){fmpz_one(di);for(int j=0;j<8;j++){fmpq_poly_get_denominator(dj,S[i][j].v);fmpz_lcm(di,di,dj);}fmpz_mul(den,den,di);
  for(int j=0;j<8;j++){fmpq_poly_scalar_mul_fmpz(scaled,S[i][j].v,di);fmpq_poly_get_denominator(dj,scaled);check(fmpz_is_one(dj),"integer row clearing");fmpq_poly_get_numerator(fmpz_poly_mat_entry(Z,i,j),scaled);}
 }
 fmpz_poly_t det;fmpz_poly_init(det);fmpz_poly_mat_det_fflu(det,Z);P own;fmpq_poly_set_fmpz_poly(own.v,det);fmpq_poly_scalar_div_fmpz(own.v,own.v,den);own.sync();check(fmpq_poly_equal(own.v,delta.v),"independent fixed8 determinant equals saved");
 for(int i=0;i<8;i++)for(int j=0;j<8;j++){P l,r;for(int k=0;k<8;k++){l=pa(l,pm(S[i][k],adj[k][j]));r=pa(r,pm(adj[i][k],S[k][j]));}P expected=i==j?delta:P();check(fmpq_poly_equal(l.v,expected.v)&&fmpq_poly_equal(r.v,expected.v),"both128 adjugate polynomial products");}
 std::string extra;check(!(std::cin>>extra),"no trailing source input");
 std::cout<<"{\"all_passed\":true,\"source84\":true,\"joint84_inverse49_bothproducts\":true,\"all1470_products_and686_NF_outputs\":true,\"all784_LCD_and64_selected_source_coefficients\":true,\"independent_fixed8_det_equal\":true,\"both128_adjugate_products\":true,\"det_degree\":"<<fmpq_poly_degree(delta.v)<<",\"peak_polynomial_slots\":"<<PEAK<<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
 fmpz_poly_clear(det);fmpz_poly_mat_clear(Z);fmpz_clear(den);fmpz_clear(di);fmpz_clear(dj);fmpq_poly_clear(scaled);return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
