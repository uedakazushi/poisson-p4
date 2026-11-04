// Only denominator-unit and one scalar identity checks; no large DAG expansion.
#define main independent_Lframe_inverse_audit_main
#include "s61_nonsquare_ell_Lframe_review_20261104_v63.cpp"
#undef main
struct Ratio{Poly n,d;Ratio(){d=constant(1);}};
int main(int argc,char**argv){assert(argc==2);fmpq_mpoly_ctx_init(qc,3,ORD_LEX);{
 auto started=chrono::steady_clock::now();ifstream f(argv[1]);int na,nn,root,theta;f>>na>>nn>>root>>theta;assert(na==1012&&nn==18697);
 Poly a,e,d;fmpq_mpoly_gen(a.x,0,qc);fmpq_mpoly_gen(e.x,1,qc);fmpq_mpoly_gen(d.x,2,qc);
 Poly F=scale(a,500);add(F,scale(multiply(e,e),90));add(F,scale(e,3));add(F,constant(-24));
 Poly B=scale(multiply(a,d),50);add(B,scale(multiply(multiply(e,e),d),18));add(B,scale(multiply(e,d),-9));add(B,scale(d,-2));add(B,scale(e,60));add(B,constant(-40));
 Poly C=scale(a,100);add(C,scale(multiply(e,e),18));add(C,scale(e,-9));add(C,constant(-2));
 Poly J=scale(multiply(a,d),5);add(J,scale(e,-6));add(J,constant(-4));
 vector<Ratio>atoms(na);long nonconstant=0;
 for(auto&v:atoms){load(f,v.n);fmpq_mpoly_zero(v.d.x,qc);load(f,v.d);assert(!zero(v.d));Poly rem=v.d;nonconstant+=fmpq_mpoly_total_degree_si(rem.x,qc)>0;
  for(auto&factor:vector<Poly>{a,d,F,B,C,J})while(fmpq_mpoly_total_degree_si(rem.x,qc)>0){Poly quotient;if(!fmpq_mpoly_divides(quotient.x,rem.x,factor.x,qc))break;rem=quotient;}
  assert(!zero(rem)&&fmpq_mpoly_is_fmpq(rem.x,qc));
 }
 vector<array<int,3>>nodes(nn);for(int i=0;i<nn;i++){
  int op,x,y;f>>op>>x>>y;assert(f&&op>=0&&op<=4);nodes[i]={op,x,y};
  if(op==0)assert(i==0);else if(op==1)assert(x>=0&&x<na);else assert(x>=0&&x<i&&y>=0&&y<i);
 }
 string tail;assert(!(f>>tail));map<int,Ratio>done;
 auto value=[&](auto&&self,int i)->Ratio{
  if(done.count(i))return done.at(i);auto[op,x,y]=nodes[i];Ratio r;
  if(op==0){}else if(op==1)r=atoms[x];else{
   Ratio l=self(self,x),v=self(self,y);r.d=multiply(l.d,v.d);
   if(op==4)r.n=multiply(l.n,v.n);else{r.n=multiply(l.n,v.d);add(r.n,scale(multiply(v.n,l.d),op==2?1:-1));}
  }assert(!zero(r.d));done[i]=r;return r;
 };
 Ratio z=value(value,root);assert(zero(z.n));Ratio th=value(value,theta);
 Poly expected=multiply(multiply(a,multiply(C,multiply(C,C))),J);fmpq_t q;fmpq_init(q);fmpq_set_si(q,-1024000,243);fmpq_mpoly_scalar_mul_fmpq(expected.x,expected.x,q,qc);fmpq_clear(q);
 assert(same(th.n,multiply(expected,th.d)));
 cout<<"{\"all_passed\":true,\"all1012_atom_denominators_only_FULLU_units\":true,\"nonconstant_atom_denominators\":"<<nonconstant<<",\"all18697_node_references_topological\":true,\"one_scalar_separator_identity_exact_QQ\":true,\"theta_is_literal_known_unit\":true,\"evaluated_nodes\":"<<done.size()<<",\"no_other_scalar_roots_expanded\":true,\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-started).count()<<"}"<<endl;
 }fmpq_mpoly_ctx_clear(qc);flint_cleanup();return 0;}
