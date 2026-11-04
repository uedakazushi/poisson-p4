#!/usr/bin/env python3
"""Independent oldpoint finite-free10 arithmetic and balanced primitive-element audit."""
import argparse,gzip,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_ZN_PRIMITIVE_KRYLOV_20261104_v63_REVIEW';P=32003
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def mod(x):
 q=Q(str(x));assert q.denominator%P;return q.numerator*pow(q.denominator,-1,P)%P
def polyadd(a,b):
 d=dict(a)
 for e,c in b.items():d[e]=d.get(e,Q(0))+c
 return {e:c for e,c in d.items()if c}
def polyscale(a,c):return {e:v*c for e,v in a.items()if v*c}
def polymul(a,b):
 d={}
 for e,c in a.items():
  for f,v in b.items():
   g=tuple(x+y for x,y in zip(e,f));d[g]=d.get(g,Q(0))+c*v
 return {e:c for e,c in d.items()if c}
def main(out):
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
 folder=ROOT/'reports/Q_NEW_FAMILY_ZN_PRIMITIVE_KRYLOV_20261104_v63';ar=read(folder/'report.json');ac=read(folder/'certificate.json')
 assert ar['all_passed']and ar['ZN_primitive_degree10_proved']and ar['old_saved_beta_point_index']==0 and ar['prime']==P
 for rel,h in ar['input_sha256'].items():assert sha(ROOT/rel)==h,rel
 for rel,h in ar['evidence_sha256'].items():assert sha(folder/rel)==h,rel
 pointp=ROOT/'reports/Q_NEW_FAMILY_OCTIC_BASE_DOMINANCE_20261103_v63/certificate.json.gz';point=read(pointp);assert point['same_existing_point_index']==0
 spinp=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz';spin=read(spinp)
 beta=list(map(mod,point['beta16']));assert beta==ac['beta16_modp']and point['beta16']==ac['beta16_QQ']
 mon=spin['quadratic_monomials'];quadratic=[beta[i]*beta[j]%P for i,j in mon]
 mats=[]
 for name in ['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136','moment_transpose_coefficients16x136']:
  vals=[sum(mod(a)*b for a,b in zip(row,quadratic))%P for row in spin[name]];mats.append([vals[4*i:4*i+4]for i in range(4)])
 C,D,M=mats;assert C==ac['C4_modp']and D==ac['D4_modp']
 assert all((M[i][j]+C[i][j]-2*D[i][j])%P==0 for i in range(4)for j in range(4))
 def mm(a,b):return [[sum(a[i][k]*b[k][j]for k in range(4))%P for j in range(4)]for i in range(4)]
 f4=sum(mm(M,M)[i][i]for i in range(4))%P;assert f4==ac['F4_modp']and f4
 # Independent Newton determinant of C+tD, unlike the author's24 permutations.
 def pmatrix(a,b):
  c=[[[0]*4 for _ in range(4)]for _ in range(len(a)+len(b)-1)]
  for k,A in enumerate(a):
   for l,B in enumerate(b):
    AB=mm(A,B)
    for i in range(4):
     for j in range(4):c[k+l][i][j]=(c[k+l][i][j]+AB[i][j])%P
  return c
 assert sum(C[i][i]for i in range(4))%P==sum(D[i][i]for i in range(4))%P==0
 second=pmatrix([C,D],[C,D]);fourth=pmatrix(second,second)
 tr2=[sum(a[i][i]for i in range(4))%P for a in second]
 tr4=[sum(a[i][i]for i in range(4))%P for a in fourth]
 tr2sq=[sum(tr2[i]*tr2[k-i]for i in range(3)if 0<=k-i<3)%P for k in range(5)]
 pencil=[(a*pow(8,-1,P)-b*pow(4,-1,P))%P for a,b in zip(tr2sq,tr4)]
 assert pencil==ac['pencil_p0_to_p4_modp']and pencil[0]
 assert(f4*f4-400*pencil[0]-1600*pencil[4])%P==0
 invf=pow(f4*f4%P,-1,P)
 base={'R':(pencil[3]+4*pencil[4])*pow(pencil[0],-1,P)%P,'Z0':pencil[0]*invf%P,'Z1':pencil[1]*invf%P,'Z2':pencil[2]*invf%P}
 assert base==ac['coefficient_base_modp']
 gp=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_CONTENT_GCD_20261104_v63/primitive.txt';g=[0]*11
 with gp.open()as f:
  assert f.readline().split()==['6','891']
  for line in f:
   c,zero,r,R,z0,z1,z2=map(int,line.split());assert zero==0
   g[r]=(g[r]+c*pow(base['R'],R,P)*pow(base['Z0'],z0,P)*pow(base['Z1'],z1,P)*pow(base['Z2'],z2,P))%P
 assert g==ac['specialized_G_coefficients_ascending']and g[-1]
 gg=[c*pow(g[-1],-1,P)%P for c in g];assert gg==ac['monic_G_coefficients_ascending']
 def reduce(a):
  b=[c%P for c in a]
  while len(b)>10:
   k=len(b)-11;c=b[-1]
   for i in range(11):b[k+i]=(b[k+i]-c*gg[i])%P
   assert b[-1]==0;b.pop()
  return b+[0]*(10-len(b))
 zero=[0]*10;one=[1]+[0]*9;rr=[0,1]+[0]*8
 def add(a,b):return[(x+y)%P for x,y in zip(a,b)]
 def scale(a,c):return[(x*c)%P for x in a]
 def mul(a,b):
  c=[0]*19
  for i,x in enumerate(a):
   if x:
    for j,y in enumerate(b):c[i+j]=(c[i+j]+x*y)%P
  return reduce(c)
 def power(a,n):
  z=one
  for _ in range(n):z=mul(z,a)
  return z
 cp=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_PRS_CONTENT_PREFLIGHT_20261104_v63/C1-primitive.txt.gz';c1=[]
 with gzip.open(cp,'rt')as f:
  assert f.readline().split()==['5','2']
  for _ in range(2):
   v=[0]*32
   for _ in range(int(f.readline())):
    c,z2,z1,z0,R,r=map(int,f.readline().split())
    v[r]=(v[r]+c*pow(base['Z2'],z2,P)*pow(base['Z1'],z1,P)*pow(base['Z0'],z0,P)*pow(base['R'],R,P))%P
   c1.append(v)
  assert not f.read().strip()
 assert c1==ac['C1_evaluated_coefficients_constant_then_linear']
 snum=scale(reduce(c1[0]),-1);sden=reduce(c1[1]);assert scale(snum,-1)==ac['C1_reduced_constant_B']and sden==ac['C1_reduced_linear_A']
 hp=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_HOMOGENIZATION_20261104_v63/certificate.json';h=read(hp)
 assert h['rational_coefficient_pairs']['E1']['constant']=='1/20'and h['rational_coefficient_pairs']['E2']['constant']=='5/4'
 def load(rec):return {(r,s,R,0):Q(int(a),int(b))for a,b,r,s,R in rec['coefficients']}
 j=load(h['factors']['J']);fn=load(h['factors']['Fn']);n1=load(h['rational_coefficient_pairs']['E1']['numerator']);n2=load(h['rational_coefficient_pairs']['E2']['numerator'])
 dc={(3,0,0,0):Q(1),(2,0,0,0):Q(-2),(1,0,0,0):Q(1)};z0={(0,0,0,1):Q(1)}
 # Expand the small universal s-polynomials FIRST, then homogenize directly.
 # This is distinct from the author's factorized JH/FH/NH construction.
 numerator=polyadd(polymul(polymul(polymul(z0,dc),j),n1),polyscale(polymul(polymul(z0,z0),polymul(fn,n2)),25))
 denominator=polyscale(polymul(polymul(dc,fn),polymul(j,j)),20)
 assert max(e[1]for e in numerator)<=4 and max(e[1]for e in denominator)==4
 rp=[power(rr,i)for i in range(max(e[0]for e in set(numerator)|set(denominator))+1)]
 sp=[power(snum,i)for i in range(5)];tp=[power(sden,i)for i in range(5)]
 def homogeneous(p,d):
  value=zero
  for(r,s,R,z),c in p.items():
   cc=mod(c)*pow(base['R'],R,P)*pow(base['Z0'],z,P)%P
   value=add(value,scale(mul(rp[r],mul(sp[s],tp[d-s])),cc))
  return value
 U=homogeneous(numerator,4);V=homogeneous(denominator,4);assert U==ac['U']and V==ac['V']
 for name,p,d in [('JH',j,1),('FH',fn,2),('NH1',n1,3),('NH2',n2,2)]:assert homogeneous(p,d)==ac['intermediate_residues'][name]
 ups=[power(U,i)for i in range(10)];vps=[power(V,i)for i in range(10)]
 cols=[mul(ups[i],vps[9-i])for i in range(10)];matrix=[[cols[j][i]for j in range(10)]for i in range(10)]
 assert matrix==ac['Krylov_matrix_rows']
 # Independent fraction-free Bareiss condensation modulo the same prime.
 a=[r[:]for r in matrix];sign=1;prev=1;pivots=[]
 for k in range(9):
  row=next(i for i in range(k,10)if a[i][k])
  if row!=k:a[row],a[k]=a[k],a[row];sign=-sign
  p=a[k][k];pivots.append(p);inv=pow(prev,-1,P)
  for i in range(k+1,10):
   for j in range(k+1,10):a[i][j]=(p*a[i][j]-a[i][k]*a[k][j])*inv%P
  for i in range(k+1,10):a[i][k]=0
  prev=p
 det=sign*a[9][9]%P;assert det==21349==ac['Krylov_det_mod32003']
 cert={'prime':P,'same_old_beta_point_index':0,'old_beta_modp':beta,'C':C,'D':D,'moment_transpose':M,
       'F4':f4,'pencil_coefficients_Newton_not_permutations':pencil,'base':base,'monicG':gg,'C1_constant_and_linear_raw':c1,
       'balanced_U':U,'balanced_V':V,'all100_matrix_entries':matrix,'independent_Bareiss_pivots':pivots,'determinant':det,
       'small_symbolic_s_pair_termcounts':[len(numerator),len(denominator)],
       'finite_free_model':'Z_(32003)[R,Z0,Z1,Z2] localized at the frozen residue maximal ideal and G-leading coefficient; monicG gives free rank10',
       'positive_nonvanishing_argument':'nonzero residue determinant implies nonzero characteristic-zero coefficient determinant; no emptiness or flatness lifting',
       'generic_V_nonzero':'the actual source chart and generic C1-linear recovery make r,r-1,J,Fn,C1_linear nonzero in the irreducible degree10 field',
       'special_V_unit_not_assumed':True,
       'primitive_element_argument':'generically divide every column by V^9; then 1,eta,...,eta^9 are independent, eta=U/V=ZN, so minimal degree10'}
 ep=out/'independent-certificate.json';ep.write_text(json.dumps(cert,indent=2)+'\n')
 inputs=[Path(__file__),ROOT/'research/notes'/(STEM+'.md'),folder/'report.json',folder/'certificate.json',pointp,spinp,gp,cp,hp,
         ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_CONTENT_PRIMITIVITY_20261104_v63_REVIEW/report.json',ROOT/'reports/Q_NEW_FAMILY_GENERAL_SOURCE_SCALE_20261104_v63_REVIEW/report.json']
 report={'id':STEM,'all_passed':True,'same_oldpoint_only':0,'prime_only':P,'no_new_beta_point_or_prime':True,
         'independent_original_moment_and_Newton_pencil_base_values':True,'independent_full_G_and_C1_specialization_reduction':True,
         'independent_expanded_degree4_balanced_U_V_without_special_denominator_inverse':True,'all100_finitefree_matrix_coefficients_match':True,
         'independent_Bareiss_determinant':det,'ZN_primitive_degree10_independently_proved':True,'special_V_unit_not_needed':True,
         'full_target_field_application_separate_no_NEW_equation_claim':True,
         'no_generic_eta_norm_PRS_resultant_source26_or_author_producer_rerun':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in inputs},'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started,'frozen':True}
 assert ep.stat().st_size<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k]for k in ['all_passed','independent_Bareiss_determinant','ZN_primitive_degree10_independently_proved','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/STEM);a=p.parse_args();main(a.output.resolve())
