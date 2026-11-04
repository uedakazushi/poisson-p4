#!/usr/bin/env python3
"""Lossless monic98 Schur and known38 transport, as finite scalar/linear circuits."""
import argparse,copy,hashlib,json,resource,signal,time
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
from sympy.polys.rings import ring
from s61_nonsquare_free281_circuit_20261104_v63 import DAG,canon
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_JOINT_LELL_FREE183_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def readmat(p):
 ls=p.read_text().splitlines();r,c=map(int,ls[0].split());assert len(ls)==1+r*c;out=[]
 for line in ls[1:]:
  t=line.split();n=int(t[0]);assert len(t)==1+4*n;out.append([[list(map(int,t[i+1:i+4])),t[i]]for i in range(1,len(t),4)])
 return[[out[i*c+j]for j in range(c)]for i in range(r)]
class Circuit(DAG):
 def __init__(self,prefix):
  self.atoms=copy.deepcopy(prefix['atoms_num_den']);self.nodes=copy.deepcopy(prefix['nodes']);self.atom_ids={json.dumps(v,separators=(',',':')):i for i,v in enumerate(self.atoms)};self.node_ids={tuple(v):i for i,v in enumerate(self.nodes)};self.one=self.const([[[0,0,0],'1']],[[[0,0,0],'1']]);self.l=[['z']];self.li={('z',):0}
 def lv(self,i):
  key=('v',i)
  if key not in self.li:self.li[key]=len(self.l);self.l.append(['v',i])
  return self.li[key]
 def lin(self,terms):
  z={}
  for c,n in terms:
   if c and n:z[n]=self.add(z.get(n,0),c)
  if not z:return 0
  if len(z)==1:
   n,c=next(iter(z.items()))
   if c==self.one:return n
  ts=tuple((c,n)for n,c in sorted(z.items()));key=('l',ts)
  if key not in self.li:self.li[key]=len(self.l);self.l.append(['l',[list(t)for t in ts]])
  return self.li[key]
 def linear_apply_scalar(self,roots,values):
  cache={0:0}
  def value(i):
   if i in cache:return cache[i]
   nd=self.l[i]
   if nd[0]=='v':v=values[nd[1]]
   else:assert nd[0]=='l';v=self.summands((c,value(j))for c,j in nd[1])
   cache[i]=v;return v
  return[value(i)for i in roots]
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 pre=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63';pr=read(pre/'report.json');assert pr['all_exact_joint_coordinate_and_quadric_checks_passed']
 for name,h in pr['evidence_sha256'].items():assert sha(pre/name)==h
 pflight=read(pre/'preflight.json');stage=pflight['stages'][1];mons=stage['target_monomials'];mi={tuple(m):j for j,m in enumerate(mons)};newm2=[m for i,m in stage['four_quadric_product_rows'][:29]];assert len(newm2)==29
 qbase=ROOT/'reports/S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63/certificate.json';qd=read(qbase);olddefs=[qd['source_rows'][i]for i in qd['retained_row_indices']];oldm2=[d[2]for d in olddefs[:29]]
 fbase=ROOT/'reports/S61_NONSQUARE_FREE281_CIRCUIT_20261104_v63';fr=read(fbase/'report.json');fd=read(fbase/'certificate.json');freview=ROOT/'reports/S61_NONSQUARE_FREE281_CIRCUIT_20261104_v63_REVIEW/report.json';assert read(freview)['all_passed']and fr['same_actual319_cokernel']
 assert sha(fbase/'certificate.json')==fr['certificate_sha256'];o=Circuit(fd['circuit']);prefix_nodes=len(o.nodes)
 R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field();enc=lambda p:[[list(m),str(c)]for m,c in p.items()]
 def atom(c):c=K(c);return o.const(enc(c.numer),enc(c.denom))
 one=o.one;negative=atom(-1);zeroPoly=[[[0,0,0],'1']]
 C=100*a+18*e**2-9*e-2;delta=K(3200)/K(81)*K(a*C*d**2)
 Nraw=readmat(pre/'coordinate-inverse-numerators7x7.txt');Fraw=readmat(pre/'coordinate-forward7x7.txt');Qraw=readmat(pre/'six-prime-quadrics-numerators.txt');D=readmat(pre/'six-prime-quadrics-denominators.txt')
 # delta has rational coefficients, and enc(delta.numer)/enc(delta.denom) may have a constant denominator.
 # Use exactly Nraw/delta, not Nraw/delta.numer.
 N=[[o.mul(o.const(p,zeroPoly),atom(1/delta))for p in row]for row in Nraw]
 F=[[o.const(p,zeroPoly)for p in row]for row in Fraw];quad=[[o.const(Qraw[i][j],D[i][0])for j in range(29)]for i in range(6)]
 oldlight=[0,1,2,3,4,6,7]
 # Images of oldv in new y and of new y in oldv. H stays the heavy coordinate7 on the new side.
 B=[[0]*8 for _ in range(8)];A=[[0]*8 for _ in range(8)]
 for v,x in[(0,0),(1,1),(2,2),(3,6),(4,3),(7,4)]:
  for j in range(7):B[v][j]=N[x][j]
 for j in range(7):B[6][j]=o.sub(o.add(N[5][j],N[1][j]),o.mul(atom(3),N[0][j]))
 B[5][7]=one
 for y in range(7):
  A[y][0]=o.add(F[y][0],o.mul(atom(3),F[y][5]));A[y][1]=o.sub(F[y][1],F[y][5]);A[y][2]=F[y][2];A[y][3]=F[y][6];A[y][4]=F[y][3];A[y][6]=F[y][5];A[y][7]=F[y][4]
 A[7][5]=one
 def monomial_map(m,images):
  p={(0,)*8:one}
  for v,n in enumerate(m):
   for _ in range(n):
    q={}
    for ex,c in p.items():
     for j,z in enumerate(images[v]):
      if z:
       mm=list(ex);mm[j]+=1;mm=tuple(mm);q[mm]=o.add(q.get(mm,0),o.mul(c,z))
    p=q
  return p
 newmi={tuple(m):j for j,m in enumerate(newm2)};oldmi={tuple(m):j for j,m in enumerate(oldm2)}
 TB=[{newmi[m]:c for m,c in monomial_map(m,B).items()}for m in oldm2]
 TA=[{oldmi[m]:c for m,c in monomial_map(m,A).items()}for m in newm2]
 # qprime_i=sum_j Change[i,j]*qold_j. Its explicit inverse has polynomial entries.
 change=[[K.zero]*6 for _ in range(6)];inv=[[K.zero]*6 for _ in range(6)]
 change[0][1]=K(8100);change[1][1]=K(81*e)/K(4*a);change[1][3]=-K(81)/K(4*a);change[2][1]=K(10*a+3*e**2-3*e)/K(a);change[2][3]=-K(3*(e-1))/K(a);change[2][4]=K(-10);change[3][1]=-K(d*e+8);change[3][2]=K(8);change[3][3]=K(d);change[4][0]=change[5][5]=K.one
 inv[0][4]=inv[5][5]=K.one;inv[1][0]=K.one/8100;inv[3][0]=K(e)/8100;inv[3][1]=-K(4*a)/81;inv[2][0]=K.one/8100;inv[2][1]=K(a*d)/162;inv[2][3]=K.one/8;inv[4][0]=K.one/8100;inv[4][1]=K(2*(e-1))/135;inv[4][2]=-K.one/10
 assert all(sum((change[i][k]*inv[k][j]for k in range(6)),K.zero)==(K.one if i==j else K.zero)for i in range(6)for j in range(6))
 cm=[[atom(c)for c in row]for row in change];im=[[atom(c)for c in row]for row in inv]
 def transform(values,direction,linear):
  # Cofactor contragredient: old_to_prime uses inverseq and oldv(new); prime_to_old uses changeq and new(oldv).
  monmap=TB if direction=='old_to_prime'else TA;coeff=im if direction=='old_to_prime'else cm;summ=o.lin if linear else o.summands
  tmp=[[0]*29 for _ in range(6)]
  for qi in range(6):
   terms=[[]for _ in range(29)]
   for j,row in enumerate(monmap):
    for k,c in row.items():terms[k].append((c,values[qi*29+j]))
   tmp[qi]=[summ(ts)for ts in terms]
  return[summ((coeff[old][new],tmp[old][j])for old in range(6)if coeff[old][new])for new in range(6)for j in range(29)]
 # All174 prime quadric product generators, exact coefficient roots (no polynomial parameter expansion).
 qrows=[]
 for qi in range(6):
  for m in newm2:
   row={}
   for n,c in zip(newm2,quad[qi]):
    if c:row[mi[tuple(x+y for x,y in zip(m,n))]]=c
   qrows.append(row)
 selected=stage['monic_pivot_source_rows'];targets=stage['monic_pivot_target_columns'];targetD=stage['retained_target_columns'];assert len(selected)==len(targets)==98 and len(targetD)==141
 rest=[i for i in range(174)if i not in selected];assert len(rest)==76 and all(i in rest for i in stage['all_unselected_overlap_rows_retained'])
 U=[[qrows[r].get(c,0)for c in targets]for r in selected];RD=[[qrows[r].get(c,0)for c in targetD]for r in selected];NC=[[qrows[r].get(c,0)for c in targets]for r in rest]
 assert all(U[i][i]==one for i in range(98))and all(U[i][j]==0 for i in range(98)for j in range(i))
 # Generic linear operators use disjoint variable ranges; this keeps large K/P compositions factored.
 old_input=list(range(174));prime_input=list(range(174,348));rest_input=list(range(348,424));k_input=list(range(424,446));target_input=list(range(446,685))
 old2prime=transform([o.lv(i)for i in old_input],'old_to_prime',True);prime2old=transform([o.lv(i)for i in prime_input],'prime_to_old',True)
 xlin=[o.lv(i)for i in rest_input];b=[o.lin((NC[e][j],xlin[e])for e in range(76)if NC[e][j])for j in range(98)];beta=[]
 for j in range(98):beta.append(o.lin([(negative,b[j])]+[(o.sub(0,U[i][j]),beta[i])for i in range(j)if U[i][j]]))
 restored=[0]*174
 for i,r in enumerate(selected):restored[r]=beta[i]
 for i,r in enumerate(rest):restored[r]=xlin[i]
 restoredold=transform(restored,'prime_to_old',True)
 pbase=ROOT/'reports/S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63';Pd=read(pbase/'global-P22.json');pc=Pd['circuit'];sc={0:0};lc={0:0}
 def ps(i):
  if i in sc:return sc[i]
  nd=pc['scalar_nodes'][i]
  if nd[0]=='c':v=o.const(*pc['atoms_num_den'][nd[1]])
  elif nd[0]=='a':v=o.add(ps(nd[1]),ps(nd[2]))
  elif nd[0]=='s':v=o.sub(ps(nd[1]),ps(nd[2]))
  else:assert nd[0]=='m';v=o.mul(ps(nd[1]),ps(nd[2]))
  sc[i]=v;return v
 def pl(i):
  if i in lc:return lc[i]
  nd=pc['linear_nodes'][i]
  if nd[0]=='v':v=restoredold[nd[1]]
  else:assert nd[0]=='l';v=o.lin((ps(c),pl(j))for c,j in nd[1])
  lc[i]=v;return v
 Prest=[pl(i)for i in Pd['P22_linear_output_roots_in_frozen_original_kernel_order']]
 Kpath=ROOT/'reports/S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63/kernel22.json';Kd=read(Kpath);kr=[[]for _ in range(174)]
 for j,col in enumerate(Kd['sparse_columns']):
  for i,n,de in col:kr[i].append((o.const(n,de),o.lv(k_input[j])))
 Koldtimes=[o.lin(ts)for ts in kr];KprimeTimes=transform(Koldtimes,'old_to_prime',True);KrestTimes=[KprimeTimes[i]for i in rest]
 # Whole source319 known38 transport: ONLY its quadratic basis changes. Cubic/33quartic source labels remain original.
 oldx=fd['new_quadratic174_roots'];xprime=o.linear_apply_scalar(old2prime,dict(zip(old_input,oldx)));x76=[xprime[i]for i in rest]
 u22=fd['new_B22_roots'];h108=fd['new_C108_roots'];Px=o.linear_apply_scalar(Prest,dict(zip(rest_input,x76)));k22=[o.sub(Px[j],u22[j])for j in range(22)]
 Kk=o.linear_apply_scalar(KrestTimes,dict(zip(k_input,k22)));z76=[o.sub(x76[j],Kk[j])for j in range(76)];kernel184=z76+h108
 Bnew=fd['source304_new_block_B22'];Cnew=fd['source304_new_block_C108'];p,q=[76+Cnew.index(g)for g in[182,186]];assert[p,q]==[76,77]
 theta=fd['known_unit_separator_root'];sep=o.sub(o.add(kernel184[p],kernel184[q]),theta);assert sep==fd['separator_identity_zero_root']and fd['separator_identity_exact_QQ_checked']
 w=kernel184[:];w[p]=theta;norm=[o.mul(v,atom(-K(243)/(K(1024000)*K(a*C**3*(5*a*d-6*e-4)))))for v in w];norm[p]=one
 # Target normal-form operator, in239independent formal coefficient inputs.
 c=[o.lv(target_input[j])for j in targets];mult=[]
 for j in range(98):mult.append(o.lin([(one,c[j])]+[(o.sub(0,U[i][j]),mult[i])for i in range(j)if U[i][j]]))
 NF=[o.lin([(one,o.lv(target_input[d]))]+[(o.sub(0,RD[i][j]),mult[i])for i in range(98)if RD[i][j]])for j,d in enumerate(targetD)]
 record={'coefficient_domain':'FULLU=QQ[a,e,d,1/(a*d*F*B*C*J)]','old_source319_definitions':olddefs,
 'new_target239_monomials':mons,'new_prime_quadric_source174_definitions':[[i,m]for i in range(6)for m in newm2],
 'only_quadratic174_source_basis_changed':True,'cubic112_keep_old_light_multipliers':True,'all33_quartic_source_generators_retained_before_exact15_deletions':True,
 'quadric_generator_coefficient_roots6x29':quad,'quadric_product_coefficient_sparse_rows174':[[[j,c]for j,c in sorted(row.items())]for row in qrows],
 'unit98_source_rows':selected,'unit98_target_columns':targets,'retained76_quadratic_source_rows':rest,'retained141_target_columns':targetD,
 'all18_overlap_source_rows_retained':stage['all_unselected_overlap_rows_retained'],
 'U98_scalar_roots':U,'RD98x141_scalar_roots':RD,'NC76x98_scalar_roots':NC,'U_exact_unit_upper_triangular':True,
 'target_projection_operator':{'input_variables':target_input,'output141_linear_roots':NF,'recipe':'row_D-row_C*U^{-1}*R_D; triangular recursion; no unknown-factor division'},
 'old_to_prime_quad_cofactor_operator':{'input_variables':old_input,'output174_linear_roots':old2prime},
 'prime_to_old_quad_cofactor_operator':{'input_variables':prime_input,'output174_linear_roots':prime2old},
 'quad_restoration_operator':{'input_variables':rest_input,'prime174_linear_roots':restored,'old174_linear_roots':restoredold,'recipe':'beta_R=-U^{-T}*NC^T*x76; beta_E=x76'},
 'descended_P22_operator':{'input_variables':rest_input,'output22_linear_roots':Prest,'recipe':'P22_old after oldquad restoration'},
 'descended_K22_product_operator':{'input_variables':k_input,'output76_linear_roots':KrestTimes,'recipe':'restrictE*old_to_prime*K22_old'},
 'descended_PK_proof':'The pure-quadric identity M*K22_old=0 and unitU imply restore(restrictE*old_to_prime*K22_old)=old_to_prime*K22_old. Restoring the sourcebasis gives K22_old; P22_old*K22_old=I22, hence Pnew*Knew=I22.',
 'coordinate_inverse_proof_dependency':str((pre/'report.json').relative_to(ROOT)),'six_quadric_change_K6':[[str(c)for c in row]for row in change],'inverse_six_quadric_change':[[str(c)for c in row]for row in inv],
 'monomial_old_to_new_scalar_roots':[[[j,c]for j,c in sorted(row.items())]for row in TB],'monomial_new_to_old_scalar_roots':[[[j,c]for j,c in sorted(row.items())]for row in TA],
 'old15_transport':'Restrict the transformed quadratic sourcecofactors toE76; all nonquadratic coefficients are untouched. Same quartic11 and cubic4 unit pivots delete15 source rows after the98 target pivots.',
 'Bnew22_original319_labels':Bnew,'Cnew108_original319_labels':Cnew,'new_quad76_kernel_after_old15_roots':x76,'new_B22_kernel_after_old15_roots':u22,'new_C108_kernel_after_old15_roots':h108,
 'Pnew_on_actualx76_roots':Px,'kernel22_coordinate_roots':k22,'Knew_times_actualk22_roots':Kk,
 'actual184_kernel_roots':kernel184,'separator184_indices':[p,q],'separator_known_unit_root':theta,'separator_known_unit':'-1024000*a*C^3*J/243','separator_zero_check_root':sep,
 'separator_identity_is_literal_same_roots_as_reviewed_free281':True,'changed184_kernel_roots':w,'normalized_unit_kernel_roots':norm,
 'actual_free183_shape':[183,141],'source183_changed184_indices':[i for i in range(184)if i!=p],
 'actual_free183_generator_recipe':'On the exact221x141 Schur presentation, deleteold15 by their unit pivots, givingquad76+surplus130. Split surplus130=Bnew22+Cnew108 and use[A76+Bnew*Pnew,Cnew]free184. Change generator_q to generator_q-generator_p, deletep using the savedtheta-unitkernel. No fixed183oldsubmatrix claim.',
 'source319_lazy_target_oracle':'Apply the old acceptedTn/Qhat source319 polynomial generator to oldv(newy)=inverseN/delta, Hunchanged. Quadratic rows use savedprimequadrics/newweightedmonomials; cubic112 use oldgenuineG_i times OLDv_j evaluatedatnewy; all33quartic rows use the old full135 oracle. Apply the explicit NF239->141 operator. Never discardS/overlap rows or higherG6.',
 'rank_equivalence':'Full actual319 rank239 iff the new free183->141 map is surjective; coker equal on FULLU. Surjectivity remains unproved.',
 'allpoint_rank141_full239_and_genericbeta_exhaustion_open':True,'no_new_point_field_GB_matrix_rank_or_expanded_parametric_quartics':True,
 'scalar_circuit':{'atoms_num_den':o.atoms,'nodes':o.nodes},'linear_nodes':o.l}
 (out/'certificate.json').write_text(json.dumps(record,separators=(',',':'))+'\n')
 assert len(o.nodes)+len(o.l)<200000 and (out/'certificate.json').stat().st_size<32*1024*1024
 inputs=[Path(__file__),ROOT/'code/s61_nonsquare_free281_circuit_20261104_v63.py',ROOT/'research/notes'/(ID+'.md'),pre/'report.json',pre/'preflight.json',pre/'circuit.json',pre/'coordinate-inverse-numerators7x7.txt',pre/'coordinate-forward7x7.txt',pre/'six-prime-quadrics-numerators.txt',pre/'six-prime-quadrics-denominators.txt',fbase/'report.json',fbase/'certificate.json',freview,pbase/'global-P22.json',Kpath,qbase,ROOT/'reports/S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63/adapter.json']
 report={'id':ID,'finite_exact_unit98_and_known38_transport_constructed':True,'actual_lossless_Schur_shape':[221,141],'actual_free_source_target_shape':[183,141],
 'U98_unit_upper_triangular_exactly_checked':True,'all18_overlap_rows_and_all33_quartic_generators_retained_at_Schur_stage':True,
 'only_quad_source174_basis_changed_nonquadratic_pivots_preserved':True,'descended_P22K22_identity_by_exact_restoration_proof':True,
 'actual_source319_cokernel_preserved_on_FULLU':True,'no_rank141_or_full239_claim':True,'known38_missing42_unchanged':True,
 'scalar_prefix_nodes_reused':prefix_nodes,'scalar_nodes':len(o.nodes),'linear_nodes':len(o.l),'core_atoms':len(o.atoms),
 'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in inputs},'certificate_sha256':sha(out/'certificate.json'),'artifact_bytes':(out/'certificate.json').stat().st_size,'seconds':time.monotonic()-began,'frozen':True,'independent_review_pending':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
