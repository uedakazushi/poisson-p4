#!/usr/bin/env python3
"""Three coherent tauAB quartic residuals, compiled but never evaluated."""
import hashlib,json,math,signal,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C3_TAUAB_QUARTIC3_SOURCE_20261105_v63'
PATHS={
 'frame':'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63/actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json',
 'frame_peer':'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63_REVIEW/report.json',
 'B4':'reports/AMBIENT_B18_C3_TAUA_B4_SOURCE_V2_20261105_v63/four-actual-cubic-sources-two-cached-adjoints-and-rank2-completion.json',
 'B4_peer':'reports/AMBIENT_B18_C3_TAUA_B4_SOURCE_V2_20261105_v63_REVIEW/report.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'C2roots':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json',
 'C2_peer':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json',
 'raw33':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
 'ell':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63/all164-cached-pure-ell-readbacks-and-physical12-second-values.json',
 'ell_peer':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63_REVIEW/report.json',
 'right_source':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63/two-coherent-HB-quartic-sources-and-unpaired-affine-operator.json',
 'right_source_peer':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63_REVIEW/report.json',
 'right_roots':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_CACHED_VALUES_20261105_v63/phase1-six-root-existing-cache-readbacks.json',
 'right_value':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_CACHED_VALUES_20261105_v63/phase2-actual-RIGHT2-quartics-two-factors-and-rank.json',
 'right_value_peer':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_CACHED_VALUES_20261105_v63_REVIEW/report.json',
 'alpha':'reports/AMBIENT_B18_C3_ALPHA4_CACHED_VALUES_20261105_v63/eight-cached-even-root-values-and-four-complete-alpha-readbacks.json',
 'alpha_peer':'reports/AMBIENT_B18_C3_ALPHA4_CACHED_VALUES_20261105_v63_REVIEW/report.json',
 'graph_peer':'reports/AMBIENT_B18_C3_GRAPH5_FORCED_ZERO_20261105_v63_REVIEW/report.json',
 'cover':'reports/AMBIENT_B18_C3_GRAPH5_FORCED_ZERO_20261105_v63/actual-cover-NF-and-five-triangular-cubic-quartic-contradictions.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,(int,Q)):return 1
 if isinstance(x,str):
  try:Q(x);return 1
  except(ValueError,ZeroDivisionError):return 0
 if isinstance(x,(tuple,list)):return sum(map(count,x))
 if isinstance(x,dict):return sum(map(count,x.values()))
 return 0


class Pool:
 def __init__(self):self.nodes=[];self.seen={};self.inputs=[];self.c(Q(0));self.c(Q(1))
 def put(self,key):
  if key not in self.seen:self.seen[key]=len(self.nodes);self.nodes.append(list(key))
  return self.seen[key]
 def c(self,q):
  q=Q(q);assert max(abs(q.numerator).bit_length(),q.denominator.bit_length())<=2048
  return self.put(('c',q.numerator,q.denominator))
 def inp(self,key,value,path,locator):
  i=self.put(('input',tuple(key)))
  if not any(r['node']==i for r in self.inputs):
   self.inputs.append({'node':i,'key':list(key),'value':str(Q(value)),'source_path':path,'locator':locator})
  return i
 def add(self,a,b):
  if not a:return b
  if not b:return a
  if self.nodes[a][0]==self.nodes[b][0]=='c':
   x,y=[Q(*self.nodes[k][1:])for k in[a,b]]
   bound=max(max(abs(x.numerator).bit_length()+y.denominator.bit_length(),abs(y.numerator).bit_length()+x.denominator.bit_length())+1,x.denominator.bit_length()+y.denominator.bit_length());assert bound<=2048
   return self.c(x+y)
  return self.put(('a',*sorted([a,b])))
 def mul(self,a,b):
  if not a or not b:return 0
  if a==1:return b
  if b==1:return a
  if self.nodes[a][0]==self.nodes[b][0]=='c':
   x,y=[Q(*self.nodes[k][1:])for k in[a,b]]
   assert max(abs(x.numerator).bit_length()+abs(y.numerator).bit_length(),x.denominator.bit_length()+y.denominator.bit_length())<=2048
   return self.c(x*y)
  return self.put(('m',*sorted([a,b])))
 def neg(self,a):return self.mul(self.c(-1),a)
 def total(self,terms):
  n=0
  for a in terms:n=self.add(n,a)
  return n


def main():
 started=time.monotonic();signal.alarm(60);out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 pins={p:sha(ROOT/p)for p in PATHS.values()};d={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 assert all(v['all_passed']for k,v in d.items()if k.endswith('_peer'))
 assert d['right_value']['actual_RIGHT2_rank']==1
 v=list(map(Q,d['right_value']['rank_certificate']['pivot_column']));assert v==[Q(-3,5),Q(-121,80)]
 slope=Q(121,48);assert v[1]-slope*v[0]==0
 P=Pool();basis=d['raw33']['actual33_basis'];first=[basis.index([i+1,0,0])for i in range(7)]
 ir3,ir2l=basis.index([0,3,0]),basis.index([0,2,1]);ir2,irl=basis.index([0,2,0]),basis.index([0,1,1])
 izr=[basis.index([i+1,1,0])for i in range(7)];izl=[basis.index([i+1,0,1])for i in range(7)]
 c2={r['physical_row']:r for r in d['C2']['physical12_coefficients']};oldroots={i:dict(x)for i,x in d['C2roots']['all_requested164_root_values']}
 rows=[r['physical_row']for r in d['frame']['all12_verified_source_records']];pivot=[0,4,5];pivrows=[rows[i]for i in pivot]
 assert pivrows==[18,60,67]
 lambdas={r['physical_row']:r['delta_r_ell2']for r in d['right_roots']['all3_lambda_readbacks']}
 alpha={r['physical_row']:r['delta_r_ell2']for r in d['alpha']['all4_actual_alpha_delta_B6']};assert all(Q(v)==0 for v in alpha.values())
 betaQ=[[Q(x)for x in row]for row in d['B4']['beta2_plane_r2_rell_ell2']];active=[i for i,b in enumerate(betaQ)if any(b)];assert active==pivot
 bB2Q=[list(map(Q,r['bB_second_plane']))for r in d['frame']['all12_verified_source_records']]
 mon2=[(2,0),(1,1),(0,2)];mon3=[(3,0),(2,1),(1,2)];mon4=[(4,0),(3,1),(2,2)]
 def pad(a,b):
  c=a.copy()
  for m,n in b.items():
   z=P.add(c.get(m,0),n)
   if z:c[m]=z
   else:c.pop(m,None)
  return c
 def psc(a,q):return{m:P.mul(n,P.c(q))for m,n in a.items()if P.mul(n,P.c(q))}
 def pneg(a):return psc(a,Q(-1))
 def pmul(a,b):
  c={}
  for x,n in a.items():
   for y,k in b.items():
    m=(x[0]+y[0],x[1]+y[1])
    if m[1]<3 and sum(m)<=4:c[m]=P.add(c.get(m,0),P.mul(n,k))
  return{m:n for m,n in c.items()if n}
 def psum(polys):
  c={}
  for pp in polys:c=pad(c,pp)
  return c
 r={(1,0):1};ell={(0,1):1}
 beta=[{m:P.c(x)for m,x in zip(mon2,b)if x}for b in betaQ]
 def delta_input(physical,index):
  val=dict(c2[physical]['delta']).get(index,'0')
  return P.inp(('delta',physical,index),val,PATHS['C2'],['physical12_coefficients',physical,'delta_sparse33',index])
 def nu_input(side,index):
  rootid=[20086,40865][side];val=oldroots[rootid].get(index,'0')
  return P.inp(('nu_sparse33',rootid,index),val,PATHS['C2roots'],['all_requested164_root_values',rootid,index])
 nu=[]
 for side in range(2):
  assert Q(oldroots[[20086,40865][side]].get(0,'0'))==0
  assert Q(oldroots[[20086,40865][side]].get(basis.index([0,1,0]),'0'))==0
  assert Q(d['ell']['actual_nuA_B_pure_ell3'][side][1])==0
  e2=P.inp(('nu_ell2',side),d['ell']['actual_nuA_B_pure_ell3'][side][2],PATHS['ell'],['actual_nuA_B_pure_ell3',side,2])
  nu.append({mon2[0]:nu_input(side,ir2),mon2[1]:nu_input(side,irl),mon2[2]:e2})
 nuAB=pad(nu[0],pneg(psum([psc(beta[i],Q(1))if False else{m:P.mul(n,nu_input(0,first[i]))for m,n in beta[i].items()}for i in active])))
 D=[];E=[]
 for physical in pivrows:
  la=P.inp(('lambda',physical),lambdas[physical],PATHS['right_roots'],['all3_lambda_readbacks',physical,'delta_r_ell2'])
  dp={mon3[0]:delta_input(physical,ir3),mon3[1]:delta_input(physical,ir2l),mon3[2]:la}
  gamma=psum([pmul(pad({(1,0):delta_input(physical,izr[i])},{(0,1):delta_input(physical,izl[i])}),beta[i])for i in active])
  D.append(dp);E.append(pad(dp,pneg(gamma)))
 T=[list(map(Q,row))for row in d['right_source']['constant_source_normalization_T']]
 h=[list(map(Q,row))for row in d['right_source']['source_first7_h_rows']]
 H=[[row[i]for i in active]for row in h]
 HT=[[sum((H[j][i]*T[j][k]for j in range(3)),Q(0))for k in range(3)]for i in range(3)]
 PA=[psum([psc(D[i],T[i][k])for i in range(3)])for k in range(3)]
 PB=[psum([psc(E[i],HT[i][k])for i in range(3)])for k in range(3)]
 nf=[list(map(Q,row))[:3]for row in d['cover']['whole_NF4_rows2_columns15']]
 nf=[list(map(Q,d['cover']['two_block_source_curve_NF_receipts'][i]['NF4_to_r2ell2']))[:3]for i in range(2)]
 outputs=[];proofpolys=[]
 for physical,index in zip(pivrows,pivot):
  WA=pneg(pmul(nu[1],{m:P.c(x)for m,x in zip(mon2,bB2Q[index])if x}))
  WB=pneg(pmul(nuAB,beta[index]))
  WB=pad(WB,pneg(psum([pmul(pad({(1,0):delta_input(physical,izr[i])},{(0,1):delta_input(physical,izl[i])}),E[k])for k,i in enumerate(active)])))
  c0=pad({(1,0):delta_input(physical,ir3)},{(0,1):delta_input(physical,ir2l)})
  c1={(0,1):P.inp(('lambda',physical),lambdas[physical],PATHS['right_roots'],['all3_lambda_readbacks',physical,'delta_r_ell2'])}
  WA=pad(WA,pneg(pad(pmul(c0,PA[0]),pmul(c1,PA[1]))))
  WB=pad(WB,pneg(pad(pmul(c0,PB[0]),pmul(c1,PB[1]))))
  aa=P.total(P.mul(WA.get(m,0),P.c(q))for m,q in zip(mon4,nf[0]));bb=P.total(P.mul(WB.get(m,0),P.c(q))for m,q in zip(mon4,nf[1]))
  residual=P.add(bb,P.neg(P.mul(P.c(slope),aa)))
  outputs.append({'physical_row':physical,'NF4_A_node':aa,'NF4_B_node':bb,'residual_dual_node':residual})
  proofpolys.append({'physical_row':physical,'WA4_polynomial':[[list(m),n]for m,n in sorted(WA.items())],'WB4_polynomial':[[list(m),n]for m,n in sorted(WB.items())]})
 # All nine other tauAB quartics vanish by these exact accepted source
 # premises, not by missing higher coefficients being setzero.
 nonpivot=[]
 for j,physical in enumerate(rows):
  if j in pivot:continue
  vals=dict(c2[physical]['delta']);assert all(Q(x)==0 for x in bB2Q[j])
  if j<7:assert not beta[j]
  else:assert all(Q(x)==0 for x in d['frame']['all12_verified_source_records'][j]['bA_first_z7'])
  assert all(Q(vals.get(i,'0'))==0 for i in[ir3,ir2l]+[ii for i in active for ii in[izr[i],izl[i]]])
  nonpivot.append(physical)
 assert len(nonpivot)==9
 # Scalar interval sweep only; no signed target or polynomial coefficient
 # is evaluated. Constants and supplied input sizes are exact.
 inputs={r['node']:r for r in P.inputs};bounds=[];maxH=maxOp=0
 for i,node in enumerate(P.nodes):
  op=node[0]
  if op=='c':q=Q(node[1],node[2]);U,Dn=(abs(q.numerator)+q.denominator-1)//q.denominator,q.denominator;operation=max(abs(q.numerator).bit_length(),q.denominator.bit_length())
  elif op=='input':q=Q(inputs[i]['value']);U,Dn=(abs(q.numerator)+q.denominator-1)//q.denominator,q.denominator;operation=max(abs(q.numerator).bit_length(),q.denominator.bit_length())
  else:
   a,b=[bounds[k]for k in node[1:]]
   if op=='a':U,Dn=a[0]+b[0],math.lcm(a[1],b[1]);operation=2*max(a[2],b[2])+1
   elif op=='m':U,Dn=a[0]*b[0],a[1]*b[1];operation=a[2]+b[2]
   else:raise AssertionError(op)
  height=max((U*Dn).bit_length(),Dn.bit_length());maxH=max(maxH,height);maxOp=max(maxOp,operation);assert operation<=2048 and height<=2048
  bounds.append((U,Dn,height))
 input_count=count(d);source_bound=input_count+12*len(P.nodes)+10000
 assert len(P.nodes)<=3000 and source_bound<150000
 data={'nodes':P.nodes,'provided_inputs':P.inputs,'outputs3':outputs,'quartic_polynomial_nodes':proofpolys,
   'monomial_quotient':'ell^3=0 exactly, totaldegree<=4; omittedell3channels are cut relations, notunknowncoefficientzero',
   'all9_nonpivot_tauAB_quartic_classes_zero_source_rows':nonpivot,
   'known_quartic_vector':list(map(str,v)),'residual_dual_B_minus_121A_over48':str(slope),
   'nuA_B2_formula':'nuA_plane2-sum_i(nuA_zi*beta_i2)',
   'B_primitive_full_clear_formula':'-nuA_B2*u_j2-sum_i(delta_j2_zi*r*r+delta_j2_zi*ell*ell)*E_i3; this includes allbeta2^2 Bzclearing',
   'A_primitive_formula':'-nuB_plane2*bB_j2',
   'H_cubic_clear':'subtract (delta_j30*r+delta_j21*ell)*P0,3+lambda_j*ell*P1,3 in bothA/B; delta_j03*ell3H removedbyEXACTcut',
   'unknown_highest_annihilation':'z_iH has Bpart -E_i3 atorder3; allHδ2zz gives B>=5, Hδ3zmixed gives>=5, Hδ4clearsviaquadraticpivots withA/Bcarry>=5. Primitive A/Bquartics need onlynu/b secondjets andfirstz; beta3 enters>=5. Noδ4,mixedΦ3,β3orgaugecallback.',
   'firstframe_choice_scope':'Z_i,2=diag(0,-beta_i2,0) commute, so firstframecommutatorchoicesgivezeroA/Bquartic. AlternativeHcube decompositionsdifferbyRIGHT2HBsources, whoseknownrank1imageiskilledbytheresidualdual. Noallorderfreeplane-moduleclaim.',
   'source_bound_metadata':{'nodes':len(P.nodes),'provided_scalar_inputs':len(P.inputs),'max_scalar_height_bound':maxH,'max_operation_precheck_bound':maxOp,'source_owned_upper':source_bound,'source_hard_numeric':150000},
   'all_values_not_evaluated':True,'new_physical_source_IFT_gauge_jet':False,
   'positive_completion':'Anyresidualnonzero givesaquarticindependentoftheknownHBvector, henceC3length19/HF36640. All3zeroisONLYtauAB12subsetnegative,notlength20/fullquarticrank.'}
 cp=out/'three-tauAB-quartic-residual-scalar-circuit-and-source-proof.json';cp.write_text(json.dumps(data,indent=2)+'\n')
 bp=out/'all-scalar-majorants-no-value-evaluation.json';bp.write_text(json.dumps({'bounds_U_D_height':[[str(u),str(den),h]for u,den,h in bounds],'max_height':maxH,'max_precheck':maxOp},indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'SOURCE_ONLY_three_tauAB_quartic_residuals_NOT_EVALUATED','input_sha256':pins,
   'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'evidence_sha256':{p.name:sha(p)for p in[cp,bp]},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)


if __name__=='__main__':main()
