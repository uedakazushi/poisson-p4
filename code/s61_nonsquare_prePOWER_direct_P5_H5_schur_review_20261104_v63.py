#!/usr/bin/env python3
"""Independent directP5 source/type and noncommutative block audit."""
import hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='S61_NONSQUARE_PREPOWER_DIRECT_P5_H5_SCHUR_20261104_v63_REVIEW';OUT=ROOT/'reports'/STEM
PRIMARY='S61_NONSQUARE_PREPOWER_DIRECT_P5_H5_SCHUR_20261104_v63';P4='S61_NONSQUARE_PREPOWER_DIRECT_P4_H5_SCHUR_20261104_v63';ENTRY='S61_NONSQUARE_PREPOWER_U4_T15_SCHUR_PREFLIGHT_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+'\n')
# Free noncommutative block algebra. Only inverse pairs explicitly certified below reduce.
ONE={():F(1)};ZERO={};inverse={'Pi':'P','P':'Pi','Hi':'H','H':'Hi','Ki':'K','K':'Ki'}
def word(*xs):return {tuple(xs):F(1)}
def add(p,q,c=F(1)):
 out=dict(p)
 for w,v in q.items():
  a=out.get(w,F(0))+c*v
  if a:out[w]=a
  else:out.pop(w,None)
 return out
def mul(p,q):
 out={}
 for w,c in p.items():
  for v,d in q.items():
   stack=[]
   for x in w+v:
    if stack and inverse.get(stack[-1])==x:stack.pop()
    else:stack.append(x)
   a=tuple(stack);out=add(out,{a:c*d})
 return out
def mm(A,B):
 out=[]
 for row in A:
  rr=[]
  for j in range(len(B[0])):
   p={}
   for k,a in enumerate(row):p=add(p,mul(a,B[k][j]))
   rr.append(p)
  out.append(rr)
 return out
def eye(n):return [[ONE if i==j else ZERO for j in range(n)] for i in range(n)]
def neg(p):return {w:-v for w,v in p.items()}
def fixed_block_proof():
 P,Pi,H,Hi,C=map(lambda x:word(x),['P','Pi','H','Hi','C'])
 K=[[P,ZERO],[C,H]];Ki=[[Pi,ZERO],[neg(mul(mul(Hi,C),Pi)),Hi]]
 assert mm(K,Ki)==mm(Ki,K)==eye(2)
 K,Ki,B,C,D=map(lambda x:word(x),['K','Ki','B','C2','D'])
 A=[[K,B],[C,D]];L=[[Ki,ZERO],[neg(mul(C,Ki)),ONE]];Li=[[K,ZERO],[C,ONE]]
 U=[[ONE,neg(mul(Ki,B))],[ZERO,ONE]];Ui=[[ONE,mul(Ki,B)],[ZERO,ONE]]
 assert mm(L,Li)==mm(Li,L)==mm(U,Ui)==mm(Ui,U)==eye(2)
 assert mm(mm(L,A),U)==[[ONE,ZERO],[ZERO,add(D,mul(mul(C,Ki),B),-1)]]
 return True

def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(st,name='report.json'):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 primary=load(PRIMARY);assert pins[f'reports/{PRIMARY}/report.json']=='23bc50cba8f58ddee02b00ee982b39335a9ebca5a274d85edb45f103ddff013f'
 for p,h in primary['input_sha256'].items():assert sha(ROOT/p)==h
 op=load(PRIMARY,'operator.json');proof=load(PRIMARY,'proof.json');old=load(P4,'operator.json');pre=load(ENTRY,'operator-and-cost.json')
 assert load(P4+'_REVIEW')['all_passed'];assert load('S61_NONSQUARE_PREPOWER_FIFTH_PURE_PIVOT_20261104_v63_REVIEW')['all_passed'];assert load('S61_NONSQUARE_PREPOWER_U4_T15_MODP_POLYNOMIAL_20261104_v63_REVIEW')['all_passed']
 # Actual source cofactors retain the OLD raw rho convention.
 source=pre['nodes'];assert pre['source_current15_rows']==[5,2,4] and pre['source_coefficients3']==[33,35,37]
 assert pre['known_rho_LCM_powers']==[1,0,1,0,2,0]
 assert source[33]['op']=='mul' and source[33]['args']==[17,32]
 assert source[17]['reference']['key']=='literal_pair_determinant' and source[32]['reference']=={'constant':'-81/800000'}
 assert source[35]['op']=='mul' and source[35]['args']==[24,34] and source[34]['reference']=={'constant':'-9/40000'}
 assert source[37]['op']=='mul' and source[37]['args']==[31,36] and source[36]['reference']=={'constant':'-63/160000'}
 rho=load('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63','certificate.json');rhos=rho['known_row_units_rho14']
 for i,c in [(1,'40000/9'),(3,'160000/63'),(4,'-800000/81')]:assert rhos[i]=={'constant':c,'a_d_F_B_C_J_powers':[1,0,1,0,2,0]}
 assert all(not rho['raw_compact14x5'][i][j] for i in [1,3,4] for j in [1,2,3])
 # Cofactor formulas are the universal2x2 adjugate, BEFORE any P4 inverse.
 for i,j1,j2 in [(24,20,23),(31,27,30)]:assert source[i]['op']=='sub' and source[i]['args']==[j1,j2]
 assert [source[j]['reference']['index'] for j in [18,19,21,22,25,26,28,29]]==[[4,0],[3,4],[4,4],[3,0],[4,4],[1,0],[4,0],[1,4]]
 # Record full dependency closures, allowing only named FULLU unit or Delta5 inverses.
 invclosures={}
 for export in ['input15x22','hybrid_row_unit','heavy_column_unit','H5_inverse','source_current15_to_original25']:
  seen=set();todo=[old['exports'][export]]
  while todo:
   i=todo.pop()
   if i in seen:continue
   seen.add(i);todo+=old['nodes'][i]['args']
  rows=[]
  for i in sorted(seen):
   n=old['nodes'][i]
   if 'inverse' in n['op']:
    text=json.dumps(n.get('unit_proof'));assert all(x not in text for x in ['R434','h22','D4'])
    if isinstance(n.get('unit_proof'),dict):assert set(n['unit_proof'])=={'constant','a_d_F_B_C_J_powers'}
    else:assert 'FULLU' in text or 'Delta5' in text
    rows.append([i,n['op'],n.get('unit_proof')])
  invclosures[export]=rows
 # Type-check every operator and exact actual source/target orders.
 N=op['nodes'];assert len(N)==74
 for i,n in enumerate(N):
  assert n['id']==i and all(0<=j<i for j in n['args']);sh=n['shape'];arg=[N[j]['shape'] for j in n['args']];k=n['op']
  if k=='accepted_external_export':
   r=n['reference'];assert r['path']==f'reports/{P4}/operator.json' and r['sha256']==pins[r['path']] and r['node']==old['exports'][r['export']] and sh==old['nodes'][r['node']]['shape']
  elif k=='accepted_external_scalar_DAG_node':
   r=n['reference'];assert r['path']==f'reports/{ENTRY}/operator-and-cost.json' and r['sha256']==pins[r['path']] and r['node'] in [33,35,37] and sh==[1,1]
  elif k=='matrix_product':assert len(arg)==2 and arg[0][1]==arg[1][0] and sh==[arg[0][0],arg[1][1]]
  elif k=='submatrix':assert sh==[len(n['rows']),len(n['columns'])] and all(0<=j<arg[0][0] for j in n['rows']) and all(0<=j<arg[0][1] for j in n['columns'])
  elif k in ['negative','diagonal_known_unit_inverse']:assert len(arg)==1 and sh==arg[0]
  elif k=='central_scalar_multiple':assert arg[0]==[1,1] and sh==arg[1]
  elif k=='difference':assert arg==[sh,sh]
  elif k in ['identity','zero']:assert not arg
  elif k=='block_matrix':
   assert n['args']==[j for row in n['blocks'] for j in row]
   rr=[[N[j]['shape'] for j in row] for row in n['blocks']]
   assert all(len(set(z[0] for z in row))==1 for row in rr) and all(len(set(row[j][1] for row in rr))==1 for j in range(len(rr[0])))
   assert sh==[sum(row[0][0] for row in rr),sum(x[1] for x in rr[0])]
  elif k in ['row_permutation','column_permutation']:assert sorted(n['order'])==list(range(sh[0])) and sh[0]==sh[1]
  elif k=='transpose':assert sh==arg[0][::-1]
  elif k=='single_row_linear_replacement':assert sh==[15,15] and n['replaced_row']==5 and n['source_rows']==[5,2,4] and n['coefficients']==n['args'] and arg==[[1,1]]*3
  elif k=='certified_chart_scalar_inverse':assert sh==[1,1] and arg==[[1,1]] and n['chart'] in ['D(Delta5)','D(Saux)']
  elif k=='accepted_external_polynomial':assert sh==[1,1] and n['reference']['json_pointer']=='/exact_Saux_Laurent_polynomial' and n['terms']==7022
  elif k=='universal_polynomial_cofactor_adjugate':assert sh==[5,5] and arg==[[5,5]] and n['expanded'] is False
  else:raise AssertionError(k)
 # Independently carry exact Laurent row replacement15 and its inverse (commuting coefficients).
 # Encode exponent tuples in cu,c1,c3; standalone helpers imported from no producer.
 def cm(p,q):
  out={}
  for m,a in p.items():
   for n,b in q.items():
    w=tuple(x+y for x,y in zip(m,n));v=out.get(w,F(0))+a*b
    if v:out[w]=v
    else:out.pop(w,None)
  return out
 def ca(p,q):
  out=dict(p)
  for w,v in q.items():
   c=out.get(w,F(0))+v
   if c:out[w]=c
   else:out.pop(w,None)
  return out
 def matmul(A,B):
  return [[sumc([cm(A[i][k],B[k][j]) for k in range(len(B))]) for j in range(len(B[0]))] for i in range(len(A))]
 def sumc(ps):
  x={}
  for p in ps:x=ca(x,p)
  return x
 one={(0,0,0):F(1)};E=[[one if i==j else {} for j in range(15)] for i in range(15)];A=[[dict(p) for p in row] for row in E];Ai=[[dict(p) for p in row] for row in E]
 A[5]=[{} for _ in range(15)];Ai[5]=[{} for _ in range(15)]
 for j,m in [(5,(1,0,0)),(2,(0,1,0)),(4,(0,0,1))]:A[5][j]={m:F(1)}
 for j,m,c in [(5,(-1,0,0),1),(2,(-1,1,0),-1),(4,(-1,0,1),-1)]:Ai[5][j]={m:F(c)}
 assert matmul(A,Ai)==matmul(Ai,A)==E
 # Actual order and assembly are checked as finite syntax, not by specializing real coefficients.
 ctx=op['context'];assert ctx['pivot_current15']==[0,1,3,6,5,2,4,7,13,10] and ctx['remaining_current15']==[8,9,11,12,14]
 assert ctx['pivot_original32']==[22,11,13,25,15,27,31,28,30,29] and ctx['remaining_original32']==[1,2,3,4,6,7,8,9,14,17,18,19]
 assert ctx['remaining_original25']==[18,19,21,22,24] and ctx['all8_Q9_tails_retained']
 assert op['coefficient_ring']=='FULLU additionally localize ONLY at Delta5*Saux; R434,D4,h22 never inverted.'
 assert ctx['determinant']=='detK10=detP5*Delta5=Saux*Delta5' and ctx['new_scalar_values']==0
 exp=op['exports'];rep=N[exp['row_replacement15']];ri=N[exp['row_replacement15_inverse']]
 assert [N[j]['reference']['node'] for j in rep['args']]==[33,35,37]
 assert N[ri['args'][0]]['op']=='certified_chart_scalar_inverse' and N[ri['args'][0]]['args']==[rep['args'][0]]
 # Wire-level reconstruction of the actual block recipe, beyond shape checks.
 def flat(i):
  return sum((flat(j) for j in N[i]['args']),[]) if N[i]['op']=='matrix_product' else [i]
 def sub(i,parent,rr,cc):
  q=N[i];assert q['op']=='submatrix' and q['args']==[parent] and q['rows']==rr and q['columns']==cc
 def one_ref(i,name):assert N[i]['op']=='accepted_external_export' and N[i]['reference']['export']==name
 P5node=exp['P5'];H5node=exp['H5'];MM=N[P5node]['args'][0]
 f=flat(MM);assert len(f)==6;Pr,Sunit,Rnode,Mnode,Tunit,Pc=f
 assert Rnode==exp['row_replacement15'];one_ref(Mnode,'input15x22');one_ref(Sunit,'hybrid_row_unit');one_ref(Tunit,'heavy_column_unit')
 assert N[Pr]['order']==ctx['pivot_current15']+ctx['remaining_current15']
 oldcols=[1,2,3,4,6,7,8,9,11,13,14,15,17,18,19,22,25,27,28,29,30,31]
 assert [oldcols[j] for j in N[Pc]['order']]==ctx['pivot_original32']+ctx['remaining_original32']
 sub(P5node,MM,list(range(5)),list(range(5)));sub(H5node,MM,list(range(5,10)),list(range(5,10)))
 sub(exp['pivot10'],MM,list(range(10)),list(range(10)))
 ki=N[exp['pivot10_inverse']];assert ki['op']=='block_matrix' and len(ki['blocks'])==2
 pi,z,lower,hi=ki['args'];one_ref(hi,'H5_inverse');assert N[z]['op']=='zero'
 assert N[pi]['op']=='central_scalar_multiple';isax,adj=N[pi]['args'];assert N[adj]['op']=='universal_polynomial_cofactor_adjugate' and N[adj]['args']==[P5node]
 assert N[isax]['op']=='certified_chart_scalar_inverse' and N[isax]['chart']=='D(Saux)'
 assert N[lower]['op']=='negative';lprod=flat(N[lower]['args'][0]);assert lprod[0]==hi and lprod[2]==pi
 C55=lprod[1];sub(C55,MM,list(range(5,10)),list(range(5)))
 res=N[exp['residual5x12']];assert res['op']=='difference';Dnode,correction=res['args'];sub(Dnode,MM,list(range(10,15)),list(range(10,22)))
 ckb=flat(correction);assert len(ckb)==3 and ckb[1]==exp['pivot10_inverse'];Cnode,Bnode=ckb[0],ckb[2]
 sub(Cnode,MM,list(range(10,15)),list(range(10)));sub(Bnode,MM,list(range(10)),list(range(10,22)))
 fs=flat(exp['full_source15']);assert fs[1:]==[Pr,Sunit,Rnode]
 ft=flat(exp['full_target22']);assert ft[:2]==[Tunit,Pc]
 Lnode,U=N[fs[0]],N[ft[2]];assert Lnode['op']==U['op']=='block_matrix'
 assert Lnode['blocks'][0][0]==exp['pivot10_inverse']
 minusCK=Lnode['blocks'][1][0];assert N[minusCK]['op']=='negative' and flat(N[minusCK]['args'][0])==[Cnode,exp['pivot10_inverse']]
 minusKB=U['blocks'][0][1];assert N[minusKB]['op']=='negative' and flat(N[minusKB]['args'][0])==[exp['pivot10_inverse'],Bnode]
 fsi=flat(exp['full_source15_inverse']);assert len(fsi)==4 and fsi[0]==exp['row_replacement15_inverse']
 assert N[fsi[1]]['op']=='diagonal_known_unit_inverse' and N[fsi[1]]['args']==[Sunit]
 assert N[fsi[2]]['op']=='transpose' and N[fsi[2]]['args']==[Pr]
 assert N[fsi[3]]['op']=='block_matrix' and N[fsi[3]]['blocks'][0][0]==exp['pivot10'] and N[fsi[3]]['blocks'][1][0]==Cnode
 fti=flat(exp['full_target22_inverse']);assert len(fti)==3 and N[fti[1]]['op']=='transpose' and N[fti[1]]['args']==[Pc] and N[fti[2]]['op']=='diagonal_known_unit_inverse' and N[fti[2]]['args']==[Tunit]
 assert N[fti[0]]['op']=='block_matrix' and flat(N[fti[0]]['blocks'][0][1])==[exp['pivot10_inverse'],Bnode]
 surv=flat(exp['surviving5_source_original25']);assert len(surv)==2;sub(surv[0],exp['full_source15'],list(range(10,15)),list(range(15)));one_ref(surv[1],'source_current15_to_original25')
 assert fixed_block_proof()
 # P5 determinant is a polynomial bordered identity; P4 inverse never enters any node.
 assert proof['determinant_polynomial_proof'].startswith('Universal bordered determinant identity') and proof['no_R434_inverse_in_H5_or_prefix']
 save('independent-P5-source-block-proof.json',{'typed_nodes74':True,'actual_raw_cofactors3_bound':pre['source_coefficients3'],'cu_exact':'(-81/800000)*delta2','cu_is_unit_on_Delta5':True,'raw_heavy_zero_by_pair_adj2_and_literal_middle3_zeros':True,'row_replacement_two_inverse_products15':True,'all_inverse_dependency_closures':invclosures,'P5_bordered_determinant_Saux_without_D4_inverse':True,'adj5_polynomial_unexpanded':True,'noncommuting_K10_and_full_Schur_block_inverse_identities':True,'remaining5x12_source_target_orders_and8tails':True,'chart':'FULLU D(Delta5*Saux)','R434_zero_not_discarded_no_nonempty_intersection_claim':True,'no_actual5x12_coefficient_or_rank5_evaluation':True,'7022_QQ_coefficients_not_newly_reexpanded':True})
 report={'id':STEM,'all_passed':True,'scope':'actual source-bound directP5/H5 typed Schur, no scalar coefficient job','actual_samecokernel_shape':[5,12],'parameter_localizers':['Delta5','Saux'],'R434_D4_h22_not_inverted':True,'all8_Q9_pure_tails_retained':True,'source_row_replacement_and_inverse_certified':True,'new_scalar_coefficients':0,'global57_quartic239_or_rank5_not_claimed':True,'no_producer_import_rerun_point_prime_GB_or_scalar_search':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{'independent-P5-source-block-proof.json':sha(OUT/'independent-P5-source-block-proof.json')},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};save('report.json',report);print(json.dumps(report))
if __name__=='__main__':main()
