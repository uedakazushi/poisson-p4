#!/usr/bin/env python3
"""Global H² unit deletion and COMPLETE H-degree1 support circuit on the new122 target."""
import argparse,copy,hashlib,itertools,json,resource,signal,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from s61_nonsquare_joint_Lell_free183_20261104_v63 import Circuit
from s61_nonsquare_free281_circuit_20261104_v63 import canon
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def native(p):
 ls=p.read_text().splitlines();r,c=map(int,ls[0].split());out=[]
 for line in ls[1:]:
  t=line.split();assert len(t)==1+4*int(t[0]);out.append([[list(map(int,t[i+1:i+4])),t[i]]for i in range(1,len(t),4)])
 assert len(out)==r*c;return[[out[i*c+j]for j in range(c)]for i in range(r)]
class Polynomial:
 """An exact weighted polynomial circuit over an existing scalar circuit."""
 def __init__(self,o):self.o=o;self.nodes=[['z']];self.ids={('z',):0};self.coefficient_cache={}
 def node(self,key,value):
  if key not in self.ids:self.ids[key]=len(self.nodes);self.nodes.append(value)
  return self.ids[key]
 def var(self,i):return self.node(('v',i),['v',i])
 def const(self,c):return 0 if not c else self.node(('c',c),['c',c])
 def mul(self,a,b):return 0 if not a or not b else self.node(('m',a,b),['m',a,b])
 def lin(self,terms):
  pairs=tuple((c,n)for c,n in terms if c and n)
  if not pairs:return 0
  if len(pairs)==1 and pairs[0][0]==self.o.one:return pairs[0][1]
  return self.node(('l',pairs),['l',[list(p)for p in pairs]])
 def coefficient(self,i,mon):
  """Exact finite extraction; PURE98 queries are deliberately not invoked by this preflight."""
  mon=tuple(mon);key=(i,mon)
  if key in self.coefficient_cache:return self.coefficient_cache[key]
  if not i:z=0
  else:
   n=self.nodes[i]
   if n[0]=='c':z=n[1]if not any(mon)else 0
   elif n[0]=='v':z=self.o.one if mon==tuple(int(k==n[1])for k in range(8))else 0
   elif n[0]=='l':z=self.o.summands((c,self.coefficient(j,mon))for c,j in n[1])
   else:
    assert n[0]=='m';z=0
    for left in itertools.product(*(range(k+1)for k in mon)):
     right=tuple(a-b for a,b in zip(mon,left));z=self.o.add(z,self.o.mul(self.coefficient(n[1],left),self.coefficient(n[2],right)))
  self.coefficient_cache[key]=z;return z
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 fb=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63';fr=read(fb/'report.json');fd=read(fb/'certificate.json');assert sha(fb/'certificate.json')==fr['certificate_sha256']
 freview=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63_REVIEW/report.json';assert read(freview)['all_passed']
 jb=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_FREE183_20261104_v63';jd=read(jb/'certificate.json');olddefs=jd['old_source319_definitions'];mons=jd['new_target239_monomials'];mi={tuple(m):j for j,m in enumerate(mons)};kept122=fd['remaining122_original_joint_target_columns']
 NF=fd['full_target_projection239_to122'];assert len(NF['output122_linear_roots'])==122
 o=Circuit(fd['scalar_circuit']);o.l=copy.deepcopy(fd['linear_nodes']);o.li={('l',tuple(tuple(t)for t in n[1]))if n[0]=='l'else tuple(n):i for i,n in enumerate(o.l)}
 prefix_scalar=len(o.nodes);prefix_linear=len(o.l);one=o.one;R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field();enc=lambda p:[[list(m),str(c)]for m,c in p.items()]
 def atom(c):c=K(c);return o.const(enc(c.numer),enc(c.denom))
 def parse(c):return K.from_expr(s.sympify(c))
 weightsH=[mons[j][7]for j in kept122];classes={h:[i for i,v in enumerate(weightsH)if v==h]for h in range(3)};assert list(map(len,[classes[0],classes[1],classes[2]]))==[98,23,1]
 # An exact structural support check: none of the117 quadratic pivots mixes H degree.
 support_cache={0:set()}
 def support(i):
  if i in support_cache:return support_cache[i]
  n=o.l[i]
  z={n[1]}if n[0]=='v'else set().union(*(support(j)for c,j in n[1]if c))
  support_cache[i]=z;return z
 var_index={v:j for j,v in enumerate(NF['input_variables'])}
 for i,root in enumerate(NF['output122_linear_roots']):assert all(mons[var_index[v]][7]==weightsH[i]for v in support(root))
 h2col=classes[2][0];h2input=NF['input_variables'][kept122[h2col]];assert o.l[NF['output122_linear_roots'][h2col]]==['v',h2input]
 Hroots=[NF['output122_linear_roots'][i]for i in classes[1]]
 quartic_ret=[olddefs[g][1]for g in fd['Cnew108_original319_labels']if olddefs[g][0]=='quartic'];assert len(quartic_ret)==22 and 77 in quartic_ret and 131 in quartic_ret
 C108=fd['Cnew108_original319_labels'];source319pair=[next(i for i,row in enumerate(olddefs)if row==['quartic',r])for r in[77,131]];assert source319pair==[303,318]
 assert all(olddefs[g][0]=='cubic_normalized'for g in fd['Bnew22_original319_labels'])
 # Locate literal surviving pair generators through the ACTUAL free165->164 coordinatechange/deletion.
 retained164=[i for i in range(165)if i!=fd['separator165_coordinates'][0]];pos164=[]
 for g in source319pair:
  pos165=57+C108.index(g);assert pos165 not in fd['separator165_coordinates'];pos164.append(retained164.index(pos165))
 # Reuse the already computed924 coefficient roots; no old query or source generation is repeated.
 hb=ROOT/'reports/S61_NONSQUARE_QUARTIC_HSYM2_PREFLIGHT_20261103_v63';hr=read(hb/'report.json');hd=read(hb/'preflight.json');assert sha(hb/'preflight.json')==hr['preflight_sha256']and hr['Hsym2_queries']==924
 opath=ROOT/'reports/S61_NONSQUARE_LAZY_QUARTIC_ORACLE_20261103_v63/oracle.json';oracle=read(opath);arith=hd['arithmetic'];rawarith=oracle['arithmetic'];assert arith['nodes'][:len(rawarith['nodes'])]==rawarith['nodes']and arith['atoms'][:len(rawarith['atoms'])]==rawarith['atoms']
 oreview=ROOT/'reports/S61_NONSQUARE_LAZY_QUARTIC_ORACLE_20261103_v63_REVIEW/report.json';assert read(oreview)['all_passed']
 adpath=ROOT/'reports/S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63/adapter.json';ad=read(adpath);adreview=ROOT/'reports/S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63_REVIEW/report.json';assert read(adreview)['all_passed']
 Hcache={0:0}
 def hs(i):
  if i in Hcache:return Hcache[i]
  assert arith['node_weights'][i]==0
  n=arith['nodes'][i]
  if n[0]=='c':v=o.const(*arith['atoms'][n[1]])
  elif n[0]=='m':v=o.mul(hs(n[1]),hs(n[2]))
  else:assert n[0]=='l';v=o.summands((o.const(*arith['atoms'][c]),hs(j))for c,j in n[1])
  Hcache[i]=v;return v
 qrows=hd['source_rows'];assert qrows==oracle['genuine_quartic_original_rows'];pairmons=hd['H_times_light_pair_monomials'];light=[0,1,2,3,4,6,7];D=list(map(parse,ad['Tden']));DH=D[5]
 rawH2={r:hs(i)for r,i in hd['raw_H2_coefficient_nodes']};rawH={r:[hs(i)for i in row]for r,row in zip(qrows,hd['raw33_by28_coefficient_nodes'])}
 oldm2=[r[2]for r in olddefs[:29]];oldmi={tuple(m):i for i,m in enumerate(oldm2)};newm2=[m for qi,m in jd['new_prime_quadric_source174_definitions'][:29]];TB=[dict(row)for row in jd['monomial_old_to_new_scalar_roots']]
 def Hproject(newquad):
  values={v:0 for v in NF['input_variables']}
  for m,c in zip(newm2,newquad):
   if c:
    assert m[7]==0;mm=m[:];mm[7]=1;values[NF['input_variables'][mi[tuple(mm)]]]=c
  return o.linear_apply_scalar(Hroots,values)
 rawQuarticProjected={}
 for row in qrows:
  newquad=[0]*29
  for k,monH in enumerate(pairmons):
   mon=monH[:];mon[5]=0;oldcol=oldmi[tuple(mon)];i,j=hd['light_pair_columns'][k];scaled=o.mul(rawH[row][k],atom(D[i]*D[j]/DH))
   for newcol,c in TB[oldcol].items():newquad[newcol]=o.add(newquad[newcol],o.mul(scaled,c))
  rawQuarticProjected[row]=Hproject(newquad)
 # Full source77/131 H² coefficients are exactly the old literal values: jointH unchanged andNF preservesHdegree.
 hpath=ROOT/'reports/S61_BINARY_HYPERPLANE_HEAVY_SQUARE_PIVOT_20261102.json';h=read(hpath);u,v=map(parse,h['row_combination_coefficients']);b77,b131=map(parse,h['coefficient_formulas']);B=(50*a+18*e**2-9*e-2)*d+60*e-40;C=100*a+18*e**2-9*e-2;U=K((50*a*B*C)**2);assert u*b77+v*b131==U
 source_matrix=[[u,v],[-b131,b77]];source_inverse=[[b77/U,-v/U],[b131/U,u/U]];assert all(sum((source_matrix[i][k]*source_inverse[k][j]for k in range(2)),K.zero)==(K.one if i==j else K.zero)for i in range(2)for j in range(2))
 ur,vr,b0r,b1r,Ur,Ui=list(map(atom,[u,v,b77,b131,U,1/U]));chiH=[o.add(o.mul(ur,x),o.mul(vr,y))for x,y in zip(rawQuarticProjected[77],rawQuarticProjected[131])];psiH=[o.add(o.mul(o.sub(0,b1r),x),o.mul(b0r,y))for x,y in zip(rawQuarticProjected[77],rawQuarticProjected[131])]
 quartic21labels=[['psi',77,131]]+[r for r in quartic_ret if r not in[77,131]];quartic21=[psiH]
 for row in quartic_ret:
  if row in[77,131]:continue
  quartic21.append([o.sub(x,o.mul(o.mul(rawH2[row],Ui),y))for x,y in zip(rawQuarticProjected[row],chiH)])
 assert len(quartic21)==21
 # Genuine16 cubic H-linear coefficients are the SAME frozenTn-normalizedsource57; no extra adapter factor.
 cpath=ROOT/'reports/S61_NONSQUARE_QUARTIC_122_CIRCUIT_20261103_v63/circuit.json';cc=read(cpath);cmons=read(ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json')['normal_monomials'];heavy=[];heavyterms=0
 heavy_by_basis={row[0]['basis_row']:row for row in cc['heavy_coefficient_entries']};assert sorted(heavy_by_basis)==list(range(41,57))
 for i in range(16):
  row=heavy_by_basis[41+i]
  rr=[]
  for j,entry in enumerate(row):
   assert entry['basis_row']==41+i;mon=cmons[entry['column']];assert mon[5]==1 and mon[light[j]]==1 and sum(mon)==2
   t=entry['sparse_encoding'].split();n=int(t[0]);assert len(t)==1+4*n;heavyterms+=n;rr.append(o.const([[list(map(int,t[k+1:k+4])),t[k]]for k in range(1,len(t),4)],[[[0,0,0],'1']]))
  heavy.append(rr)
 # Recreate only the light inverse entries, not any full cubic/quartic coefficient vector.
 pre=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63';N=native(pre/'coordinate-inverse-numerators7x7.txt');delta=K(3200*a*C*d**2)/81;NN=[[o.mul(o.const(z,[[[0,0,0],'1']]),atom(1/delta))for z in row]for row in N]
 oldB=[NN[0],NN[1],NN[2],NN[6],NN[3],[o.sub(o.add(NN[5][j],NN[1][j]),o.mul(atom(3),NN[0][j]))for j in range(7)],NN[4]]
 hnew=[[o.summands((row[j],oldB[j][k])for j in range(7))for k in range(7)]for row in heavy];newmi={tuple(m):i for i,m in enumerate(newm2)};cubics=[]
 for i in range(16):
  for oldj in range(7):
   coeff=[0]*29
   for k in range(7):
    for l in range(7):
     m=[0]*8;m[k]+=1;m[l]+=1;idx=newmi[tuple(m)];coeff[idx]=o.add(coeff[idx],o.mul(hnew[i][k],oldB[oldj][l]))
   cubics.append(Hproject(coeff))
 assert len(cubics)==112
 # The only remaining canonical quadric generator with H degree1 is q5*H.
 qb=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_QPAIR_UNITS_20261104_v63';qs=read(qb/'checkpoint.json');assert 57 in qs['remaining_rows']
 pairdata=read(ROOT/'reports/S61_NONSQUARE_JOINT_LELL_OVERLAP_UNIT_20261104_v63/remaining-pair-block.json');assert pairdata['source_rows'][57]==[5,[0,0,0,0,0,0,0,1]]
 assert kept122==[jd['retained141_target_columns'][pairdata['target141_remaining_columns'][i]]for i in qs['remaining_columns']]
 quadH=[o.const(*qs['matrix58x136'][57][qs['remaining_columns'][i]])for i in classes[1]]
 for i in classes[0]+classes[2]:assert not o.const(*qs['matrix58x136'][57][qs['remaining_columns'][i]])
 block=[quadH]+cubics+quartic21;assert len(block)==134 and all(len(row)==23 for row in block)
 # A whole-polynomial circuit for chi, leaving every pure98 coefficient deferred and source-bound.
 P=Polynomial(o);newvar=[P.var(i)for i in range(8)];varimages=[]
 for old in range(8):
  if old==5:varimages.append(P.lin([(atom(DH),newvar[7])]))
  else:j=light.index(old);varimages.append(P.lin((o.mul(atom(D[old]),oldB[j][k]),newvar[k])for k in range(7)))
 pcache={0:0}
 def pn(i):
  if i in pcache:return pcache[i]
  n=rawarith['nodes'][i]
  if n[0]=='v':z=varimages[n[1]]
  elif n[0]=='c':z=P.const(o.const(*rawarith['atoms'][n[1]]))
  elif n[0]=='m':z=P.mul(pn(n[1]),pn(n[2]))
  else:assert n[0]=='l';z=P.lin((o.const(*rawarith['atoms'][c]),pn(j))for c,j in n[1])
  pcache[i]=z;return z
 roots=dict(zip(oracle['genuine_quartic_original_rows'],oracle['genuine_quartic_roots']));assert rawarith['node_weights'][roots[77]]==rawarith['node_weights'][roots[131]]==4
 q77=P.lin([(atom(1/DH**2),pn(roots[77]))]);q131=P.lin([(atom(1/DH**2),pn(roots[131]))]);chiPoly=P.lin([(ur,q77),(vr,q131)])
 # No coefficient(P, puremonomial) call has occurred: the lazyoperator retains the exact missingpure98 image.
 assert not P.coefficient_cache
 kept121=[i for i in range(122)if i!=h2col];source163=[i for i in range(164)if i!=pos164[0]]
 record={'coefficient_domain':'FULLU','target122_Hdegree_classes':classes,'target117_projection_preserves_Hdegree_exact_structural_support':True,
 'literal_free164_pair_positions_for77_131':pos164,'literal_old319_pair_positions':[303,318],'neither_pairrow_deleted_or_mixed_by_previous_sourceoperations':True,
 'H2_source_pair_change_scalar_roots':[[ur,vr],[o.sub(0,b1r),b0r]],'H2_source_pair_inverse_scalar_roots':[[atom(c)for c in row]for row in source_inverse],
 'H2_literal_coefficient_roots_from_reviewed_formulas':[b0r,b1r],'H2_known_unit_root':Ur,'H2_known_unit':'(50*a*B*C)^2','source_pair_determinant_equals_knownunit_exact_QQ_checked':True,
 'H2_pivot_target122_column':h2col,'free163_source164_changed_indices':source163,'target121_original122_columns':kept121,'actual_free163_shape':[163,121],
 'actual_free163_recipe':'Inactualfree164firstreplacepair77/131bychi=uQ77+vQ131andpsi=-b131Q77+b77Q131(detUunit). Projecteveryrowto121columnsbyrow_D-(row_H2/U)*chi_D,thendeletechi. Sameactual319cokernel; no fixedold163submatrix.',
 'chi_full_joint_weight4_polynomial_root':chiPoly,'chi_polynomial_nodes':P.nodes,'chi_raw77_Qhat_polynomial_root':q77,'chi_raw131_Qhat_polynomial_root':q131,
 'chi_lazy_post117_coefficient_recipe':'UsePolynomial.coefficient(chi_root,joint_monomial)andtheacceptedfull_target_projection239_to122linearoperator. Thisisfiniteexactweightedpolynomialextraction,notapure98expansionhere.',
 'chi_H1_post117_scalar_roots':chiH,'chi_H2_post117_known_scalar_root':Ur,'pure98_coefficient_queries_performed':0,
 'retained_quartic_original_rows22':quartic_ret,'all_remaining21_quartic_H2free_labels':quartic21labels,'all_remaining21_quartic_H1_scalar_roots':quartic21,
 'all33_raw_quartic_H1_post117_scalar_roots':[[r,rawQuarticProjected[r]]for r in qrows],'all33_raw_H2_coefficient_scalar_roots':[[r,rawH2[r]]for r in qrows],
 'quarticT34_toTn_rule':'ForH*oldlight_i*oldlight_j multiplythefrozenT34coefficientbyD_i*D_j/D_H;H2coefficientsunchanged;thenjointcoordinateandNF117.',
 'genuine_cubic16_Hlinear_old_scalar_roots':heavy,'genuine_cubic16_Hlinear_joint_scalar_roots':hnew,'all112_genuine_cubic_H1_scalar_roots':cubics,'all112_cubic_source319_definitions':[row for row in olddefs if row[0]=='cubic_normalized'],
 'remaining_q5H_H1_scalar_roots':quadH,'complete_H1_image_block_shape':[134,23],'complete_H1_image_block_scalar_roots':block,
 'H1_rowspan_equivalence_proof':'Afterold15andP57complementthefullgeneratorimagesarethesameR-spanasoriginalquad/cubic/quarticgenerators. Hgradingisstrict. Onlyq5Hsurvivesinthecanonicalquadblock,all112cubicsarekept(redundant4areexactoldrelations),all22retainedquarticsundergoaninvertible77/131changeandH2pivot,leavingall21. Thereforethis134-rowblockspanstheENTIREH1projectionoffree163,notjusta112-rowsufficientsubset.',
 'H1_support_upper_bounds':[sum(bool(x)for x in row)for row in block],'H1_scalar_DAG_nonzero_root_is_not_nonvanishing_certificate':True,
 'no_H1_rank_Gaussian_or_allpoint_surjectivity_claim':True,'allpoint_full239_and_actualpure98_image_open':True,
 'input_heavy112_coefficient_terms':heavyterms,'reused_quartic_Hsym2_coefficient_queries':924,
 'scalar_circuit':{'atoms_num_den':o.atoms,'nodes':o.nodes},'linear_nodes':o.l,
 'parent_target239_to122_operator':NF,'parent_joint_target239_monomials':mons,'parent_target122_monomial_indices':kept122}
 (out/'certificate.json').write_text(json.dumps(record,separators=(',',':'))+'\n');assert len(o.nodes)+len(o.l)+len(P.nodes)<200000 and(out/'certificate.json').stat().st_size<32*1024*1024
 paths=[Path(__file__),ROOT/'code/s61_nonsquare_joint_Lell_free183_20261104_v63.py',ROOT/'research/notes'/(ID+'.md'),fb/'report.json',fb/'certificate.json',freview,jb/'certificate.json',hb/'report.json',hb/'preflight.json',opath,oreview,adpath,adreview,hpath,cpath,pre/'coordinate-inverse-numerators7x7.txt',qb/'checkpoint.json']
 report={'id':ID,'global_H2_knownunit_sourcepair_and_deletion_verified':True,'actual_free_source_target_shape':[163,121],'target_after_H2_split':{'pure_quartic':98,'H_times_light_quadric':23},
 'complete_H1_projection_shape':[134,23],'all112_cubics_and_all21_remaining_quartics_included':True,'one_residual_quadH_row_included':True,
 'reused_all924_quartic_coefficient_queries_without_regeneration':True,'T34_toTn_Hcoefficient_scaling_explicit':True,'pure98_coefficient_queries':0,'no_H1_rank_Gaussian_or_full163matrix_expansion':True,
 'no_rank121_full239_or_new_uniform_kernel_count_claim':True,'scalar_prefix_nodes':prefix_scalar,'linear_prefix_nodes':prefix_linear,'scalar_nodes':len(o.nodes),'linear_nodes':len(o.l),'weighted_polynomial_nodes':len(P.nodes),'core_atoms':len(o.atoms),
 'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'certificate_sha256':sha(out/'certificate.json'),'artifact_bytes':(out/'certificate.json').stat().st_size,'seconds':time.monotonic()-began,'frozen':True,'independent_review_pending':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
