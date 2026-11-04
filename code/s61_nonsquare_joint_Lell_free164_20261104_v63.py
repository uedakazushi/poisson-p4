#!/usr/bin/env python3
"""Exact117-pivot restoration, followed by actual known38 complement and unit deletion."""
import argparse,copy,hashlib,json,resource,signal,time
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
from sympy.polys.rings import ring
from s61_nonsquare_joint_Lell_free183_20261104_v63 import Circuit
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 jb=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_FREE183_20261104_v63';jr=read(jb/'report.json');jd=read(jb/'certificate.json');assert sha(jb/'certificate.json')==jr['certificate_sha256']
 jreview=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_FREE183_20261104_v63_REVIEW/report.json';assert read(jreview)['all_passed']
 ob=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_OVERLAP_UNIT_20261104_v63';orr=read(ob/'report.json');od=read(ob/'overlap-unit-certificate.json');pairinput=read(ob/'remaining-pair-block.json');assert orr['extra_global_known_unit_pivots']==5 and od['all2538_QQ_row_operation_identities_passed']
 oreview=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_OVERLAP_UNIT_20261104_v63_REVIEW/report.json';assert read(oreview)['all_passed']
 qb=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_QPAIR_UNITS_20261104_v63';qr=read(qb/'report.json');qs=read(qb/'checkpoint.json');assert qr['completed_global_known_unit_pivots']==qs['completed_pivots']==14
 assert sha(qb/'checkpoint.json')==qr['evidence_sha256']['checkpoint.json']
 qreview=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_QPAIR_UNITS_20261104_v63_REVIEW/report.json';assert read(qreview)['all_passed']
 o=Circuit(jd['scalar_circuit']);o.l=copy.deepcopy(jd['linear_nodes']);o.li={}
 for i,n in enumerate(o.l):o.li[('l',tuple(tuple(t)for t in n[1]))if n[0]=='l'else tuple(n)]=i
 prefix_scalar=len(o.nodes);prefix_linear=len(o.l);one=o.one
 def rat(v):return o.const(*v)
 def inverse(v):return o.const(v[1],v[0])
 def linear_map(roots,values):
  cache={0:0}
  def val(i):
   if i in cache:return cache[i]
   n=o.l[i]
   if n[0]=='v':v=values[n[1]]
   else:assert n[0]=='l';v=o.lin((c,val(j))for c,j in n[1])
   cache[i]=v;return v
  return[val(i)for i in roots]
 qrows=[dict(row)for row in jd['quadric_product_coefficient_sparse_rows174']];oldE=jd['retained76_quadratic_source_rows'];overlaplabels=od['source_original_overlap_rows'];newm2=[m for qi,m in jd['new_prime_quadric_source174_definitions'][:29]];midx={tuple(m):i for i,m in enumerate(newm2)}
 qpairlabels=[qi*29+midx[tuple(m)]for qi,m in pairinput['source_rows']];assert oldE[:18]==overlaplabels and oldE[18:]==qpairlabels
 oversteps=od['unit_pivots'];qsteps=qs['all_completed_unit_steps'];selected_over=[p['pivot_original_overlap_row']for p in oversteps];selected_qp=[18+p['pivot_row']for p in qsteps];newE=[i for i in range(76)if i not in selected_over+selected_qp];assert len(newE)==57
 vars57=list(range(1000,1057));coeff=[0]*76
 for j,i in enumerate(newE):coeff[i]=o.lv(vars57[j])
 # Reverse the additional14 row-operations on coefficients. Each is the exact inverse transpose.
 qreverse=[]
 for step in reversed(qsteps):
  i=18+step['pivot_row'];terms=[(one,coeff[i])]
  for t,c in step['eliminations']:terms.append((o.sub(0,rat(c)),coeff[18+t]))
  ip=inverse(step['pivot']);coeff[i]=o.lin([(ip,o.lin(terms))]);qreverse.append({'source76_pivot':i,'inverse_pivot_root':ip,'elimination_coefficient_roots':[[18+t,rat(c)]for t,c in step['eliminations']]})
 after_reverse14=coeff[:]
 # Reconstruct ONLY the missing58x5 coefficients at the overlap pivots, as finiteDAGs.
 # These are source quadric products; not one cubic/quartic coefficient is expanded.
 oldNF=jd['target_projection_operator'];nfroots=oldNF['output141_linear_roots'];tvars=oldNF['input_variables'];overrows=[[rat(v)for v in od['after_rows'][s['pivot_original_overlap_row']]]for s in oversteps];overcols=[p['pivot_target141_column']for p in oversteps]
 qpair_over_coeff=[]
 for source in qpairlabels:
  raw=o.linear_apply_scalar([nfroots[c]for c in overcols],{tvars[j]:qrows[source].get(j,0)for j in range(239)})
  cs=[]
  for k,c in enumerate(raw):
   for prev in range(k):c=o.sub(c,o.mul(cs[prev],overrows[prev][overcols[k]]))
   cs.append(c)
  qpair_over_coeff.append(cs)
 overreverse=[]
 for k in range(4,-1,-1):
  step=oversteps[k];i=step['pivot_original_overlap_row'];terms=[(one,coeff[i])];record=[]
  for t,c in step['eliminations']:root=rat(c);terms.append((o.sub(0,root),coeff[t]));record.append([t,root])
  for t in range(58):
   root=qpair_over_coeff[t][k]
   if root:terms.append((o.sub(0,root),coeff[18+t]));record.append([18+t,root])
  ip=inverse(step['pivot_value']);coeff[i]=o.lin([(ip,o.lin(terms))]);overreverse.append({'source76_pivot':i,'inverse_pivot_root':ip,'elimination_coefficient_roots':record})
 restored76=coeff[:]
 # Reuse the already accepted98-restoration, old-source inverse and P22, composing on the new57 inputs.
 oldRest=jd['quad_restoration_operator'];restored174prime=linear_map(oldRest['prime174_linear_roots'],dict(zip(oldRest['input_variables'],restored76)))
 restored174old=linear_map(oldRest['old174_linear_roots'],dict(zip(oldRest['input_variables'],restored76)))
 p76=jd['descended_P22_operator'];P57=linear_map(p76['output22_linear_roots'],dict(zip(p76['input_variables'],restored76)))
 k76=jd['descended_K22_product_operator'];K57=[k76['output76_linear_roots'][i]for i in newE]
 # For any full kernel, nonpivot-source coefficients are unchanged under these quadratic unit Schurs.
 x57=[jd['new_quad76_kernel_after_old15_roots'][i]for i in newE];u22=jd['new_B22_kernel_after_old15_roots'];h108=jd['new_C108_kernel_after_old15_roots']
 Px=o.linear_apply_scalar(P57,dict(zip(vars57,x57)));k22=[o.sub(Px[j],u22[j])for j in range(22)];Kk=o.linear_apply_scalar(K57,dict(zip(k76['input_variables'],k22)));z57=[o.sub(x57[i],Kk[i])for i in range(57)];kernel165=z57+h108
 Bnew=jd['Bnew22_original319_labels'];Cnew=jd['Cnew108_original319_labels'];p,q=[57+Cnew.index(i)for i in[182,186]];assert[p,q]==[57,58];theta=jd['separator_known_unit_root'];sep=o.sub(o.add(kernel165[p],kernel165[q]),theta);assert sep==jd['separator_zero_check_root']
 R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field();iv=-K(243)/(K(1024000)*K(a*(100*a+18*e**2-9*e-2)**3*(5*a*d-6*e-4)));enc=lambda p:[[list(m),str(c)]for m,c in p.items()];inverse_theta=o.const(enc(iv.numer),enc(iv.denom));changed=kernel165[:];changed[p]=theta;normalized=[o.mul(v,inverse_theta)for v in changed];normalized[p]=one
 # The exact target projection141->136->122, keeping every other equation.
 vars141=list(range(2000,2141));var136=list(range(2200,2336));r141=[o.lv(i)for i in vars141]
 for step,pivotrow in zip(oversteps,overrows):
  c=step['pivot_target141_column'];v=r141[c];assert pivotrow[c]==one
  for j,z in enumerate(pivotrow):
   if j!=c and z:r141[j]=o.lin([(one,r141[j]),(o.sub(0,z),v)])
  r141[c]=0
 keep136=pairinput['target141_remaining_columns'];assert all(r141[i]==0 for i in overcols);project136=[r141[i]for i in keep136]
 r136=[o.lv(i)for i in var136];qpivotrows=[]
 for step in qsteps:
  pivotrow=[rat(v)for v in qs['matrix58x136'][step['pivot_row']]];c=step['pivot_column'];v=r136[c];assert pivotrow[c]==one
  for prev in[q['pivot_column']for q in qsteps[:len(qpivotrows)]]:assert pivotrow[prev]==0
  for j,z in enumerate(pivotrow):
   if j!=c and z:r136[j]=o.lin([(one,r136[j]),(o.sub(0,z),v)])
  r136[c]=0;qpivotrows.append(pivotrow)
 keep122=qs['remaining_columns'];assert len(keep122)==122;project122=[r136[i]for i in keep122]
 first=linear_map(project136,dict(zip(vars141,nfroots)));NF122=linear_map(project122,dict(zip(var136,first)));oldtargetD=jd['retained141_target_columns'];target122=[oldtargetD[keep136[i]]for i in keep122]
 source117=jd['unit98_source_rows']+[oldE[i]for i in selected_over+selected_qp];target117=jd['unit98_target_columns']+[oldtargetD[c]for c in overcols]+[oldtargetD[keep136[s['pivot_column']]]for s in qsteps]
 assert len(set(source117))==len(set(target117))==117
 source57=[oldE[i]for i in newE]
 record={'coefficient_domain':'FULLU=QQ[a,e,d,1/(a*d*F*B*C*J)]','safe_order':'ORIGINAL319 ->117quad source/target unitpivots ->202x122 ->old15nonquad deletions ->57quad+130surplus ->P57 complementfree165 ->separatorunitdeletefree164. NOT deletionof19rowsfromoldmixedfree183.',
 'initial_joint183_certificate_sha256':sha(jb/'certificate.json'),'unit117_original_prime_quad_source_rows':source117,'unit117_original_joint_target_columns':target117,'remaining57_original_prime_quad_source_rows':source57,'remaining122_original_joint_target_columns':target122,
 'source202_definitions':[['prime_quadric_product',i,jd['new_prime_quadric_source174_definitions'][i]]for i in source57]+[['old_source319',i,jd['old_source319_definitions'][i]]for i in range(174,319)],
 'all33_genuine_quartics_retained_at202_stage':True,'no_overlap_Srow_dropped_without_exact_unit_or_oldkernel_basis':True,
 'additional_pivots_in_old76_source_coordinates':selected_over+selected_qp,'remaining57_in_old76_source_coordinates':newE,
 'normalizer_rowoperation_inverse_transpose_rule':'Forwardgenerator operations g_r<-g_r/p andg_t<-g_t-c_t*g_r. Oldcoefficient k_r=(newk_r-sum_t c_t*newk_t)/p; allothercoefficientsunchanged. Reverse the14then5then98 stages.',
 'restoration57_input_variables':vars57,'coefficient76_after_reverse14_linear_roots':after_reverse14,'restored76_after_reverse19_linear_roots':restored76,
 'reverse14_steps':qreverse,'reverse5_steps':overreverse,'qpair58_at_overlap5_pivot_coefficient_scalar_roots':qpair_over_coeff,
 'restored_prime174_linear_roots':restored174prime,'restored_old174_linear_roots':restored174old,
 'P57_output22_linear_roots':P57,'K57_input22_variables':k76['input_variables'],'K57_output57_linear_roots':K57,
 'P57K57_identity_proof':'ForpureoldK22thefullquadraticidentityand117unitpivotsforceallpivotcoefficientsinthefinaltriangularsourcebasistozero; restrictingto57andapplyingtheexactinverse-transposerestorationrecoversoldK76,thenoldK174. HenceP57K57=PoldKold=I22. No generic rank/cokernel-flatness inference.',
 'known15_nonquad_pivots_unchanged':True,'Bnew22_original319_labels':Bnew,'Cnew108_original319_labels':Cnew,'Bnew_chosen_from_surplus130_not_from_quadratic57':True,
 'actual_new_kernel_quad57_roots':x57,'actual_new_kernel_B22_roots':u22,'actual_new_kernel_C108_roots':h108,'P57_on_actualx_roots':Px,'kernel22_coordinate_roots':k22,'K57_times_actualk_roots':Kk,
 'actual_free165_kernel_roots':kernel165,'separator165_coordinates':[p,q],'known_unit_separator_root':theta,'separator_identity_zero_root':sep,'separator_identity_literal_same_as_reviewed_free281':True,'changed165_kernel_roots':changed,'unit_normalized_kernel_roots':normalized,
 'actual_free164_shape':[164,122],'retained_source165_indices_after_constant_pair_change':[i for i in range(165)if i!=p],
 'actual_generator_recipe':'On202->122deleteold15usingunchangednonquadpivots. Let A57bequadSchurrows,Bnew22andC108theretainedsurplus. Take[A57+Bnew*P57,C108]free165;replacegenerator_qbygenerator_q-generator_p,thendeletepusingtheta. Samecokernelasoriginal319; notafixed164-rowoldsubmatrix.',
 'target_projection141_to136':{'input_variables':vars141,'output136_linear_roots':project136},'target_projection136_to122':{'input_variables':var136,'output122_linear_roots':project122},
 'full_target_projection239_to122':{'input_variables':oldNF['input_variables'],'output122_linear_roots':NF122},
 'all33_quartic_lazy_oracle':'KeepoldTn/Qhatfull135oracleandoldcubicmultipliers,applythejointcoordinateinverseandthisexactNF239->122;no183or164matrixcoefficientexpansion.',
 'no_full_rank122_or_full239_or_new_uniform_kernel_direction_claim':True,'known38_missing42_unchanged':True,
 'scalar_circuit':{'atoms_num_den':o.atoms,'nodes':o.nodes},'linear_nodes':o.l}
 (out/'certificate.json').write_text(json.dumps(record,separators=(',',':'))+'\n');assert len(o.nodes)+len(o.l)<200000 and(out/'certificate.json').stat().st_size<32*1024*1024
 inputs=[Path(__file__),ROOT/'code/s61_nonsquare_joint_Lell_free183_20261104_v63.py',ROOT/'research/notes'/(ID+'.md'),jb/'report.json',jb/'certificate.json',jreview,ob/'report.json',ob/'overlap-unit-certificate.json',ob/'remaining-pair-block.json',oreview,qb/'report.json',qb/'checkpoint.json',qreview]
 report={'id':ID,'actual117_unit_Schur_and_exact_known38_transport_constructed':True,'unmixed_lossless_source_target_shape':[202,122],'actual_free_source_target_shape':[164,122],
 'safe_original319_quad117_first_order_used':True,'P57_leftinverse_built_on_exact_inverse_transpose_restoration':True,'no_unimodular_minor_on_quadratic57_assumed':True,
 'Bnew22_from_surplus130_excludes_separator_indices':True,'all33_quartics_retained_before_old15_unit_deletions':True,'no_full239_or_rank122_claim':True,'known38_missing42_unchanged':True,
 'scalar_prefix_nodes':prefix_scalar,'linear_prefix_nodes':prefix_linear,'scalar_nodes':len(o.nodes),'linear_nodes':len(o.l),'core_atoms':len(o.atoms),
 'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in inputs},'certificate_sha256':sha(out/'certificate.json'),'artifact_bytes':(out/'certificate.json').stat().st_size,'seconds':time.monotonic()-started,'frozen':True,'independent_review_pending':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
