#!/usr/bin/env python3
"""Finite source transport to free281, after the compact relation's peer acceptance."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_FREE281_CIRCUIT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canon(ts):
 d={}
 for ex,c in ts:
  m=tuple(ex);z=d.get(m,Q())+Q(c)
  if z:d[m]=z
  else:d.pop(m,None)
 return [[list(m),str(c)]for m,c in sorted(d.items(),reverse=True)]
class DAG:
 def __init__(self):
  self.atoms=[];self.atom_ids={};self.nodes=[['z']];self.node_ids={('z',):0};self.one=self.const([[[0,0,0],'1']],[[[0,0,0],'1']])
 def node(self,n):
  t=tuple(n)
  if t not in self.node_ids:self.node_ids[t]=len(self.nodes);self.nodes.append(n)
  return self.node_ids[t]
 def const(self,n,d):
  n,d=canon(n),canon(d);assert d
  if not n:return 0
  if n==d:n,d=[[[0,0,0],'1']],[[[0,0,0],'1']]
  if len(n)==len(d)==1 and n[0][0]==d[0][0]==[0,0,0]:n,d=[[[0,0,0],str(Q(n[0][1])/Q(d[0][1]))]],[[[0,0,0],'1']]
  key=json.dumps([n,d],separators=(',',':'))
  if key not in self.atom_ids:self.atom_ids[key]=len(self.atoms);self.atoms.append([n,d])
  return self.node(['c',self.atom_ids[key]])
 def add(self,a,b):
  if not a:return b
  if not b:return a
  return self.node(['a',a,b])
 def sub(self,a,b):
  if not b:return a
  if a==b:return 0
  return self.node(['s',a,b])
 def mul(self,a,b):
  if not a or not b:return 0
  if a==self.one:return b
  if b==self.one:return a
  return self.node(['m',a,b])
 def summands(self,terms):
  z=0
  for c,n in terms:z=self.add(z,self.mul(c,n))
  return z
def main(out,reviewpath):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 compact=ROOT/'reports/S61_NONSQUARE_ELL_COMPACT_RELATION_20261104_v63';cr=read(compact/'report.json');review=read(reviewpath)
 assert cr['all_literal_source_identity_and_outside37_gate_checks_passed']and review['all_passed']
 for name,h in cr['artifact_sha256'].items():assert sha(compact/name)==h
 # Reuse the accepted source304 and full P22. The old arbitrary B22 contains BOTH separator coordinates.
 pbase=ROOT/'reports/S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63';oldbasis=read(pbase/'global282-basis.json');atlaspath=ROOT/oldbasis['source304_atlas_path'];assert sha(atlaspath)==oldbasis['source304_atlas_sha256'];atlas=read(atlaspath)
 Ppath=pbase/'global-P22.json';assert sha(Ppath)==oldbasis['P22_sha256'];Pd=read(Ppath)
 Kpath=ROOT/'reports/S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63/kernel22.json';Kd=read(Kpath)
 PReview=ROOT/'reports/S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63_REVIEW/report.json';assert read(PReview)['all_passed']
 full=read(compact/'certificate.json');defs=full['source319_definitions'];keep=atlas['retained304_source319_indices'];assert len(keep)==304 and keep[:174]==list(range(174))
 assert Kd['actual_source174_definitions']==defs[:174]and Pd['source_variables']==defs[:174]
 separator=[182,186];assert [defs[k]for k in separator]==[['cubic_normalized',43,20,1],['cubic_normalized',43,20,6]]
 oldB=oldbasis['old_source319_indices_arbitrary_B22'];assert all(k in oldB for k in separator)
 surplus=keep[174:];B=[k for k in surplus if k not in separator][:22];C=[k for k in surplus if k not in B];assert len(C)==108 and all(k in C for k in separator)
 p,q=[174+C.index(k)for k in separator];assert p!=q
 o=DAG();one=o.one;R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field()
 def atomstr(v):
  c=K.from_expr(s.sympify(v));enc=lambda p:[[list(m),str(c)]for m,c in p.items()]
  return o.const(enc(c.numer),enc(c.denom))
 # The actual new319 vector is a polynomial input; no coefficient system is solved.
 ls=(compact/'full319-kernel.txt').read_text().splitlines();assert ls[0]=='319 1';new=[]
 for line in ls[1:]:
  t=line.split();n=int(t[0]);assert len(t)==1+4*n
  new.append(o.const([[list(map(int,t[i+1:i+4])),t[i]]for i in range(1,len(t),4)],[[[0,0,0],'1']]))
 assert len(new)==319
 # Import only weight-zero scalar roots actually needed for residual4's quadratic cofactors.
 bbase=ROOT/'reports/S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63';bp=bbase/'certificate.json';bd=read(bp);arith=bd['arithmetic'];cache={0:0}
 def bnode(i):
  if i in cache:return cache[i]
  assert arith['node_weights'][i]==0
  nd=arith['nodes'][i]
  if nd[0]=='c':v=o.const(*arith['atoms'][nd[1]])
  elif nd[0]=='m':v=o.mul(bnode(nd[1]),bnode(nd[2]))
  else:
   assert nd[0]=='l';v=o.summands((o.const(*arith['atoms'][c]),bnode(j))for c,j in nd[1])
  cache[i]=v;return v
 rbase=ROOT/'reports/S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63';rd=read(rbase/'certificate.json');rows15=[r['after35_deletions_coefficients319']for r in bd['raw15_source_lifts']];U=rd['residual_transform_rows15to4'];T=rd['pivot_row_transform'];coeff=[]
 for j in range(15):coeff.append(o.summands((atomstr(T[3][i]),atomstr(U[i][j]))for i in range(4)))
 rlast=[];ci={g:j for j,g in enumerate(rd['cubic_source319_indices'])}
 for g in range(319):
  if g<174:r=o.summands((coeff[j],bnode(rows15[j][g]))for j in range(15)if coeff[j])
  elif g in ci:r=atomstr(rd['post_pivot_cubic_matrix'][3][ci[g]])
  else:r=0 # exact all33 quartic zero is part of the accepted residual4 certificate
  rlast.append(r)
 piv=[x['source319_column']for x in rd['unit_pivot_steps']];assert piv==[176,177,237,181]
 assert all(rlast[piv[i]]==(one if i==3 else 0)for i in range(4))
 assert all(new[k]==0 for k in atlas['deleted15_source319_indices_in_pivot_order']if k!=181)
 v319=[o.sub(new[i],o.mul(new[181],rlast[i]))for i in range(319)]
 assert all(v319[k]==0 for k in atlas['deleted15_source319_indices_in_pivot_order'])
 x=v319[:174];u=[v319[k]for k in B];h=[v319[k]for k in C]
 # Apply the accepted P22 linear operator to the NEW quadratic cofactor vector.
 pc=Pd['circuit'];sc={0:0};lc={0:0}
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
  if nd[0]=='v':v=x[nd[1]]
  else:assert nd[0]=='l';v=o.summands((ps(c),pl(j))for c,j in nd[1])
  lc[i]=v;return v
 Px=[pl(i)for i in Pd['P22_linear_output_roots_in_frozen_original_kernel_order']];assert len(Px)==22
 k=[o.sub(Px[j],u[j])for j in range(22)];Krows=[[]for _ in range(174)]
 for j,col in enumerate(Kd['sparse_columns']):
  for i,n,de in col:Krows[i].append((o.const(n,de),k[j]))
 z=[o.sub(x[i],o.summands(Krows[i]))for i in range(174)];v282=z+h;assert len(v282)==282
 theta=-K(1024000)/K(243)*K(a)*K(100*a+18*e**2-9*e-2)**3*K(5*a*d-6*e-4)
 enc=lambda p:[[list(m),str(c)]for m,c in p.items()];thetaid=o.const(enc(theta.numer),enc(theta.denom));inverseid=o.const(enc(theta.denom),enc(theta.numer))
 sumcheck=o.sub(o.add(v282[p],v282[q]),thetaid)
 # This ONE actual scalar identity does not depend on the large P22(z) roots:
 # both coordinates are in the literal Cnew block. Evaluate it exactly over QQ(a,e,d).
 valcache={0:K.zero}
 def value(i):
  if i in valcache:return valcache[i]
  nd=o.nodes[i]
  if nd[0]=='c':
   num,den=o.atoms[nd[1]];load=lambda ts:R.from_dict({tuple(m):s.QQ(c)for m,c in ts});v=K(load(num))/K(load(den))
  elif nd[0]=='a':v=value(nd[1])+value(nd[2])
  elif nd[0]=='s':v=value(nd[1])-value(nd[2])
  else:assert nd[0]=='m';v=value(nd[1])*value(nd[2])
  valcache[i]=v;return v
 assert value(sumcheck)==K.zero
 w=v282[:];w[p]=thetaid # exact alias uses the verified unit identity, not a heuristic simplification
 normalized=[o.mul(v,inverseid)for v in w];normalized[p]=one
 # Constant source-coordinate change w_p=v_p+v_q,w_q=v_q; the other coordinates agree.
 # Its generator change is M'_p=M_p, M'_q=M_q-M_p. Drop the p generator using its known unit kernel entry.
 source281=[i for i in range(282)if i!=p]
 evidence={'source319_definitions':defs,'source304_retained319_indices':keep,'source304_new_block_A174':list(range(174)),'source304_new_block_B22':B,'source304_new_block_C108':C,
 'old_B22_contains_both_separator_coordinates':True,'separator319_indices':separator,'separator_new282_indices':[p,q],
 'last_residual4_source319_roots':rlast,'new_source319_roots':new,'transported_source319_roots':v319,
 'source319_to304_rule':'v=new-K181*lastResidual4; all15 deleted coordinates0',
 'new_quadratic174_roots':x,'new_B22_roots':u,'new_C108_roots':h,'P22_on_x22_roots':Px,'kernel_coordinate22_roots':k,
 'new_kernel282_roots':v282,'source304_to282_rule':'k=P22*x-u; z=x-K22*k; h=h; exact inverse follows P22*K22=I22',
 'known_unit_separator_root':thetaid,'known_unit_separator':'-1024000*a*C^3*J/243','separator_identity_zero_root':sumcheck,'separator_identity_exact_QQ_checked':True,'separator_evaluated_DAG_nodes':len(valcache),
 'separator_identity_proof':'new182+new186=theta by compact native; lastResidual4[182]+lastResidual4[186]=0 by exact residual4; both are literal Cnew coordinates',
 'constant_coordinate_change':{'forward':'w_p=v_p+v_q; w_q=v_q','inverse':'v_p=w_p-w_q; v_q=w_q','determinant':1},
 'new_kernel_after_constant_change_roots':w,'normalized_kernel_with_pivot_one_roots':normalized,
 'actual_free281_source_indices_in_changed282':source281,'actual_free281_shape':[281,239],
 'generator_recipe':'Let M282new=[A+Bnew*P22,Cnew]. Apply Mprime_p=M282new_p,Mprime_q=M282new_q-M282new_p, other generators unchanged; delete p. No281x239 coefficient expansion.',
 'new_source281_to319_lift':'For y with p omitted, put w_p=0,w_i=y_i; v_p=-w_q,v_q=w_q,other v_i=w_i. Split v=(z174,h108); fill original304 as (z,P22*z,h), then put deleted319 coordinates0.',
 'unit_deletion_proof':'The exact kernel w has w_p=theta unit. Hence Mprime_p= -sum_{i!=p}(w_i/theta)*Mprime_i. The normalized kernel and e_i for i!=p give a full source basis with determinant sign. Therefore the retained281 generators have exactly the same image/cokernel.',
 'no_full239_or_all_kernel_claim':True,'no_new_solver_pivot_rank_GB_point_field_or_expanded_quartic':True,
 'circuit':{'atoms_num_den':o.atoms,'nodes':o.nodes}}
 (out/'certificate.json').write_text(json.dumps(evidence,separators=(',',':'))+'\n')
 assert len(o.nodes)<200000 and (out/'certificate.json').stat().st_size<32*1024*1024
 inputs=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),reviewpath,compact/'report.json',compact/'certificate.json',compact/'full319-kernel.txt',pbase/'global282-basis.json',Ppath,PReview,atlaspath,Kpath,bp,bbase.parent/(bbase.name+'_REVIEW')/'report.json',rbase/'certificate.json',rbase.parent/(rbase.name+'_REVIEW')/'report.json']
 report={'id':ID,'finite_exact_source_transport_recipe_constructed':True,'new_actual_free_source_rank':281,'target_rank':239,'same_actual319_cokernel':True,'full239_still_open':True,'missing_uniform_directions_for80':42,
 'old_B22_separator_mixing_checked_and_fixed':True,'new_B22_excludes_both_separator319_indices':True,'separator_new282_indices':[p,q],
 'no_new_solver_rank_or_coefficient_expansion':True,'known_unit_factors_only':['a','d','F','B','C','J'],
 'scalar_circuit_nodes':len(o.nodes),'core_atoms':len(o.atoms),'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in inputs},
 'certificate_sha256':sha(out/'certificate.json'),'seconds':time.monotonic()-began,'artifact_bytes':(out/'certificate.json').stat().st_size,'frozen':True,'independent_review_pending':True}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);p.add_argument('--compact-review',type=Path,required=True);a=p.parse_args();main(a.output.resolve(),a.compact_review.resolve())
