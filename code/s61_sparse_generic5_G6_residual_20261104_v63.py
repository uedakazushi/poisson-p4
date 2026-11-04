#!/usr/bin/env python3
"""ONE small actualG6 pull,125cofactor identity,and factored130th residual."""
import argparse,gzip,itertools,json,resource,signal,time
from collections import defaultdict
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from q_new_family_v_zero_degree26_support_preflight_20261104_v63 import read,sha
from s61_sparse_generic5_input_preflight_20261104_v63 import add,mul,scale
from s61_sparse_generic5_full134_normal_circuit_20261104_v63 import Pool
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_G6_RESIDUAL_20261104_v63';ZERO=(0,0,0,0,0)
def addrow(a,m,p):
 a[m]=add(a.get(m,{}),p)
 if not a[m]:del a[m]
def enc(p):return[[str(c),*e]for e,c in sorted(p.items())]
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));signal.alarm(60);out=out.resolve();out.mkdir(parents=True,exist_ok=False)
 pref=ROOT/'reports/S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63';src=ROOT/'reports/S61_SPARSE_GENERIC5_FULL134_NORMAL_CIRCUIT_20261104_v63';res=ROOT/'reports/S61_SPARSE_GENERIC5_RESIDUAL129_CIRCUIT_20261104_v63';plane=ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_RANK3_CIRCUIT_20261104_v63';gp=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json';gr=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_REVIEW_20261031.json'
 assert read(res/'report.json')['all_passed'];g=read(gp);review=read(gr);assert g['all_passed']and review['all_passed']and g['sample_augmented_quadratic_span_rank']==135 and review['sample_augmented_quadratic_rank']==135
 assert g['linear_syzygy_kernel_dimension']==36 and g['canonical_Schouten_Bianchi_dimension']==35 and review['canonical_Bianchi_rank']==35
 c=read(pref/'certificate.json');V=s.symbols('H L s0 q0 t0');beta={i:{tuple(e):Q(v)for e,v in s.Poly(s.sympify(f),*V).terms()}for i,f in zip(c['historical_weight1_indices'],c['source_u16'])};beta={i:{e:v for e,v in p.items()if v}for i,p in beta.items()};qids=c['historical_weight2_indices'];qi={j:i for i,j in enumerate(qids)};assert qids==g['q_to_u_indices_zero_based'];data=read(pref/'pulled134-input.json.gz');F=[]
 for eq in data['equations']:
  row={}
  for m,e,v in eq:row.setdefault(tuple(m),{})[tuple(e)]=Q(v)
  F.append(row)
 target={};target_upper=set();alive=0
 for mon,n,d in g['universal_new_equation']:
  p={ZERO:Q(n,d)};qm=[]
  for i in mon:
   if i in qi:qm.append(qi[i])
   else:p=mul(p,beta[i])
  m=tuple(sorted(qm));assert len(m)<=2
  if p:alive+=1
  target_upper.update((m,e)for e in p);addrow(target,m,p)
 rhs_upper=set();contributions=0;cof=[]
 for index,value in g['selected_non_Bianchi_syzygy']:
  q,i=divmod(index,70);cof.append([i,q,str(value)])
  for m,p in F[i].items():rhs_upper.update((tuple(sorted((*m,q))),e)for e in p);contributions+=len(p)
 totalupper=len(target_upper)+len(rhs_upper);assert totalupper<=150000
 pre=dict(G6_global_terms=9257,surviving_global_terms=alive,G6_support_upper=len(target_upper),rhs_support_upper=len(rhs_upper),TOTAL_new_explicit_support_upper=totalupper,hard150k_total_guard=True,rhs_term_contributions=contributions,rhs_q3_positions=sum(len(m)==3 for m,e in rhs_upper),no126inverse_GB_point_prime_G28expansion=True)
 (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n')
 rhs={};peak=0
 for i,q,value in cof:
  for m,p in F[i].items():addrow(rhs,tuple(sorted((*m,q))),scale(p,Q(value)))
  peak=max(peak,sum(map(len,rhs.values())));assert peak+sum(map(len,target.values()))<=150000
 assert rhs==target and all(len(m)<=2 for m in rhs)
 pull=dict(parameter_variables=c['parameter_variables'],G6_q45_coefficients=[[list(m),enc(p)]for m,p in sorted(target.items())],linear_q_cofactors125=cof,identity='G6=sum_i<70 K_i(q)*F_i, allq3coefficientscancelEXACTLY',all_q3_positions_zero=True,q_to_u_indices_zero_based=qids)
 (out/'pulled-G6-and-cofactors.json.gz').write_bytes(gzip.compress(json.dumps(pull,separators=(',',':')).encode(),mtime=0))
 direct=read(src/'full134-direct-circuit.json.gz');blocks=read(src/'coefficient-block-circuit.json.gz');R=read(res/'residual129-factored-circuit.json.gz');Wdata=read(plane/'plane-W-original.json.gz');oldrows=blocks['rows134'];oldnodes=direct['nodes']+blocks['appended_nodes'];Nrows=direct['normal_q_rows42'];ni={i:k for k,i in enumerate(Nrows)};J=R['selected_original134_rows'];non=R['residual_original134_rows'];pool=Pool();oldcache={};rescache={}
 def prior(i,degree=0):
  if(i,degree)not in oldcache:
   if oldnodes[i][0]=='qq':oldcache[i,degree]=pool.const(oldnodes[i][1])
   else:oldcache[i,degree]=pool.node('external_full134_node',[i],degree,None)
  return oldcache[i,degree]
 def rr(i,degree=0):
  if(i,degree)not in rescache:rescache[i,degree]=pool.node('external_residual129_node',[i],degree,None)
  return rescache[i,degree]
 def neg(x):return pool.prod(pool.const(-1),x)
 E=[[prior(i)for i in row]for row in direct['E45x3']];qf=[prior(i,1)for i in direct['forward_q45']];G=pool.summ([pool.prod(pool.atom(p),*[qf[i]for i in m])for m,p in target.items()])
 C=[{tuple(t[1:]):Q(p['content'])*int(t[0])for t in p['terms']}for p in Wdata['center45']];known=read(pref/'known-five-input.json.gz');d=[{tuple(e):Q(v)for e,v in p}for p in known['Q52_minus_center']];Ppoly={tuple(t[1:]):Q(t[0])for t in Wdata['common_P']};P=pool.atom(Ppoly);W=[[pool.node('external_original_W',[i,j],0,True)for j in range(2)]for i in range(45)]
 c0=dict(target.get((),{}));L=[dict(target.get((i,),{}))for i in range(45)];Q2={m:p for m,p in target.items()if len(m)==2}
 for i in range(45):c0=add(c0,mul(target.get((i,),{}),C[i]))
 for (i,j),p in Q2.items():c0=add(c0,mul(p,mul(C[i],C[j])));L[i]=add(L[i],mul(p,C[j]));L[j]=add(L[j],mul(p,C[i]))
 SG=[pool.prod(P,pool.atom(c0))]
 for h in range(2):SG.append(pool.summ([pool.prod(pool.atom(L[i]),W[i][h])for i in range(45)]))
 for h in range(2):SG.append(pool.summ([pool.prod(pool.atom(p),pool.summ([pool.prod(W[i][h],pool.atom(d[j])),pool.prod(W[j][h],pool.atom(d[i]))]))for (i,j),p in Q2.items()]))
 lam=[pool.prod(rr(R['inverseDelta5_node']),pool.summ([pool.prod(SG[h],rr(R['adjugate5_nodes'][h][j]))for h in range(5)]))for j in range(5)];Fs=[prior(direct['F134'][j],2)for j in J];RG=pool.summ([G,*[neg(pool.prod(v,f))for v,f in zip(lam,Fs)]])
 lin=[pool.summ([pool.atom(L[i]),*[neg(pool.prod(v,prior(oldrows[j]['linear_n42'][k])))for v,j in zip(lam,J)]])for k,i in enumerate(Nrows)];crossparts=[[]for _ in range(126)];nn=[[]for _ in range(42)];Z=[pool.node('z',[i],1,True)for i in range(3)];N=[pool.node('n',[i],1,True)for i in range(42)]
 for (i,j),p in Q2.items():
  a=pool.atom(p)
  if i in ni:
   for h in range(3):crossparts[42*h+ni[i]].append(pool.prod(a,E[j][h]))
  if j in ni:
   for h in range(3):crossparts[42*h+ni[j]].append(pool.prod(a,E[i][h]))
  if i in ni and j in ni:
   a0,b0=sorted((ni[i],ni[j]));nn[a0].append(pool.prod(a,N[b0]))
 cross=[pool.summ([pool.summ(p),*[neg(pool.prod(v,prior(oldrows[j]['cross_z_n126'][k])))for v,j in zip(lam,J)]])for k,p in enumerate(crossparts)];qparts=list(map(pool.summ,nn));A=[pool.summ([lin[k],*[pool.prod(cross[42*h+k],Z[h])for h in range(3)],qparts[k],*[neg(pool.prod(v,rr(R['source_normal2_factor134x42'][j][k],1)))for v,j in zip(lam,J)]])for k in range(42)]
 K=[pool.z]*70
 for i,q,value in cof:K[i]=pool.summ([K[i],pool.prod(pool.const(Q(value)),qf[q])])
 # The original125-cofactor identity transported through the actualrowchange.
 termsR=[pool.prod(K[i],rr(R['residual_direct_F129'][non.index(i)],2))for i in range(70)if i not in J and K[i]!=pool.z];Bcof=[]
 for h,j in enumerate(J):Bcof.append(pool.summ([K[j],neg(lam[h]),*[pool.prod(K[i],rr(R['lambda129x5'][non.index(i)][h]))for i in range(70)if i not in J and K[i]!=pool.z]]))
 cofrhs=pool.summ([*termsR,*[pool.prod(v,f)for v,f in zip(Bcof,Fs)]]);assert pool.deg[RG]<=2 and all(pool.deg[x]<=1 for x in A)and pool.deg[cofrhs]<=3
 atoms=sum(map(len,pool.atoms));assert atoms+peak+sum(map(len,target.values()))<=150000
 obj=dict(nodes=pool.nodes,scalar_parameter_atoms=pool.atoms,source_files=dict(full134_direct=str((src/'full134-direct-circuit.json.gz').relative_to(ROOT)),full134_blocks=str((src/'coefficient-block-circuit.json.gz').relative_to(ROOT)),residual129=str((res/'residual129-factored-circuit.json.gz').relative_to(ROOT)),W=str((plane/'plane-W-original.json.gz').relative_to(ROOT))),external_full134_node_semantics='Concatenateoldnodes; samez/nleaves, referencedqqcoeffdegree0/qforward1/F2',external_residual129_node_semantics='Samez/nleaves and originalDAG; degree0selector/degree1Qfactor/degree2residual',G6_direct_node=G,G6_five_functional_column_nodes=SG,lambda_G6_5=lam,R_G6_direct_node=RG,R_G6_linear_n42=lin,R_G6_cross126=cross,R_G6_linear_normal_cofactor42=A,G6_K70_linear_q=K,transported_selected_F5_coefficients=Bcof,transported_cofactor_rhs_node=cofrhs,identity='R_G6=sum_i<70,unselected K_i(q)*R_i+sum_j Bcof_j(q)*Fselected_j; equalssum_k n_k*A_k',all45unknowns_retained=True,allnormal2_retained=True,only_parameter_inverses_from_prior=True,no_extra_affine_condition=True)
 raw=json.dumps(obj,separators=(',',':')).encode();assert len(raw)<=32*1024*1024;(out/'G6-normal-residual-circuit.json.gz').write_bytes(gzip.compress(raw,mtime=0))
 structural=dict(full135_covariance_independent_premise='S61_HOMOGENIZED_EQUIVARIANCE_20261101: syzygyS36/B35 trivialSp4line, fullBianchi givesbeta²F4+betaF5, beforeqhomogenization; qnondependenttransformationcoefficients',old_qhighest_rank135_nonzero=True,new_explicit135minor=False,sparse_source_rank135_by_full135_covariance_and_dominance=True,plane_restriction_rank5=True,cross_target_dimension126_and_RAW134_rank126=True,normal_only_quadratic_rowspace_dimension=135-5-126,explicit_basis_of_four_normal_only_quadrics_computed=False,canonical_leftinverse_recipe='ForX130x126 overformallyrealK, L=(X^tX)^(-1)X^t existsbyCauchyBinet/sumofsquares. ProjectorI130-XL hasrank4; itsproductwithpure-normalpartU(n) givesfour-dimensionalnormal-onlyquadraticspace.',Gram_leftinverse_computed=False,finite_residual_goal='FornormalidealJ ofpure-normalquadruple+crosscompatcubics+clearedselectedplanequartics, Bempty iff J:(n)^infinity=(1), equivalentlysome(n)^d subsetJ; no suchmembershipcertificatenow',not_exhaustion=True)
 (out/'fourth-normal-relation-gates.json').write_text(json.dumps(structural,indent=2)+'\n');assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<=32*1024*1024
 paths=[Path(__file__),ROOT/f'research/notes/{ID}.md',gp,gr,ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md',ROOT/'research/notes/S61_UNIVERSAL_WEIGHT6_REVIEW_20261031.md',pref/'report.json',pref/'certificate.json',pref/'pulled134-input.json.gz',pref/'known-five-input.json.gz',res/'report.json',res/'residual129-factored-circuit.json.gz',src/'report.json',src/'full134-direct-circuit.json.gz',src/'coefficient-block-circuit.json.gz',plane/'plane-W-original.json.gz',ROOT/'reports/S61_SPARSE_GENERIC5_CROSS_RANK_COVARIANCE_20261104_v63_REVIEW/report.json']
 report=dict(id=ID,all_passed=True,G6_pulled_terms=sum(map(len,target.values())),G6_qdegrees=sorted(set(map(len,target))),all125_cofactor_identity_QQ_exact=True,all_q3_cancelled=True,rhs_peak_terms=peak,new_scalar_atom_terms=atoms,new_TOTAL_support_bound=atoms+peak+sum(map(len,target.values())),nodes=len(pool.nodes),raw_circuit_bytes=len(raw),G6_plane_residual_and42normal_factors_complete=True,allnormal2_preserved=True,pure_normal_quadratic_kernel_dimension4_structural=True,no_Gram_inverse_basis_or_residualunit_computed=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,G6terms=report['G6_pulled_terms'],newTOTAL=report['new_TOTAL_support_bound'],nodes=len(pool.nodes),seconds=report['seconds'])))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
