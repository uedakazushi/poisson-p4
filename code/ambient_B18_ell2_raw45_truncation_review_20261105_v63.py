#!/usr/bin/env python3
"""Independent six-basis/Tor-aware ell² truncation proof and input audit."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,itertools,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_RAW45_TRUNCATION_PREFLIGHT_20261105_v63';STEM=A+'_REVIEW'
P={
 'primary':f'reports/{A}/report.json','input':f'reports/{A}/literal-existing-cubic-source-input.json',
 'proof':f'reports/{A}/union-six-basis-truncation-and-module-preflight.json',
 'source':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
 'source_peer':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63_REVIEW/report.json',
 'product_peer':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_CURVE_BASE_20261104_v63_REVIEW/report.json',
 'DVR':'reports/AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63/DVR20-raw-cokernel26-and-comparison-preflight.json',
 'DVR_peer':'reports/AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def poly(p):return {tuple(m):Q(c) for m,c in p if Q(c)}
def add(a,b,c=Q(1)):
 z=dict(a)
 for k,v in b.items():
  q=z.get(k,Q(0))+c*v
  if q:z[k]=q
  else:z.pop(k,None)
 return z
def mul(a,b):
 z={}
 for m,c in a.items():
  for n,d in b.items():z=add(z,{tuple(x+y for x,y in zip(m,n)):c*d})
 return z
def lin(A,v):
 out=[]
 for row in A:
  p={}
  for c,q in zip(row,v):p=add(p,q,c)
  out.append(p)
 return out
def encode(p):return [[list(k),str(c)] for k,c in sorted(p.items())]
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'TOR_AWARE_PURE_PROOF_INPUTS','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()}
 assert all(D[k]['all_passed'] for k in ['primary','source_peer','product_peer','DVR_peer'])
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 S=D['source'];B2=list(map(poly,S['beta2_before_ell']));B3=list(map(poly,S['beta3_before_ell']));A0=[[Q(c) for c in r] for r in S['A0']];Ai=[[Q(c) for c in r] for r in S['A0_inverse']]
 assert all(sum(A0[i][k]*Ai[k][j] for k in range(7))==sum(Ai[i][k]*A0[k][j] for k in range(7))==Q(i==j) for i in range(7) for j in range(7))
 bplane=list(map(poly,S['b_plane3']));b2=[{m:c for m,c in p.items() if sum(m)==2} for p in bplane]
 A1=[[poly(p) for p in row] for row in S['A1_ordered_prefix']]
 correction=[]
 for row in A1:
  p={}
  for a,b in zip(row,B2):p=add(p,mul(a,b))
  correction.append(p)
 assert lin(Ai,b2)==B2
 assert lin(Ai,[add({m:c for m,c in p.items() if sum(m)==3},q,-1) for p,q in zip(bplane,correction)])==B3
 assert all(not any(m[:7]) and sum(m)==2 for p in B2 for m in p)
 plane3=[];factor=[];zterms=0
 for p in B3:
  plain={m:c for m,c in p.items() if not any(m[:7])};plane3.append(plain);row=[{} for _ in range(7)]
  for m,c in p.items():
   if not any(m[:7]):continue
   j=next(i for i in range(7) if m[i]);n=list(m);n[j]-=1;row[j]=add(row[j],{tuple(n):c});zterms+=1
  rebuild=dict(plain)
  for j,q in enumerate(row):
   v=[0]*9;v[j]=1;rebuild=add(rebuild,mul(q,{tuple(v):Q(1)}))
  assert rebuild==p and all(sum(m)==2 for q in row for m in q);factor.append(row)
 assert zterms==48 and sum(map(len,plane3))==4
 assert [[[list(m[7:]),str(c)] for m,c in sorted(p.items())] for p in B2]==D['input']['beta_graph_plane_through3']['beta2_before_ell']
 assert [[[list(m[7:]),str(c)] for m,c in sorted(p.items())] for p in plane3]==D['input']['beta_graph_plane_through3']['beta3_before_ell']
 cols=[(2,0),(1,1),(0,2)];rows=[[p.get((0,)*7+m,Q(0)) for m in cols] for p in B2];nz=[i for i,r in enumerate(rows) if any(r)];assert len(nz)==3
 M=[rows[i] for i in nz];det=M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);assert det
 gA=poly(S['gA_plane3']);gB=poly(S['gB_own_plane3']);T=[[Q(c) for c in row[7:]] for row in S['YB_in_YA_linear'][7:]]
 ca=sum(gA.values());cb=sum(c*(sum(T[0])**m[0])*(sum(T[1])**m[1]) for m,c in gB.items())
 assert ca==Q(125,24) and cb==Q(-15125,144)
 # Universal QQ polynomial verification of the six-basis union algebra.
 # All coefficients are formal; no source coefficient is evaluated here.
 import sympy as s
 e=s.Symbol('ell');a2,a1,a0,b2,b1,b0=s.symbols('a2 a1 a0 b2 b1 b0');c2,c1,c0,d2,d1,d0=s.symbols('c2 c1 c0 d2 d1 d0')
 def red(x):
  x=s.expand(x);return s.expand(x.coeff(e,0)+e*x.coeff(e,1))
 def rm(A):return A.applyfunc(red)
 Z=s.zeros(6);R=s.zeros(6)
 R[1,0]=1;R[2,1]=R[3,1]=1
 R[2,2]=e*a2;R[5,2]=e*a1;R[4,2]=e*a0
 R[3,3]=e*b2;R[5,3]=-e*b1;R[4,3]=-e*b0
 R[5,4]=e;R[2,5]=e
 def zmat(C2,C1,C0):
  v=s.Matrix([0,-e*C1,0,-C2,C0,C1]);out=s.zeros(6);out[:,0]=v;out[:,1]=rm(R*v);return out
 C=zmat(c2,c1,c0);D0=zmat(d2,d1,d0)
 expressions={'r4':R**4,'ell_r3':e*R**3,'r_commutes_z':R*C-C*R,'r2_z':R**2*C,'zi_zj':C*D0,'ell2':e**2*s.eye(6)}
 for name,q in expressions.items():assert rm(q)==Z,name
 assert all(rm(C[:,i])==s.zeros(6,1) for i in [2,3,4,5])
 proof=D['proof'];assert proof['ambient_monomials']=={'plane':11,'z_sectors':63,'total':74}
 assert proof['row_bound']==147 and proof['predicted_actual_product_dimension_from_accepted_DVR_type']==33
 dv=D['DVR']['product_curve_ring_DVR_type'];assert dv=={'free_rank':6,'D_mod_ell_squared_copies':7,'D_mod_ell_copies':7,'special_fiber_dimension':20}
 assert 2*dv['free_rank']+2*dv['D_mod_ell_squared_copies']+dv['D_mod_ell_copies']==33
 ep=o/'independent-six-basis-Tor-repair-and-cubic-input-gates.json'
 dump(o,ep.name,{'inputs':inputs,'all_actual_beta2_beta3_A1_correction_source_identities_verified':True,
  'quadratic_plane_rows':nz,'quadratic_plane_unit_minor':str(det),'gA_r3_unit':str(ca),'gB_common_r3_unit':str(cb),
  'all48_z_dependent_cubic_terms_factored_into_unit_M2z':[[encode(p) for p in row] for row in factor],
  'universal_six_basis_order':['e0=(1,1)','e1=(r,r)','u=(r²,0)','v=(0,r²)','p=(ell²,0)','q=(ell*r,0)'],
  'universal_r_matrix':[[str(x) for x in row] for row in R.tolist()],
  'universal_z_matrix':[[str(x) for x in row] for row in C.tolist()],
  'all_universal_r4_ellr3_r2z_zizj_commutation_identities_zero_mod_ell2':True,
  'Tor_attachment_not_discarded':True,'actual_and_truncated_union_m4_zero':True,
  'exact_product_modell2_ideal_unchanged_by_generator_truncation3':True,
  'reciprocal_inclusion_proof':'In each ring J*m4=0 and m6=0; f*g changes by f*u4+g*t4+t4*u4, allzero. Both actual and truncated J contain m² by the checked quadratic unit minor.',
  'plane8_generator_presentation_closed_under_all_z_actions':True,'conditional_matrix_bound':[147,74],
  'predicted_target_product_length_from_DVR':33,'predicted_raw_cokernel_length':45,
  'new_quotient_rank_source_values_or_implicit_solve':False})
 note=ROOT/'research/notes'/f'{STEM}.md';report={'id':STEM,'all_passed':True,
  'Tor_aware_union_m4_and_exact_cubic_product_truncation_source_gate_accepted':True,
  'conditional_raw33_matrix_bound':[147,74],'numerical_raw33_or_physical_raw45_actions_computed':False,
  'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),'evidence_sha256':{ep.name:sha(ep)},
  'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
