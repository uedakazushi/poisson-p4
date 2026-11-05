#!/usr/bin/env python3
"""Exact known planeG1 and12-source-adjoint input; no Φ4/M2 values."""
import argparse,json,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,srec,mv
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_GRAPH_PLANE34_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 graph=read('reports/AMBIENT_B18_ELL2_T2_LINEAR_GRAPH4_SOURCE_20261105_v63/actual23-linear-graph-four-common-contractions-and-six-adjoint-input.json')
 neg=read('reports/AMBIENT_B18_ELL2_T2_GRAPH_INVARIANT_ADJOINT6_20261105_v63/actual-two-graph-invariant-six-adjoint-certificates.json')
 read('reports/AMBIENT_B18_ELL2_T2_GRAPH_INVARIANT_ADJOINT6_20261105_v63_REVIEW/report.json')
 t=read(H.PATHS['transport']);pref=read(H.PATHS['active']);res=read(H.PATHS['resolution']);physical=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json')
 A=[list(map(H.parse,v)) for v in physical['C2_coordinate_actions9_sparse_columns']]
 lifts=[]
 for col in res['physical_lift_m0_columns']:lifts.append({(i,tuple(m)):F(c) for i,p in col for m,c in p if F(c)})
 def shifted(p,m):return {(i,(q[0]+m[0],q[1]+m[1])):c for (i,q),c in p.items() if q[1]+m[1]<2}
 active=pref['active_variables'];plane=[]
 for rawj in t['active_F1_free17']:
  block,j=divmod(rawj,14);var=8 if block==0 else active[block-1];p=shifted(lifts[j],(1,0) if var==7 else (0,1)) if var in [7,8] else {}
  for i,c in A[var][j].items():p=add(p,lifts[i],-c)
  assert all(2<=sum(m)<=3 for i,m in p),'planeG1 hasonlydegrees2/3';plane.append(p)
 plane += [{} for _ in range(9)];assert len(plane)==26
 degrees=[sorted({sum(m) for i,m in p}) for p in plane]
 # G1linear source12 is zero, so e115 is an obstruction for ALL132quadratic outputs.
 G1lin=list(map(H.parse,graph['actual_G1linear3_by26_columns']))
 for r in neg['results']:
  d=H.parse(r['restricted_support_obstruction234']);q={}
  for k,c in d.items():
   source,var=divmod(k,9)
   assert not G1lin[source]
  assert r['nonzero_value']!='0'
 # Domain F1[26]*(r,ell) +F2[26]*(r²,r*ell), all104 independent nuisance slots.
 # Output Φ3 andΦ4, eachphysical3*(r^degree,r^(degree-1)*ell):12 coordinates.
 columns=[{} for _ in range(12)]
 for source,p in enumerate(plane):
  for (target,m),c in p.items():
   degree=sum(m)
   for v,mon in enumerate([(1,0),(0,1)]):
    product=(m[0]+mon[0],m[1]+mon[1])
    if product[1]<2:
     output=(degree-2)*6+2*target+product[1];columns[output]=add(columns[output],{2*source+v:c})
   if degree==2:
    for v,mon in enumerate([(2,0),(1,1)]):
     product=(m[0]+mon[0],m[1]+mon[1])
     if product[1]<2:columns[6+2*target+product[1]]=add(columns[6+2*target+product[1]],{52+2*source+v:c})
 psi=[]
 for functional in graph['two_invariant_functionals_on_all234_F1_coefficients']:
  q={}
  for k,c in H.parse(functional).items():
   source,var=divmod(k,9);assert var in [7,8];q[2*source+(var-7)]=c
  psi.append(q)
 primitive='reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63/phase2-all105-equations-and-primitive56-source-columns.json';hs[primitive]=H.sha(ROOT/primitive)
 slots=sum(map(len,plane+columns+psi))+4*12*12+2*104+10000;assert slots<20000
 data={'exact_plane_G1_polynomial_columns':[[[i,list(m),str(c)] for (i,m),c in sorted(p.items())] for p in plane],
  'plane_degrees_per_source_column':degrees,'source_formula':'G1new_active_free17=p0*d1raw(Efree17), sincep0*d1(a11)=0 andd1*d2=0. p0 sendsraw14basistoactualphysicalmonomiallifts; noformalunitinverseTaylorquery.',
  'all_G1_plane_linear_coefficients_zero':True,'all_source12_to16_quadratic_V2_invariant_obstructions_extend_to_full132':True,
  'exact_equations':'F=Mnew*H0v,F0=0. Φ3=G1plane2*F1; Φ4=G1plane3*F1+G1plane2*F2. NoM2valueallowed: rho4mustannihilateALL52F2nuisancecoordinates.',
  'domain104_order':'F1 source26*(r,ell), thenF2 source26*(r²,r*ell)',
  'output12_order':'Φ3physical3*(r³,r²ell),thenΦ4physical3*(r⁴,r³ell)',
  'twelve_source_adjoint_columns104':[srec(v) for v in columns],
  'two_invariant_functionals_on104_with_F2_zero':[srec(v) for v in psi],
  'cached_Phi3_source':'primitive V56 sourcecombo Wv=(23:-3375/11,29:-50,30:25/2); B56containsgenuinecommon105bivariateplane3. NoCBzero7copyreplacement.',
  'Phi4_source_callbacks':'sameprimitive3combo ordinaryr⁴/r³ell, NOTzero-filledfromB56. Needadjointfirst, thenexplicitoverlapwithold9physicalφ4r⁴/functionals orraw33functional-factorproof beforeanynewsourcequery.',
  'all_Phi4_M1_M2_values_not_computed':True,'twelve_adjoint_solve_not_started':True,
  'whole_twelve_adjoint_slot_bound':slots,'limits':{'seconds':60,'memory_bytes':2_000_000_000,'raw_artifact_bytes':32*1024*1024,'QQ_coefficients':150000,'bits':384}}
 ep=out/'actual-plane-G1-degree23-and-twelve-adjoint-source-input.json';H.dump(ep,data)
 report={'id':STEM,'all_passed':True,'exact_known_G1_plane_polynomial':True,'plane_degree_support':sorted({d for ds in degrees for d in ds}),
  'twelve_adjoint_solve_started':False,'new_physical_or_M1_M2_Phi4_values':False,'whole_twelve_adjoint_slot_bound':slots,
  'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','plane_degree_support','whole_twelve_adjoint_slot_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
