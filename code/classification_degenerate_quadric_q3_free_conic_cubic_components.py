#!/usr/bin/env python3
"""A single full150/actualR6 witness for each complete free conic+cubic source."""
import hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_three_hyperplanes_quadric_independent import ambient_matrices
from classification_rank3_mixed_component_review import minor_certificate
resource.setrlimit(resource.RLIMIT_CPU,(85,86));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
start=time.monotonic();path=Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_FREE_CONIC_CUBIC_SOURCES_20261107.json');src=json.loads(path.read_text());xx=s.symbols('a b c r z');mm=s.symbols('m0:6');witness=dict(zip(mm,[1,2,3,5,7,11]));dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(s.QQ);records=[]
def enc(P):return {str(ij):str(s.expand(P[ij])) for ij in itertools.combinations(range(5),2) if P[ij]}
for case in src['cases']:
 if case['quadratic_syzygy_kernel_dimension']!=2:continue
 P=s.zeros(5)
 for ij in itertools.combinations(range(5),2):P[ij]=s.sympify(case['full_unimodular_brackets'][''.join(map(str,ij))]);P[ij[::-1]]=-P[ij]
 scale=s.ilcm(1,*[s.denom(a) for ij in itertools.combinations(range(5),2) for a in s.Poly(P[ij],*xx,*mm).coeffs()]);P*=scale;point=P.subs(witness);dirs=[P.diff(m) for m in mm]
 D,V,M,F=ambient_matrices(point,xx,dirs);mats={'full_D':D,'vertical':V,'group_vertical_scale':M,'all_six_source':F};ranks={k:len(dm(A).rref()[1]) for k,A in mats.items()};assert D*F==s.zeros(350,F.cols)
 rec={'case':case['case'],'source_scale':str(scale),'source_parameters':list(map(str,mm)),'witness':{str(k):str(v) for k,v in witness.items()},'point_brackets':enc(point),'all_six_parameter_directions':[enc(R) for R in dirs],'ranks_QQ':ranks,'certificates':{k:minor_certificate(A,ranks[k]) for k,A in mats.items()},'matrix_sha256':{k:hashlib.sha256(json.dumps(A.tolist(),default=str).encode()).hexdigest() for k,A in mats.items()},'projective_family_dimension_lower':ranks['all_six_source']-25,'projective_tangent_dimension_upper':125-ranks['full_D'],'true_normal_gap':150-ranks['full_D']-ranks['all_six_source'],'all_actual_source_directions_in_kernel':True}
 records.append(rec);print(json.dumps({k:rec[k] for k in ['case','ranks_QQ','projective_family_dimension_lower','projective_tangent_dimension_upper','true_normal_gap']}),flush=True)
report={'task_id':'E_reduced_nonnormal','all_passed':True,'input_sha256':{str(path):hashlib.sha256(path.read_bytes()).hexdigest()},'cases':records,'producer_source':Path(__file__).read_text(),'seconds':time.monotonic()-start,'scope':'One full ambient rank and all six true source directions for each irreducible complete source; only equal upper/lower dimensions prove componenthood. No known-component membership follows from these ranks alone.'}
out=Path('reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_FREE_CONIC_CUBIC_COMPONENTS_20261107.json')
with out.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
