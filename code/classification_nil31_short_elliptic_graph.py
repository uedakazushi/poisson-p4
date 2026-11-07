#!/usr/bin/env python3
"""Weakest full-graph recognition of the dense short2/3 source orbit."""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as s
from s61_model import S61,X,pushforward
from classification_rank3_twzero_graph_diagnostic import darboux
ROOT=Path(__file__).resolve().parents[1]
def main():
 started=time.monotonic();h,x,y,z,w=X;P=s.zeros(5);br=[h*z-h*y+2*z*w-y*y/3,h*z-2*y*z/3,h*h+y*w/3,-2*z*z/3,(h*y+z*w)/3,2*h*z/3]
 for (i,j),a in zip(itertools.combinations(range(1,5),2),br):P[i,j]=a;P[j,i]=-a
 # The one dense orbit anchor t=d=1, merely choosing a symplectic point on it.
 T=s.eye(5);T[2,0]=T[3,0]=-1;translated=pushforward(P,T);value=translated.subs(dict(zip(X,(1,0,0,0,0))))[1:,1:];assert value.det()!=0
 frame=s.eye(5);frame[1:,1:]=darboux(value);framed=pushforward(translated,frame);model=S61();normalized=model.normalize(framed);model.verify_point(normalized,check_equations=True);point=normalized['parameters'];us=s.symbols('u1:62');records={}
 paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/classification_rank3_twzero_graph_diagnostic.py',ROOT/'reports/CLASSIFICATION_NIL31_SHORT_TWO_THIRDS_CANONICAL_20261107.json']
 for name,file in [('T39','Q51_UNIVERSAL_SECTION_20261101.json'),('T40','Q52_UNIVERSAL_SECTION_20261101.json')]:
  path=ROOT/'reports'/file;paths.append(path);graph=json.loads(path.read_text());idx=graph['weight2_u_indices_zero_based'];residual=[s.cancel(point[k]-s.sympify(q).subs(dict(zip(us,point)),simultaneous=True))for k,q in zip(idx,graph['quadratic_section_q'])];assert len(residual)==45
  records[name]=dict(in_entire_reduced_component_graph=all(a==0 for a in residual),all45_residuals=list(map(str,residual)),nonzero_residuals={str(k):str(a)for k,a in zip(idx,residual)if a!=0})
 out=dict(task_id='E_multiple_hyperplane',all_passed=True,anchor='t=d=1 of the canonical entire nil31 short2/3 source',original_symplectic_evaluation=[1,0,1,1,0],translation=[[str(a)for a in row]for row in T.tolist()],evaluation=[[str(a)for a in row]for row in value.tolist()],evaluation_det=str(value.det()),Darboux_frame=[[str(a)for a in row]for row in frame.tolist()],perspective=list(map(str,normalized['perspective'])),all61_S61_parameters=list(map(str,point)),all150_reconstruction_and134_equations_checked=True,graph_membership=records,whole_source_implication='An exact positive full graph result plus the accepted dense single PGL/scalar orbit td nonzero covers the whole six-parameter source by closedness; no component uniqueness or smoothness is needed. Negative results exclude only T39/T40.',independent_review_required_if_positive=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_bytes=Path(__file__).read_text(),seconds=time.monotonic()-started)
 path=ROOT/'reports/CLASSIFICATION_NIL31_SHORT_ELLIPTIC_GRAPH_20261107.json';assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(graphs={k:v['in_entire_reduced_component_graph']for k,v in records.items()},seconds=out['seconds'])))
if __name__=='__main__':main()
