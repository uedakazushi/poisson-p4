#!/usr/bin/env python3
"""Natural torus boundary of the complete doubleA magnetic source."""
from pathlib import Path
import hashlib,json,time
import sympy as s
from classification_graded_poisson_frame import graded_frame
ROOT=Path(__file__).resolve().parents[1]
def main():
 started=time.monotonic();xx=s.symbols('h x1 x2 x3 x4');params=s.symbols('m0:6');env={str(v):v for v in xx+params};sourcepath=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json';source=json.loads(sourcepath.read_text())['cases'][0];assert source['name']=='nil22_image_double_A';whole=s.zeros(5)
 for ij,text in source['universal_brackets'].items():
  i,j=map(int,ij.split(','));whole[i,j]=s.sympify(text,locals=env);whole[j,i]=-whole[i,j]
 values=dict.fromkeys(params,0);P=whole.subs(values);directions=[whole.diff(v).subs(values)for v in params];frame=graded_frame(P,directions,xx,(0,1,1,0,0),-1);frame.update(source_name=source['name'],all6_source_directions_retained=True,canonical_magnetic_boundary='F0=0',finite_strict_positive_IFT_available=all(w>0 for w in frame['true_normal_relative_weights']))
 inputs=[Path(__file__),ROOT/'code/classification_graded_poisson_frame.py',ROOT/'code/classification_rank3_mixed_component_review.py',sourcepath,ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_MAGNETIC_NORMALFORMS_20261107.json'];result=dict(task_id='E_multiple_hyperplane',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=Path(__file__).read_text(),frame=frame,scope='Theuniquegradedmagneticboundary ofacceptedwhole6source, notrandompoint. IfalltrueNweightspositive fullfiniteIFTcanbeused; ifzero/mixed doNOTtruncate. No componenthoodor membership fromframealone.',seconds=time.monotonic()-started)
 out=ROOT/'reports/CLASSIFICATION_NIL22_DOUBLE_A_BOUNDARY_FRAME_20261107.json'
 if out.exists():raise FileExistsError(str(out))
 out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:frame[k]for k in ['D_rank','actual_G_rank','actual_full_source_F_rank','true_normal_relative_weights','finite_strict_positive_IFT_available','seconds']}))
if __name__=='__main__':main()
