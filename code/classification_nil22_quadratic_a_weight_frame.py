#!/usr/bin/env python3
"""The existinglambda7 quadraticA orbit in its exact torus-normalized coordinates."""
from pathlib import Path
import hashlib,json,time
import sympy as s
from classification_graded_poisson_frame import graded_frame
ROOT=Path(__file__).resolve().parents[1]
def main():
 started=time.monotonic();xx=s.symbols('h x1 x2 x3 x4');params=s.symbols('m0:6')+(s.Symbol('lam'),);env={str(z):z for z in xx+params};sourcepath=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json';packet=json.loads(sourcepath.read_text())['cases'][7];assert packet['name']=='nil22_image_quadratic_A';whole=s.zeros(5)
 for ij,text in packet['universal_brackets'].items():
  i,j=map(int,ij.split(','));whole[i,j]=s.sympify(text,locals=env);whole[j,i]=-whole[i,j]
 values=dict(zip(params,[1,0,0,0,0,0,7]));P=whole.subs(values);directions=[whole.diff(v).subs(values)for v in params];frame=graded_frame(P,directions,xx,(0,1,-1,1,-1),0)
 frame.update(name=packet['name'],canonical_parameter_point={str(p):str(v)for p,v in values.items()},all7_source_directions_retained=True,anchor_relation='Sameexistinglambda7genericwitness orbit: acceptedaffinefibertranslation normalizesmagneticinvariantI=m0+m4=8 thencommonfiberscale givesF0=1. lambdaunchanged. Sourceclassnotnewpointfarm.')
 inputs=[Path(__file__),ROOT/'code/classification_graded_poisson_frame.py',ROOT/'code/classification_rank3_mixed_component_review.py',sourcepath,ROOT/'research/classification/NIL22_IMAGE_SOURCES_AND_COMPONENTS.md'];out=dict(task_id='E_multiple_hyperplane',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=Path(__file__).read_text(),frame=frame,scope='Fullsourcegradedframeonly. Ifallnormalweights0, torusequivarianttruecutIFT forcesallnonzero-weight implicitvariables identicallyzero; exactzero-weightlocalideal canbecomputedwithoutjets. No componenthoodfromframealone.',seconds=time.monotonic()-started)
 dest=ROOT/'reports/CLASSIFICATION_NIL22_QUADRATIC_A_WEIGHT_FRAME_20261107.json'
 if dest.exists():raise FileExistsError(str(dest))
 dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:frame[k]for k in ['D_rank','actual_G_rank','actual_full_source_F_rank','normal_count','true_normal_relative_weights','seconds']}))
if __name__=='__main__':main()
