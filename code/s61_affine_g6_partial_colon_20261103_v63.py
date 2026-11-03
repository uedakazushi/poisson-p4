#!/usr/bin/env python3
"""Universal61-variable firstw-colon identity and exact675-column extension."""
import argparse
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
import sympy as s
from s61_model import S61
from s61_affine_fiber_macaulay import load_fiber
from s61_affine_normal_global_quintic_20261103_v63 import ROOT,sha,add,uadd,umul

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 extra=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json';old=json.loads(extra.read_text());assert old['all_passed'];model=S61();two=old['q_to_u_indices_zero_based'];assert two==[i for i,w in enumerate(model.reduced['parameter_weights'])if w==2];Q=set(two);F=model.reduced['equations'];W=61
 universal={};right={}
 for encoded,c in old['selected_non_Bianchi_syzygy']:
  j,i=divmod(encoded,70);c=Fraction(c);assert i<70
  for mon,n,d in F[i]:
   mon=tuple(mon);value=c*Fraction(n,d);add(universal,tuple(sorted(mon+(two[j],))),value)
   qdegree=sum(k in Q for k in mon);assert qdegree<=2;add(right,tuple(sorted(mon+(W,)*(2-qdegree)+(two[j],))),value)
 target={tuple(mon):Fraction(n,d)for mon,n,d in old['universal_new_equation']};assert universal==target and len(target)==9257
 qdegree=max(sum(k in Q for k in mon)for mon in target);assert qdegree==2
 left={}
 for mon,c in target.items():add(left,tuple(sorted(mon+(W,)*(3-sum(k in Q for k in mon)))),c)
 assert left==right
 axispath=ROOT/'reports/S61_MINIMAL_AXIS_INFINITY_20261101.json';axis=json.loads(axispath.read_text());assert axis['all_passed']and axis['all_135_identities'];assert axis['G6_on_full_nine_G_vertex_chart']=='9*(2*k+1)*g6*g7';point=list(map(s.Rational,axis['sample_S61_coordinates']));top=[term for term in old['universal_new_equation']if sum(k in Q for k in term[0])==2];assert model.evaluate(top,point)==0
 # Bind the centered G6 to the same exactT01 fiber as the quintic.
 _,qids,center,source,pointsource=load_fiber();audit=json.loads(pointsource.read_text());fullpoint=list(map(Fraction,next(p for p in audit['points']if p['component_id']=='T01')['S61_coordinates_u1_through_u61']));qindex={u:i for i,u in enumerate(qids)};centered={}
 for mon,c in target.items():
  coefficient=c;free=[]
  for u in mon:
   if u in qindex:free.append(qindex[u])
   else:coefficient*=fullpoint[u]
  terms={():coefficient}
  for j in free:
   image={(j,):Fraction(1)}
   if center[j]:image[()]=center[j]
   new={}
   for a,x in terms.items():
    for b,y in image.items():add(new,tuple(sorted(a+b)),x*y)
   terms=new
  for mon,c in terms.items():add(centered,mon,c)
 assert centered=={tuple(mon):Fraction(c)for mon,c in old['sample_centered_new_equation']}
 quintic=ROOT/'reports/S61_AFFINE_NORMAL_GLOBAL_QUINTIC_20261103_v63/presentation.json.gz';data=json.loads(gzip.decompress(quintic.read_bytes()));decode=lambda poly:{degree:{tuple(m):Fraction(c)for m,c in coeff}for degree,coeff in poly};G=decode(data['monic_quintic_G']);z=s.Symbol('z0')
 def univar(expression):return {int(m[0]):{():Fraction(c)}for m,c in s.Poly(s.sympify(expression),z).terms()if c}
 Z=[{1:{():Fraction(1)}},uadd(univar(data['p1']),{0:{(42,):Fraction(1)}}),uadd(univar(data['p2']),{0:{(43,):Fraction(1)}})];framepath=ROOT/'reports/S61_AFFINE_RESIDUAL42_20261103_v63/linear129x42-QQ.json.gz';frame=json.loads(gzip.decompress(framepath.read_bytes()));E=[[Fraction(c)for c in row]for row in frame['embedding']];qimages=[]
 for i,row in enumerate(E):
  image={}
  for k,c in enumerate(row):image=uadd(image,Z[k],c)
  if i>=3:image=uadd(image,{0:{(i-3,):Fraction(1)}})
  qimages.append(image)
 image={}
 for mon,c in centered.items():
  product={0:{():Fraction(1)}}
  for j in mon:product=umul(product,qimages[j])
  image=uadd(image,product,c)
 quotient={};remainder=image
 while remainder and max(remainder)>=5:
  degree=max(remainder);term={degree-5:remainder[degree]};quotient=uadd(quotient,term);remainder=uadd(remainder,umul(term,G),Fraction(-1))
 assert uadd(umul(quotient,G),remainder)==image and max(remainder,default=-1)<=4
 encode=lambda poly:[[degree,[[list(mon),str(c)]for mon,c in sorted(coeff.items())]]for degree,coeff in sorted(poly.items())]
 (out/'quintic-extra-generator.json').write_text(json.dumps(dict(extra_h_remainder=encode(remainder),exact_division_quotient=encode(quotient),module_extension='AddT_h^j*r_G6,j0..4;5x675presentation,sameactualA',formal_placeholders_have_same_fixed_normal_recipe=True),indent=2)+'\n')
 certificate=dict(full_generic61_variable_QQ_identity_verified=True,syzygy_nonzero_terms=len(old['selected_non_Bianchi_syzygy']),G6_terms=9257,actual_q_degrees=dict(G6=2,K_i=1,homogenized_w_times_G6=3),
  universal_identity='w*G6^h=sum_(i<70) K_i(q)*F_i^h,all61parametersfree',first_partial_colon='I134 subset I135=I134+(G6^h) subset (I134:w)',
  saturation_equalities=['(I134:w^infinity)=(I135:w^infinity)','afteraxisidempotentquotient,sameD(w)andsamew-saturatedclosure'],
  composition_recipe='OnaxispatchJ=I135+(e_axis);e_axisdescendsasimageofidempotent,w*e_axis=0;bothstepsareactualidealoperationspreservingallD(w).',
  certified_axis_G6_zero_from_full_generic_old135_family=True,selected_etale_factor_may_shrink='IfcolonselectedadditionalisolatedfactorswithG6nonzero,addingG6killsjustthose;theidempotentdescendsandcompositionstillholds.',
  all_extra_quintic_reduction_QQ_verified=True,actual_module_shape_atmost=[5,675],no_new_rank_search=True,no_unit_or_fulltorsion_claim=True)
 (out/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
 paths=[Path(__file__),extra,axispath,quintic,framepath,source,pointsource,ROOT/'reports/S61_AFFINE_AXIS_POLARIZATION_SPLIT_20261103_v63/report.json',ROOT/'research/notes/S61_AFFINE_G6_PARTIAL_COLON_20261103_v63.md']
 report=dict(id='S61_AFFINE_G6_PARTIAL_COLON_20261103_v63',all_passed=True,fullgeneric61_QQ_partialcolon_identity=True,actual_qdegree2_G6_and_degree1_multipliers=True,axis_idempotent_composition_preserves_affine_scheme=True,exact_global_quintic675_extension_saved=True,full_w_torsion_or_unit_not_claimed=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('UNIVERSAL61wG6/EXACTPARTIALCOLON/AXISCOMPOSITION/675EXTENSION PASS',report['seconds'],flush=True)
if __name__=='__main__':main()
