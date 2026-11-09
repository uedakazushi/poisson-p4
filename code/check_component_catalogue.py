#!/usr/bin/env python3
"""Validate the main ID crosswalk and all exact serialized general-family maps.

Sample evaluations check serialization and the coordinate adapter, not universal
Jacobi identities, density, componenthood, or the final classification theorem.
"""
from pathlib import Path
from fractions import Fraction
import argparse,gzip,hashlib,json,re
import sympy as s
from s61_model import S61,jacobi_projective
ROOT=Path(__file__).resolve().parents[1]
def sha(b):return hashlib.sha256(b).hexdigest()
def run(output):
 records=json.loads((ROOT/'research/COMPONENTS.json').read_text())['components']
 presentation=json.loads((ROOT/'code/data/component_manuscript_catalogue.json').read_text())['records']
 assert [r['id']for r in records]==[r['id']for r in presentation]
 source=(ROOT/'papers/component_family_catalogue.tex').read_text()
 assert len(set(re.findall(r'\\label\{catalogue:([^}]+)\}',source)))==len(records)==147
 main='\n'.join(p.read_text()for p in (ROOT/'papers').glob('*.tex'))
 all_labels=set(re.findall(r'\\label\{([^}]+)\}',main))
 for p,c in zip(presentation,records):
  assert p['construction_label']in all_labels and p['component_label']in all_labels and p['distinctness_label']in all_labels and p['main_label']in all_labels
  assert p['projective_dimension']==c['projective_dimension']and p['generic_rank']==c['generic_rank']
 maps=json.loads((ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json').read_text())['records']
 model=S61();results=[]
 K=s.QQ.algebraic_field(s.sqrt(3))
 for r in maps:
  path=ROOT/r['map']['path'];raw=path.read_bytes();assert sha(raw)==r['map']['sha256'];decoded=gzip.decompress(raw);assert sha(decoded)==r['decoded_sha256']
  data=json.loads(decoded);sample=data.get('exact_rational_sample',data.get('exact_sample'));assert sample
  quadratic=r['component_id']=='R4_DPM'
  def number(n):return K.from_sympy(s.sympify(n))if quadratic else Fraction(str(n))
  inputs=[number(sample[k])for k in data['input_names']];values=[];inverses=0
  for node in data['nodes']:
   op,arg=node
   if op=='q':v=number(arg)
   elif op=='v':assert 0<=arg<len(inputs);v=inputs[arg]
   elif op in ['+','*']:
    assert all(0<=j<len(values)for j in arg)
    v=number('0'if op=='+'else'1')
    for j in arg:v=v+values[j]if op=='+'else v*values[j]
   elif op=='inv':
    assert 0<=arg<len(values)and values[arg]!=number('0');v=number('1')/values[arg];inverses+=1
   else:raise ValueError(op)
   values.append(v)
  output_indices=data['output_u1_through_u61'];assert len(output_indices)==61
  u=[K.to_sympy(values[j])if quadratic else s.Rational(values[j].numerator,values[j].denominator)for j in output_indices]
  P=model.reconstruct(u)
  residual=jacobi_projective(P);assert not any(residual)
  result={'id':r['component_id'],'map':str(path.relative_to(ROOT)),'encoded_sha256':sha(raw),'decoded_sha256':sha(decoded),'nodes':len(values),'nonzero_inverse_arguments':inverses,'sample_projective_jacobi':'PASS','sample_scope':'exact implementation/convention check only; not the general-source proof'}
  results.append(result);print(r['component_id'],'PASS',flush=True)
 report={'scope':__doc__,'id_count':len(records),'one_to_one_labels':'PASS','all_construction_references_resolve':'PASS','maps':results,'general_proof_dependencies':['research/notes/MPS_S61_CATALOGUE_20261103_v62.md','research/notes/KNOWN51_RANK4_COORDINATE_MAPS_20261103_v62.md'],'script_sha256':sha(Path(__file__).read_bytes())}
 output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);run(p.parse_args().output)
