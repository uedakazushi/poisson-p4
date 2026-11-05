#!/usr/bin/env python3
"""Tiny cached quartic-dual gate; no ranksearch or cubevalue evaluation."""
import hashlib,itertools,json,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_SOURCE_PREFLIGHT_20261105_v63'
PATHS={
 'cubic':'reports/AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_20261105_v63_REVIEW/independent-homogeneous600-closure-and-r-reduction-proof.json',
 'cubicpeer':'reports/AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_20261105_v63_REVIEW/report.json',
 'natural':'reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63/actual-natural-Told-to-Tnil-source-map-and-origin.json',
 'naturalpeer':'reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63_REVIEW/report.json',
 'oldLoewy':'reports/AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63/report.json',
 'newLoewypeer':'reports/AMBIENT_B18_TNIL_EXACT_LOEWY16_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(v):
 if isinstance(v,Q):return str(v)
 if isinstance(v,dict):return {k:enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def field(path,key):
 text=path.read_text();s=json.dumps(key);assert text.count(s)==1;i=text.index(':',text.index(s))+1
 while text[i].isspace():i+=1
 return json.JSONDecoder().raw_decode(text,i)[0]
def main():
 t=time.monotonic();out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);C=json.loads((ROOT/PATHS['cubic']).read_text());assert all(json.loads((ROOT/PATHS[k]).read_text())['all_passed']for k in ['cubicpeer','naturalpeer','oldLoewy','newLoewypeer'])
 q=[{int(i):Q(c)for i,c in row.items()}for row in C['known_cubic_quotient3_by600']];E=C['r_multiplication_degree3_representative_rows600'];assert E==[119,479,599]
 assert [[row.get(i,Q(0))for i in E]for row in q]==[[Q(i==j)for j in range(3)]for i in range(3)]
 origin=field(ROOT/PATHS['natural'],'actual_origin_data');M=[[Q(c)for c in col]for col in origin['actual_origin_map_columns5']];assert len(M)==5 and all(len(col)==5 for col in M)
 m2=list(itertools.combinations_with_replacement(range(8),2));m3=list(itertools.combinations_with_replacement(range(8),3));m4=list(itertools.combinations_with_replacement(range(8),4));idx3={m:i for i,m in enumerate(m3)};basis_heads=[0,3,4]
 p4=[];witness=[]
 for target in range(5):
  for mon in m4:
   triple=120*target+idx3[mon[:3]];var=mon[3];a=[row.get(triple,Q(0))for row in q];value=[Q(0)]*3;steps=[]
   for f,c in enumerate(a):
    if not c:continue
    second=120*basis_heads[f]+idx3[tuple(sorted((7,7,var)))];b=[row.get(second,Q(0))for row in q]
    value=[x+c*y for x,y in zip(value,b)];steps.append([f,c,second])
   p4.append([[i,c]for i,c in enumerate(value)if c]);witness.append([triple,var,steps])
 assert len(p4)==1650
 e4=[330*h+m4.index((7,7,7,7))for h in basis_heads];assert [[dict(col).get(i,Q(0))for i in range(3)]for col in[p4[j]for j in e4]]==[[Q(i==j)for i in range(3)]for j in range(3)]
 pure=[[row.get(120*h+idx3[(7,7,7)],Q(0))for h in range(5)]for row in q]
 data={'source_known_cubic_quotient3_by600_sparse':[list(row.items())for row in q],'known_pure_r3_basis_rows600':E,'known_pure_r4_basis_rows1650':e4,
  'all1650_constructive_quartic_reduction_columns3_sparse':p4,'all1650_two_cubic_relation_witnesses':witness,'all9_r4_unit_embedding_checks':True,
  'source_q3_restricted_to_five_pure_r3_heads':pure,'source_actual_natural_origin_columns5':M,
  'future_three_by_five_relation_operator':'q3_pure3_by5 * transpose(M_columns5); do notrankinSOURCEbuilder',
  'actual_relation_proof':'Oldm4Told=0exactly. ApplyactualS-linear naturalmap: eachmappedoldhead*r4=0. Itsdegree4initialisM0head*r4, sinceallhighermapcoefficientscontribute>=5. Constructivep4givesfiverelationsamongthreepure-r4representatives.',
  'constructive_surjectivity_proof':'Foranydegree4monomial,firstq3reductiondiffersbyknownK3*x;thenx*(headf*r3)=r*(headf*r2*x),secondq3reductiondiffersbyr*knownK3. Thusactualgr4isgeneratedbythreepure-r4classes. p4neednotbeassumedwell-definedonlargerknownquotient.',
  'exact_phi4_completion_criterion':'If3x5rank2, actualgr4dimension1 forcesitskernelEXACTfiverelationspan. Aleftannihilatorlambda3thenlambda3*p4isgenuinecompleteactualquarticdual. Rank<2isinconclusive; rank3contradictsknownHF4andmuststop.',
  'future_numeric_cap':24000,'future_static_numeric_envelope':18000,'future_hard_seconds':60,'future_memory_bytes':2<<30,'future_input_plus_output_bytes':32<<20,'future_bit_guard':2048,
  'rank_or_actual_phi4_or_cube_values_computed':False,'no_new_physical_source_gauge_or_jet':True}
 ep=out/'tiny-five-relations-and-constructive-quartic-dual-source.json';ep.write_text(json.dumps(enc(data),separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'status':'SOURCE_ONLY_tiny_cached_quartic_dual_rank_gate','input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'actual_rank_dual_or_cube_computed':False,'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
