#!/usr/bin/env python3
"""Combine exact modular, horizontal and vertical cubic certificates."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations
from sympy import GF
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1]

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def verify():
 start=time.time();names=['S61_BINARY_HYPERPLANE_CUBIC_MODULAR_20261102','S61_BINARY_HYPERPLANE_CUBIC_HORIZONTAL_20261102','S61_BINARY_HYPERPLANE_CUBIC_VERTICAL_20261102','S61_BINARY_HYPERPLANE_GENERIC_RANKS_20261102']
 paths=[ROOT/'reports'/f'{name}.json' for name in names];reports=[json.loads(p.read_text()) for p in paths];primary={}
 def audit(obj):
  if isinstance(obj,dict):
   for name,expected in obj.get('input_sha256',{}).items():
    assert sha(ROOT/name)==expected
    if name in primary:assert primary[name]==expected
    primary[name]=expected
   for key,value in obj.items():
    if key!='input_sha256':audit(value)
  elif isinstance(obj,list):
   for value in obj:audit(value)
 for report in reports:assert report['all_passed'];audit(report)
 modular,horizontal,vertical,generic=reports;candidate_path=ROOT/'code/certificates/s61_binary_hyperplane_cubic_first_minor_20261102.json';candidate=json.loads(candidate_path.read_text());matrixhash=candidate['source_schur_sha256']
 assert modular['source_matrix_sha256']==horizontal['source_matrix_sha256']==vertical['schur_matrix_sha256']==matrixhash
 assert horizontal['first_minor']==candidate
 assert modular['choices'][0]['rows']==candidate['rows'] and modular['choices'][0]['columns']==candidate['columns']
 assert modular['maximal_minor_gcd_unit_on_known_open'] and modular['each_unknown_factor_irreducible_modp']
 assert modular['first_known_orders_a_F_B_C']==[37,10,12,26]
 assert len(modular['unknown_common_factors'])==len(modular['full_rank_witnesses'])==3
 assert horizontal['hungarian_upper']==107 and horizontal['d_kernel_coranks']==[14,9,5]
 assert horizontal['exact_d_degree']==79 and horizontal['d_profiles'][-1]['rank']==36
 assert horizontal['B_initial']['rank']==11 and horizontal['B_projected']['rank']==4
 assert horizontal['B_partial_factor_audit']['net_B_order']==2
 assert horizontal['B_partial_order_lower']==7+3==10 and horizontal['B_total_order_lower']==12
 assert vertical['raw_matrix_coefficients_p_integral'] and vertical['lex_e_a_leading_coefficient_is_p_unit']
 assert vertical['specialized_d']==3 and vertical['e_degree']==213 and vertical['leading_monomial_a_e']==[39,213]
 assert {k:v['valuation'] for k,v in vertical['known_factor_valuations'].items()}=={'a':37,'F':10,'C':26}
 assert vertical['matrix_reconstruction']['known_unit_pivots']==21
 assert generic['all_core_cubic_rank_upper']==57
 # Primitive known factors over the residue ring, without guessing
 # irreducibility from numerical points. F,C are linear in a with constant
 # unit coefficient; B is primitive and linear in d.
 p=32003;R,a,e,d=ring('a,e,d',GF(p));F=500*a+90*e**2+3*e-24;B=(50*a+18*e**2-9*e-2)*d+60*e-40;C=100*a+18*e**2-9*e-2
 assert (50*a+18*e**2-9*e-2).gcd(60*e-40)==R.one
 assert F.degree(a)==C.degree(a)==1 and B.degree(d)==1
 assert all(x.gcd(y)==R.one for x,y in combinations([a,F,B,C],2))
 note=ROOT/'research/notes/S61_BINARY_HYPERPLANE_CUBIC_CODIMENSION_20261102.md'
 inputs=[Path(__file__),note,candidate_path,*paths]
 return dict(id='S61_BINARY_HYPERPLANE_CUBIC_CODIMENSION_20261102',all_passed=True,prime=p,source_matrix_sha256=matrixhash,full_cubic_matrix_shape=[63,91],known_unit_pivots=21,residual_matrix_shape=[42,70],generic_cubic_rank=57,all_core_cubic_rank_upper=57,maximal_minor_gcd_unit_over_QQ_on_aFBC_open=True,cubic_rank_drop_codimension_lower=2,cubic_rank_drop_dimension_upper=1,known_factors_pairwise_coprime_modp=True,B_primitive_and_linear_in_d=True,horizontal_degree=79,horizontal_B_order=12,vertical_e_degree=213,vertical_lex_leading_monomial_a_e=[39,213],vertical_orders_a_F_C=[37,10,26],primary_input_sha256=primary,primary_inputs_checked=len(primary),proof_note=str(note.relative_to(ROOT)),input_sha256={str(path.relative_to(ROOT)):sha(path) for path in inputs},limitations=['Nuisance zero on the three-core open a*F*B*C != 0.','Rank-drop curves and points remain; all-point rank57 is not asserted.','No quartic surjectivity, relative length65, ambient t-nilpotence, or S61 exhaustion follows from this certificate alone.'],seconds=time.time()-start)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();assert not args.output.exists();result=verify()
 with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps({k:result[k] for k in ['all_passed','cubic_rank_drop_dimension_upper','primary_inputs_checked','seconds']}))
