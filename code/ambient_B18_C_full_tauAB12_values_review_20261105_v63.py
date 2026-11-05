#!/usr/bin/env python3
"""Independent expanded-polynomial tauAB values, fixed minors and e3 proof."""
import hashlib,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
import ambient_B18_C_full_eta4_cache_values_review_20261105_v63 as ar
import ambient_B18_C_full_eta4_cache_source_review_20261105_v63 as met
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_FULL_TAUAB12_VALUES_20261105_v63';ID=PRI+'_REVIEW'
PRIMARY='627a847c1ff78872b5c214231bd20d75569a11f0eb7a6bb241d1c52dfebd4ef5'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};raw=0
 def pin(p,h=None):
  nonlocal raw
  p=pathlib.Path(p);v=sha(p)
  if h:assert v==h,(str(p),v,h)
  key=str(p.relative_to(ROOT))
  if key not in pins:raw+=p.stat().st_size
  pins[key]=v
  return json.loads(p.read_text()) if p.suffix=='.json' else None
 pr=pin(ROOT/'reports'/PRI/'report.json',PRIMARY);assert pr['all_passed']
 docs={p:pin(ROOT/p,h) for p,h in pr['input_sha256'].items()}
 peer=next(v for p,v in docs.items() if p.endswith('_REVIEW/report.json'))
 assert peer['full_formal_source_code_interface_whole_budget_PASS']
 source=next(v for p,v in docs.items() if p.endswith('FULL_TAUAB12_SOURCE_20261105_v63/report.json'))
 circuit=next(v for p,v in docs.items() if p.endswith('full-tauAB12-quartic-scalar-circuit-and-source-proof.json'))
 E=pin(ROOT/'reports'/PRI/'all-full-tauAB12-values-and-fixed-known-rank-C4-certificates.json',
   'ba47fe0be5110583db0b792f507f873f418c9ed7699de900cf3adf81e260c18b')
 assert sha(ROOT/'code/ambient_B18_C_full_tauAB12_values_20261105_v63.py')==pr['code_sha256']
 assert sha(ROOT/'research/notes'/f'{PRI}.md')==pr['note_sha256']
 formal=pin(ROOT/'reports/AMBIENT_B18_C_FULL_TAUAB12_SOURCE_20261105_v63_REVIEW/independent-full-tauAB-scalar-symbol-polynomial-proof.json',
   '8cc36264674b8bb2f33541a5a2cd9f6457d89f798cffeed2699793edfa4a1480')
 symbols=circuit['provided_inputs'];groups={}
 for j,row in enumerate(symbols):groups.setdefault(row['source_path'],[]).append((j,row))
 values=[None]*39;readbacks=[];peak=0
 for path,rows in sorted(groups.items()):
  d=pin(ROOT/path,source['input_sha256'][path]);peak=max(peak,met.nc(d))
  for j,row in rows:
   kind,*args=row['key']
   if kind=='delta':q=dict(next(r['delta'] for r in d['physical12_coefficients'] if r['physical_row']==args[0])).get(args[1],'0')
   elif kind=='nu_sparse33':q=dict(next(v for i,v in d['all_requested164_root_values'] if i==args[0])).get(args[1],'0')
   elif kind=='nu_ell2':q=d['actual_nuA_B_pure_ell3'][args[0]][2]
   elif kind=='lambda':q=next(r['delta_r_ell2'] for r in d['all3_lambda_readbacks'] if r['physical_row']==args[0])
   elif kind=='mu03':q=next(r['delta_ell3'] for r in d['all3_pivot_mu03_readbacks'] if r['physical_row']==args[0])
   else:raise AssertionError(kind)
   assert Q(q)==Q(row['value']);values[j]=Q(q)
   readbacks.append([row['node'],row['key'],str(Q(q)),path,row['locator']])
  del d
 assert readbacks==E['literal39_cached_input_checks'] and all(v is not None for v in values)
 # New prior bound for the DIFFERENT expanded-polynomial evaluation.
 # All terms and partials have a positive absolute sum and common D;
 # signed responses are not evaluated until this receipt is saved.
 polys=[p for rows in formal['all3_fullNF_symbol_polynomials']+formal['all15_residual_symbol_polynomials'] for p in rows]
 maxH=0
 for terms in polys:
  U=0;D=1
  for mon,c in terms:
   q=Q(c);tu=met.ceiling(q);td=q.denominator
   for j in mon:tu*=met.ceiling(values[j]);td*=values[j].denominator
   U+=tu;D=math.lcm(D,td)
   maxH=max(maxH,(U*D).bit_length(),D.bit_length())
 assert 2*maxH+1<=2048
 owned=sum(met.nc(d) for d in docs.values())+met.nc(formal)+met.nc(E)+met.nc(pr)+peak+39+16000
 assert owned<60000 and raw<2<<20
 save(out/'pre-arithmetic-expanded-polynomial-value-bound.json',{'method':'independent expanded39-symbol coefficient polynomials, not producer954-node DAG',
  'owned_upper':owned,'positive_polynomial_height_upper':maxH,'preop_upper':2*maxH+1,
  'all39_literal_inputs_verified':True,'hard_bits':2048,'hard_numeric':60000,'global_numeric':150000,
  'hard_seconds':60,'hard_memory_bytes':2<<30,'RAW_read':raw})
 ar.premax=ar.valmax=0
 values=list(map(ar.parse,map(str,values)))
 def evaluate(terms):
  total=Q(0)
  for mon,c in reversed(terms):
   value=ar.parse(c)
   for j in reversed(mon):value=ar.mul(value,values[j])
   total=ar.add(total,value)
  return total
 norm=[[evaluate(p) for p in row] for row in formal['all3_fullNF_symbol_polynomials']]
 residual=[[evaluate(p) for p in row] for row in formal['all15_residual_symbol_polynomials']]
 assert [[str(q) for q in row] for row in norm]==E['full6_by3_tauAB_columns_as3columns']
 assert [[str(q) for q in row] for row in residual]==E['all15_residual_values_as3_by5']
 v=list(map(ar.parse,circuit['known_full_HB_vector']))
 assert residual[0][2]==Q(1331,9216)
 W=norm[0];rows=[0,3];block=[[v[0],W[0]],[v[3],W[3]]]
 det=ar.add(ar.mul(block[0][0],block[1][1]),-ar.mul(block[0][1],block[1][0]))
 assert det==ar.mul(v[0],residual[0][2])==-Q(1331,15360)
 def divide(a,b):
  assert b;ar.check(max(abs(b.numerator).bit_length(),b.denominator.bit_length()))
  return ar.mul(a,ar.note(Q(b.denominator,b.numerator)))
 inverse=[[divide(block[1][1],det),divide(-block[0][1],det)],[divide(-block[1][0],det),divide(block[0][0],det)]]
 for A,B in [(block,inverse),(inverse,block)]:
  assert [[ar.add(ar.mul(A[i][0],B[0][j]),ar.mul(A[i][1],B[1][j])) for j in range(2)] for i in range(2)]==[[1,0],[0,1]]
 cert=E['fixed_rank_and_all4column_certificate']
 assert cert['source_pivot_index']==0 and cert['NF_rows']==rows and cert['minor']==str(det)
 assert cert['inverse']==[[str(q) for q in row] for row in inverse]
 # Independently derive all four coordinates from the fixed inverse;
 # verify every one of their six original target entries.
 factors=[]
 for X in [v]+norm:
  ab=[ar.add(ar.mul(inverse[i][0],X[0]),ar.mul(inverse[i][1],X[3])) for i in range(2)]
  assert all(ar.add(ar.mul(v[i],ab[0]),ar.mul(W[i],ab[1]))==X[i] for i in range(6))
  factors.append(list(map(str,ab)))
 assert factors==cert['all4_actual_source_columns_factor']==[['1','0'],['0','1'],['275/32','25/4'],['55/24','1']]
 # Literal four-source cut minor, using all lower target entries; no
 # Schur determinant assumption or truncation of unknown lower blocks.
 C4=[[v[i],W[i],Q(i==2),Q(i==5)] for i in [0,3,2,5]]
 dp={0:Q(1)}
 for i in range(4):
  nxt={}
  for mask,z in dp.items():
   for j in range(4):
    if mask>>j&1:continue
    sign=(-1)**sum(mask>>k&1 for k in range(j+1,4));nm=mask|1<<j
    nxt[nm]=ar.add(nxt.get(nm,Q(0)),ar.mul(z,ar.mul(Q(sign),C4[i][j])))
  dp=nxt
 assert dp[15]==det
 assert E['C4_rank4_certificate_if_positive']['universal_block_minor']==str(det)
 assert pr['exact_comparison_C_e3_proved'] and pr['actual_known_HB_plus_tauAB_rank']==2
 proof={'all39_inputs_rebound_and_no_source_producer_replay':True,
  'different_method_expanded_symbol_polynomial_evaluation':True,
  'all18_fullNF_and15_residual_values_exact':True,
  'full6_by3_columns':[[str(q) for q in row] for row in norm],
  'all15_residuals':[[str(q) for q in row] for row in residual],
  'positive_actual_source':'physical18 tauA tauB with full mu03ellP2 correction',
  'residual_fullNFindex3':'1331/9216','fixed_rows0_3_minor':str(det),
  'both_inverse_products_and_all4_six_target_factorbacks':True,'all4_factors':factors,
  'literal_C4_four_source_full_matrix':[[str(q) for q in row] for row in C4],
  'literal_C4_subset_DP_minor':str(dp[15]),
  'mathematical_completion':[
   'Two genuine full quartic initials independent in knownNF6 imply actualHF4<=4. D-rank lower4 givesHF4=4, hence b2 by accepted M2-type(2,b).',
   'Together with literalell4A/B, the actual four-source minor is nonzero. KnownC4NF6 loses4dimensions, so HF4C4<=2.',
   'C4prefix3666 and exactm5C4zero give length<=23. AcceptedDVRlower23 forces length23, HF36662 and e3.',
   'Thus fullcomparisonC has D-type free4 plus(1,1,2,3), torsionlength7. Its exactfullHF is3666444... by the ell-reduction/tail theorem.'
  ],
  'actual_full_HF':[3,6,6,6,'4 at every degree>=4'],
  'actual_full_series':'(3+3t-2t^4)/(1-t)',
  'actual_comparison_C_DVR_type':{'free_rank':4,'torsion_exponents':[1,1,2,3],'torsion_length':7},
  'actual_C4_length':23,'actual_C4_HF':[3,6,6,6,2,0],
  'no_original_H0_nil_cube_nil_index_or_exhaustion_transfer':True,
  'own_before_operation_bits':ar.premax,'own_stored_QQ_bits':ar.valmax,'owned_upper':owned,
  'pre_arithmetic_positive_bound':maxH,'hard_bits':2048,'hard_numeric':60000,'hard_global_numeric':150000}
 ep=out/'independent-expanded-tauAB-values-fixed-minors-and-exact-C-type1123.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'input_sha256':pins,'code_sha256':sha(__file__),
  'evidence_sha256':{ep.name:sha(ep)},'all18NF_15residuals_and_certificates_verified':True,
  'exact_comparison_C_e3_type1123_and_C4_length23_accepted':True,
  'actual_full_C_HF3666444_accepted':True,'elapsed_seconds':time.monotonic()-start}
 save(out/'report.json',report);signal.alarm(0)
 print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':report['elapsed_seconds']}))
if __name__=='__main__':main()
