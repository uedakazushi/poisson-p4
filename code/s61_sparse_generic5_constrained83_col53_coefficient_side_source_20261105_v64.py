#!/usr/bin/env python3
"""SOURCE: collect row multipliers before the next fixed raw-column contraction."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_COEFFICIENT_SIDE_SOURCE_20261105_v64'
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def variable(x):return{((x,1),):Q(1)}
def scale(A,c):return{e:v*c for e,v in A.items()if v*c}
def add(A,B):
 C=dict(A)
 for e,c in B.items():
  C[e]=C.get(e,Q(0))+c
  if not C[e]:del C[e]
 return C
def mul(A,B):
 C={}
 for e,c in A.items():
  for f,d in B.items():
   z=dict(e)
   for x,k in f:z[x]=z.get(x,0)+k
   g=tuple(sorted((x,k)for x,k in z.items()if k));C[g]=C.get(g,Q(0))+c*d
 return{e:c for e,c in C.items()if c}
def mon(c,**v):return{tuple(sorted((x,k)for x,k in v.items()if k)):Q(c)}
def times(*A):
 C={():Q(1)}
 for a in A:C=mul(C,a)
 return C
def encode(A):return[[[[x,k]for x,k in e],str(c)]for e,c in sorted(A.items())]
def main(a):
 assert a.freeze
 start=time.monotonic();out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);inputs={}
 def load(stem,name):
  p=ROOT/'reports'/stem/name;inputs[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 op=load('S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_TENTH_UNIT_SOURCE_V2_20261105_v64','actual-F10-pivot-and-full167x14-typed-operator.json')
 oldstem='S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64'
 old=load(oldstem,'staged-literal-col56-program-and-positive-profiles.json')
 api=load('S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_COEFFICIENTS_20261105_v64','actual15-col56-and-once-shared-row35.json');I={x['row']:x for x in api['actual15_ninth_col56_numerators']}
 c8stem='S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_COEFFICIENTS_20261105_v64'
 c8=load(c8stem,'actual-cleared-col82-and-commonQ3RE-representation.json');B={x['row']:x for x in c8['actual16_col82_numerator_values']}
 Qp,P7,E,H,s,q,t=map(variable,['Q','P7','E','H','s','q','t'])
 rowsymbols=['self','M231','M283','M214','M235','M2104','A51','A11','A44','A52']
 # Verify the unchanged stage8 identity BEFORE introducing Plucker functions.
 def weights(tag):
  U,V,C,F,N,M=[variable(x+tag)for x in['U','V','C','F','N','M']]
  aa=times(E,P7,Qp,Qp);b=scale(times(E,P7,Qp,U),-1)
  cc=scale(times(E,P7,Qp,Qp,add(times(mon(Q(1,40),H=1,t=-2),U),times(mon(Q(405,8),t=-2),V))),-1)
  d=times(mon(Q(-6075,2),H=-1),E,P7,Qp,C)
  ff=times(mon(-765450000,H=-1,q=-1),E,P7,F)
  gg=times(mon(Q(3645,64),H=-2,s=1,t=1),E,N)
  hh=times(mon(675,H=-1),E,Qp,N)
  kk=times(mon(Q(-273375,32768),H=-3,t=2),M)
  return[aa,b,cc,d,scale(d,6),ff,gg,hh,hh,kk]
 wi=weights('i');basis=[variable(x)for x in rowsymbols]
 U,V,C,F,N,M=[variable(x+'i')for x in['U','V','C','F','N','M']]
 h4=add(add(add(times(Qp,basis[0]),scale(times(U,basis[1]),-1)),scale(times(Qp,add(times(mon(Q(1,40),H=1,t=-2),U),times(mon(Q(405,8),t=-2),V)),basis[2]),-1)),times(mon(Q(-6075,2),H=-1),C,add(basis[3],scale(basis[4],6))))
 f6=add(times(Qp,h4),times(mon(-765450000,H=-1,q=-1),F,basis[5]))
 b51=add(basis[6],times(mon(Q(320,27),H=1,s=-1,t=-1),Qp,add(basis[7],basis[8])))
 n7=add(times(P7,f6),times(mon(Q(3645,64),H=-2,s=1,t=1),N,b51))
 n8=add(times(E,n7),times(mon(Q(-273375,32768),H=-3,t=2),M,basis[9]))
 contracted={}
 for x,y in zip(wi,basis):contracted=add(contracted,mul(x,y))
 assert n8==contracted
 rawrows=['i',58,31,83,14,35,104,51,11,44,52];m2basis=['M2i','M258','M231','M283','M214','M235','M2104','A51','A11','A44','A52']
 def vector(tag,selfidx):
  w=weights(tag);v=[{}for _ in range(11)];v[selfidx]=add(v[selfidx],w[0])
  for index,c in zip(range(2,11),w[1:]):v[index]=add(v[index],c)
  return v
 vi=vector('i',0);v58=vector('58',1);v35=vector('35',5)
 assert v35[5]==add(weights('35')[0],weights('35')[4])
 X,Y,Z=map(variable,['UF10','C10_Ni56','C9_C10_alpha_i'])
 combined=[add(add(mul(X,a),scale(mul(Y,b),-1)),mul(Z,c))for a,b,c in zip(vi,v58,v35)]
 # Expand M2's ONLY on the coefficient side, not by building a large column.
 raw=[{}for _ in range(11)]
 for k,w in enumerate(combined):
  if k>=7:raw[k]=add(raw[k],w);continue
  row=rawrows[k];r=str(row);raw[k]=add(raw[k],w)
  a49=variable('A'+r+'_49');a51=variable('A'+r+'_51')
  s_coef=times(mon(324000,s=-1,t=-1),a49)
  delta_coef=add(times(mon(36000,s=-1,t=-1),a49),times(mon(Q(-675,8),H=-1,q=-1),a51))
  raw[9]=add(raw[9],mul(w,add(s_coef,scale(delta_coef,-8))))
  raw[8]=add(raw[8],mul(w,delta_coef))
 # Owner map: A_i56 is a fixed existing function callback, NOT a saved value.
 owners=[]
 for row in op['selected14_row_ids']+[58]:
  node=f'row{row}_N8_sum_m';assert node in old['positive_support_and_height_profiles']
  ir=I[row];br=B[row]
  for stem,rec in [('S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_COEFFICIENTS_20261105_v64',ir),(c8stem,br)]:
   p=ROOT/'reports'/stem/rec['file'];assert sha(p)==rec['sha256'];inputs[str(p.relative_to(ROOT))]=sha(p)
  owners.append({'row':row,'N9_i56_actual_binary':ir,'B_i_actual_N8_i82_binary':br,'N8_i56_function_callback':{'literal_program':f'reports/{oldstem}/staged-literal-col56-program-and-positive-profiles.json','node':node,'profile':old['positive_support_and_height_profiles'][node],'NOT_a_new_numeric_cache':True}})
 ev={'all_passed':True,'source_only':True,'stage8_vector_identity_exact':True,'W8i_weights10':[encode(v)for v in wi],'mixed_basis11':m2basis,'W35_alias_M235_equals_a_plus6d':True,'Plucker_coefficient_side11':[encode(v)for v in combined],'raw_row_ids11':rawrows,'raw_coefficient_side11':[encode(v)for v in raw],'functional_symbols':{'UF10':'Q^4*R^2*E*F10','C10_Ni56':'353894400000*H*t0*N9_i56','C9_C10_alpha_i':'C9*C10*(B58*N8_i56-Bi*N8_58,56)'},'prefactor':'P7/(U^2*F10)','all_function_owners15':owners,'source_profile_pending':True,'no_actual_new_scalar_or_matrix_coefficients_evaluated':True,'no_numeric_callback_guess_or_F9_division':True,'all167x14_and87plus66_RHS_P290once_historicalF9_F10_complements':True}
 p=out/'exact-11-row-coefficient-operator-and-callback-owners.json';p.write_bytes(enc(ev))
 for path in[Path(__file__),ROOT/f'research/notes/{ID}.md']:inputs[str(path.relative_to(ROOT))]=sha(path)
 symbolic_scalar_upper=25000+64*sum(map(len,[n8,contracted,*wi,*combined,*raw]))+12000
 assert symbolic_scalar_upper<150000
 rawbytes=sum((ROOT/p).stat().st_size for p in inputs)+p.stat().st_size+65536;assert rawbytes<32*1024*1024
 report={'id':ID,'all_passed':True,'source_only':True,'stage8_10_weight_identity':True,'Plucker_before_column_contraction':True,'raw_rows11':True,'W35_self_alias_preserved':True,'actual_function_owners_bound':True,'actual_K5_coefficient_support_profile_pending':True,'no_new_matrix_or_scalar_VALUE_or_GO':True,'symbolic_numeric_upper':symbolic_scalar_upper,'READ_output_RAW_upper':rawbytes,'input_sha256':inputs,'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-start}
 (out/'report.json').write_bytes(enc(report));print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--freeze',action='store_true');main(p.parse_args())
