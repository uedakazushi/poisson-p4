import importlib.util,json,math,hashlib,signal,resource
from pathlib import Path
from fractions import Fraction as Q

root=Path(__file__).resolve().parents[1]
signal.alarm(60)
resource.setrlimit(resource.RLIMIT_AS,(2147483648,2147483648))
p=root/'code/ambient_B18_C3_right2_lambda3_cache_source_20261105_v63.py'
sp=importlib.util.spec_from_file_location('alpha_source_bounds_only',p);S=importlib.util.module_from_spec(sp);sp.loader.exec_module(S)
c=json.loads((root/'reports/AMBIENT_B18_C3_RIGHT2_LAMBDA3_CACHE_SOURCE_20261105_v63/compiled-six-even-roots-and-provided-joint105-plane3-cache.json').read_text())
provided={i:S.input_bound(v,c['provided_monomials10'])for i,v in c['provided105_originalID_r_s_coefficients10']}
def height(v):
 return max(max(abs(u.numerator).bit_length(),u.denominator.bit_length())if isinstance(u,Q)else u.bit_length()+d.bit_length()for u,d in v)
def conv(a,b):
 v=S.zero();v[0]=(a[0][0]*b[0][0],1)
 for m in range(1,6):
  for i,j in S.PAIRS[m]:v[m]=S.bound_sum(v[m],S.bound_term(a[i],b[j]))
 return v
cache={};inv=[];diagH=0
for r in c['operator_records_originalID_op_args']:
 i,op=r[:2]
 if op=='input':v=provided[i]
 elif op=='c':v=S.zero();v[0]=(Q(r[2],r[3]),1)
 elif op=='a':v=[S.bound_sum(a,b)for a,b in zip(cache[r[2]],cache[r[3]])]
 elif op=='n':v=cache[r[2]].copy();v[0]=(-v[0][0],1)
 elif op=='m':v=conv(cache[r[2]],cache[r[3]])
 elif op=='i':
  x=cache[r[2]];assert x[0][0]==1;v=S.zero();v[0]=(Q(1),1)
  for m in range(1,6):
   for a,b in S.PAIRS[m]:
    if a:v[m]=S.bound_sum(v[m],S.bound_term(x[a],v[b]))
  dv=conv(x,v);h=height(dv);diagH=max(diagH,h);inv.append([i,r[2],h,2*h+1])
 else:raise AssertionError(op)
 cache[i]=v
delta=[]
for r in c['root_pairs_physical3']:
 a,b=cache[r['aA_node']],cache[r['aB_node']]
 d=[S.bound_sum(x,y)for x,y in zip(a,b)];h=height(d);delta.append([r['physical_row'],h,2*h+1])
# Scalar-expression majorants for provided cache conversion and all
# lower-overlap computations, without signed numerical evaluation.
def leaf(s):q=Q(s);return(S.ceil_abs(q),q.denominator)
def scalarH(a):return a[0].bit_length()+a[1].bit_length()
hist=[]
def plus(a,b):
 z=(a[0]+b[0],math.lcm(a[1],b[1]));hist.append(2*scalarH(z)+1);return z
def times(a,b):
 z=(a[0]*b[0],a[1]*b[1]);hist.append(2*scalarH(z)+1);return z
parseH=0
for i,v in c['provided105_originalID_r_s_coefficients10']:
 a=list(map(leaf,v));parseH=max(parseH,*(max(abs(Q(s).numerator).bit_length(),Q(s).denominator.bit_length())for s in v))
 low=[a[0],plus(a[1],a[2]),times(a[2],leaf('1/5')),plus(plus(a[3],a[4]),a[5]),times(plus(a[4],times(leaf('2'),a[5])),leaf('1/5')),times(a[5],leaf('1/25'))]
 last=times(plus(a[8],times(leaf('3'),a[9])),leaf('1/25'))
result={'all_passed':True,'inverse_diagnostic_bounds':inv,'delta_subtraction_bounds':delta,'input_conversion_operation_bound':max(hist),'all_input_parse_max':parseH,'overall_supplement_operation_bound':max([max(hist)]+[x[3]for x in inv]+[x[2]for x in delta]),'no_positive_root_readback':True}
assert result['overall_supplement_operation_bound']<=79
result['all_temporary_scope_bound_including_original_DAG']=79
result['source_primary_sha256']=hashlib.sha256((root/'reports/AMBIENT_B18_C3_RIGHT2_LAMBDA3_CACHE_SOURCE_20261105_v63/report.json').read_bytes()).hexdigest()
result['code_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
result['source_bound_helper_sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
result['compiled_input_sha256']=hashlib.sha256((root/'reports/AMBIENT_B18_C3_RIGHT2_LAMBDA3_CACHE_SOURCE_20261105_v63/compiled-six-even-roots-and-provided-joint105-plane3-cache.json').read_bytes()).hexdigest()

# Post-root affine/rank bounds. No prospective lambda value is computed.
op_path=root/'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63/two-coherent-HB-quartic-sources-and-unpaired-affine-operator.json'
op=json.loads(op_path.read_text())
c2_path=root/'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json'
c2all=json.loads(c2_path.read_text());c2={r['physical_row']:dict(r['delta'])for r in c2all['physical12_coefficients']}
lam={}
for rr in c['root_pairs_physical3']:
 aa,bb=cache[rr['aA_node']][5],cache[rr['aB_node']][5]
 lam[rr['physical_row']]=(aa[0]+bb[0],math.lcm(aa[1],bb[1]))
def affine_bound(terms):
 U=0;D=1;parse_bound=0
 for key,ss in terms:
  w=Q(ss);tag,physical,index=key
  if tag=='lambda':u,d=lam[physical]
  elif tag=='cache':
   value=Q(c2[physical].get(index,'0'));u,d=S.ceil_abs(value),value.denominator
  elif tag=='constant':u,d=1,1
  else:raise AssertionError(tag)
  U+=S.ceil_abs(w)*u;D=math.lcm(D,w.denominator*d)
  parse_bound=max(parse_bound,abs(w.numerator).bit_length(),w.denominator.bit_length())
 h=max((U*D).bit_length(),D.bit_length())
 return{'U':str(U),'D':str(D),'height':h,'coefficient_parse_bound':parse_bound,'term_count':len(terms)}
bounds=[[affine_bound(terms)for terms in row]for row in op['matrix2_rows_A4_B4_columns_HB0_HB1_linear_leaf_operators']]
H=max(x['height']for row in bounds for x in row)
factor_forms=[[[['constant',0,0],'-209831/2160'],[['lambda',60,-1],'80/3'],[['lambda',67,-1],'40/3']],
 [[['constant',0,0],'61/6'],[['lambda',18,-1],'100/3'],[['lambda',60,-1],'-40/3'],[['lambda',67,-1],'50']]]
fbounds=[affine_bound(t)for t in factor_forms]
post=max(12*H+3,2*max(x['height']for x in fbounds)+1)
assert post<=2048
result['lambda3_positive_coefficient_majorants']=[[r,u,d]for r,(u,d)in lam.items()]
result['four_affine_entry_majorants']=bounds
result['two_affine_factor_forms']=factor_forms
result['two_affine_factor_majorants']=fbounds
result['max_affine_entry_height']=H
result['all_postroot_affine_minor_inverse_rank1_operation_bound']=post
result['whole_root_and_post_operation_bound']=max(79,post)
result['post_bound_proof']='Each partial affine sum has abs<=U and denominator dividingD. If entries haveheight<=H, determinant prechecks<=4H+1, inverseentries<=5H+1, bothinverse-product sums<=12H+3; rank1 factors<=6H+1. Explicit factor-form checks fit the recorded smaller bounds. Every actual operation will still beprechecked at2048.'
result['no_lambda_or_quartic_values_computed']=True
result['right2_operator_sha256']=hashlib.sha256(op_path.read_bytes()).hexdigest()
result['C2_input_sha256']=hashlib.sha256(c2_path.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,(int,float,Q)):return 1
 if isinstance(x,str):
  try:Q(x);return 1
  except(ValueError,ZeroDivisionError):return 0
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 if isinstance(x,dict):return sum(map(count,x.values()))
 return 0
source_bound=count(c)+13*len(c['operator_records_originalID_op_args'])+20000
assert source_bound<150000
result['source_whole_numeric_envelope_including_boundcache']=source_bound
result['hard_source_numeric']=150000
result['future_VALUE_static_cap']=75000
result['future_VALUE_hard_numeric']=100000
result['future_RAW_total_cap']=1000000
outdir=root/'reports/AMBIENT_B18_C3_RIGHT2_WHOLE_VALUE_GATE_20261105_v63'
outdir.mkdir(exist_ok=False)
report=outdir/'report.json';report.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result));signal.alarm(0)
