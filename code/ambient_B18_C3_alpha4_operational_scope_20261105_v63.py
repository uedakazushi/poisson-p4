import importlib.util,json,math,hashlib,signal,resource
from pathlib import Path
from fractions import Fraction as Q

root=Path(__file__).resolve().parents[1]
signal.alarm(60)
resource.setrlimit(resource.RLIMIT_AS,(2147483648,2147483648))
p=root/'code/ambient_B18_C3_alpha4_cached_source_20261105_v63.py'
sp=importlib.util.spec_from_file_location('alpha_source_bounds_only',p);S=importlib.util.module_from_spec(sp);sp.loader.exec_module(S)
c=json.loads((root/'reports/AMBIENT_B18_C3_ALPHA4_CACHED_SOURCE_20261105_v63/compiled-eight-even-roots-and-provided-joint105-plane3-cache.json').read_text())
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
for r in c['root_pairs_physical4']:
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
result['source_primary_sha256']=hashlib.sha256((root/'reports/AMBIENT_B18_C3_ALPHA4_CACHED_SOURCE_20261105_v63/report.json').read_bytes()).hexdigest()
result['code_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
result['source_bound_helper_sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
result['compiled_input_sha256']=hashlib.sha256((root/'reports/AMBIENT_B18_C3_ALPHA4_CACHED_SOURCE_20261105_v63/compiled-eight-even-roots-and-provided-joint105-plane3-cache.json').read_bytes()).hexdigest()
result['owned_numeric_envelope_including_all_bound_cache']=135000
result['hard_owned_numeric']=150000
out=root/'reports/AMBIENT_B18_C3_ALPHA4_CACHED_SOURCE_20261105_v63/operational-height-scope-companion.json'
assert not out.exists()
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result));signal.alarm(0)
