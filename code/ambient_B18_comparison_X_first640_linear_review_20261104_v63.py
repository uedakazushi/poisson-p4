#!/usr/bin/env python3
"""Independent sparse dual algebra and actual first640 source gate."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,resource,time

ROOT=Path(__file__).resolve().parents[1]
AUTHOR='AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63'
STEM=AUTHOR+'_REVIEW'

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def const(q):return{0:F(q)}if q else{}
def plus(a,b,scale=1):
 c=dict(a)
 for i,q in b.items():
  x=c.get(i,F(0))+scale*q
  if x:c[i]=x
  else:c.pop(i,None)
 return c
def times(a,b):
 a0=a.get(0,F(0));b0=b.get(0,F(0));out={}
 if a0*b0:out[0]=a0*b0
 for i in (set(a)|set(b))-{0}:
  x=a0*b.get(i,F(0))+b0*a.get(i,F(0))
  if x:out[i]=x
 return out
def full_to_bar(p):
 assert len(p)==10
 out={i:F(x)for i,x in enumerate(p[:8])if F(x)}
 r=F(p[8])+F(p[9])
 if r:out[8]=r
 return out
def unpack(p):return{i:F(q)for i,q in p if F(q)}
def encode(p):return[[i,str(q)]for i,q in sorted(p.items())]

def main():
 t=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_CPU,(60,60))
 out=ROOT/'reports'/STEM;assert not(out/'report.json').exists();out.mkdir(parents=True,exist_ok=True)
 rp='reports/'+AUTHOR+'/report.json';R=read(rp);assert R['all_passed'];hashes={rp:sha(rp)}
 for p,h in R['input_sha256'].items():assert sha(p)==h,p;hashes[p]=h
 for name,h in R['evidence_sha256'].items():
  p='reports/'+AUTHOR+'/'+name;assert sha(p)==h,p;hashes[p]=h
 for p,h in[('code/ambient_B18_comparison_X_first640_linear_20261104_v63.py',R['code_sha256']),('research/notes/'+AUTHOR+'.md',R['note_sha256'])]:assert sha(p)==h;hashes[p]=h
 peer_stems=['AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63_REVIEW',
  'AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63_REVIEW',
  'AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW',
  'AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_20261104_v63_REVIEW',
  'AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63_REVIEW']
 for st in peer_stems:
  p='reports/'+st+'/report.json';d=read(p);assert d['all_passed'];hashes[p]=sha(p)
  for name in ['input_sha256','checked_dependency_sha256']:
   for a,h in d.get(name,{}).items():
    if a in R['input_sha256']:assert R['input_sha256'][a]==h,(st,a)
 def get(st,f):return read('reports/'+st+'_20261104_v63/'+f)
 B=get('AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY','actual-UP-QB-640-typed-operator.json')
 C=get('AMBIENT_B18_T0_HONEST_C0','honest-C0-D0-D4-DAG.json')
 J=get('AMBIENT_B18_T3_FULL9_FIRSTJET','full9-firstjet-certificate.json')
 Q=get('AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC','quotient-firsts-and-actual-D0-firstjets.json')
 S=get('AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION','actual-source105-Artin20-values.json')
 H=get('AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20','actual-special-fiber20-coefficients.json')
 X0=get('AMBIENT_B18_COMPARISON_X_ORIGIN','actual-X0-and-unit-cancellation.json')
 atom_saved=get('AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR','phase1-actual-atom-firsts.json')
 expected=get('AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR','actual-first640-X-mod-m2.json')
 ops=B['operator_scalar_nodes'];roots=[node for col in B['B48_by640_sparse_columns']for row,node in col]
 needed=set();pending=roots[:]
 while pending:
  i=pending.pop()
  if i in needed:continue
  needed.add(i);n=ops[i]
  if n[0]=='sum':pending.extend(j for j,a,b in n[1])
  elif n[0]=='product':pending.extend(n[3])
  else:assert n[0]in['c','actual_C0_scalar_function']
 atoms=sorted({ops[i][1]for i in needed if ops[i][0]=='actual_C0_scalar_function'})
 assert(len(roots),len(needed),len(atoms))==(2808,1302,79)
 assert 9*(len(roots)+len(needed)+len(atoms))<150000
 (out/'preflight.json').write_text(json.dumps({'typed_ancestors':len(needed),'actual_scalar_atoms':len(atoms),'B_roots':len(roots),
  'no_native128_or_sourcejet_query':True,'input_sha256':hashes},indent=2)+'\n')
 vals={int(i):full_to_bar(p)for key in['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']for i,p in J[key].items()}
 basis=H['quotient_basis_monomials'];assert basis==S['A20_basis']
 def artin_first(raw):
  q={i:F(c)for i,c in raw};result=const(q.get(0,F(0)))
  for j in range(8):
   m=[0]*8;m[j]=1;c=q.get(basis.index(m),F(0))
   if c:result[j+1]=c
  return result
 vals.update({i:artin_first(p)for i,p in S['all_requested_original_source_root_values']})
 for i,p in enumerate(S['actual_physical_R0_readback'],2):vals[i]=artin_first(p)
 for j,i in enumerate(C['actual_C0_basis8_scalar_nodes']):vals[i]={j+1:F(1)}if j<7 else{}
 # Exact ordered prefix quotient, rebuilt by selecting monomial SUPPORT,
 # not by evaluating a derivative or duplicating a square coefficient.
 prefixes=[];words=J['source2_monomials']
 for record in C['all35_regular_coordinate_divisions']:
  f=record['original_pure_CA_source_function'];source=J['both57_source2']['CA']['pure_function2jets'][str(f)]
  assert len(words)==len(source)==55
  for j,q in enumerate(record['all7_regular_coordinate_divisions']):
   jet={}
   for word,c in zip(words,source):
    c=F(c)
    if not c or j not in word or any(k>j and k<7 for k in word):continue
    remaining=list(word);remaining.remove(j);assert len(remaining)<=1
    index=0 if not remaining else(remaining[0]+1 if remaining[0]<7 else 8)
    jet=plus(jet,{index:c})
   assert jet.get(0,F(0))==F(q['origin']);vals[q['node']]=jet
   prefixes.append({'node':q['node'],'actual_source_function':f,'prefix_index':j,'first9_bar':encode(jet),
    'retained_z_indices':list(range(j+1)),'always_retained_r_s':True})
 assert prefixes==atom_saved['all35_ordered_CA_prefix_firsts']
 # Source4 was independently reconstructed in the accepted prior peer. Here
 # recheck EVERY coefficient of numerator4=g3*quotient1 and derive its full9
 # physical adapter instead of trusting the eight saved override vectors.
 mon=[tuple(m)for m in Q['monomials_bivariate_through4']];adapter=[[F(x)for x in row]for row in J['YB_to_YA_linear']]
 derived_grad=[];g_identity_coefficients=0
 for family in['CB','CA']:
  f=Q['families'][family];g={m:F(q)for m,q in zip(mon,f['g_all15_coefficients'])if F(q)}
  assert all(sum(m)>=3 for m in g)
  for numerator,linear in zip(f['four_numerator_all15_coefficients'],f['quotient_firsts4_by2']):
   assert all(F(q)==0 for m,q in zip(mon,numerator)if sum(m)<=3)
   result={};a,b=map(F,linear)
   for e,q in g.items():
    for axis,c in enumerate((a,b)):
     n=list(e);n[axis]+=1;n=tuple(n)
     if sum(n)<=4:result[n]=result.get(n,F(0))+q*c
   assert all(result.get(m,F(0))==F(q)for m,q in zip(mon,numerator));g_identity_coefficients+=15
   gradient=[F(0)]*9
   if family=='CB':gradient=[a*x+b*y for x,y in zip(adapter[7],adapter[8])]
   else:gradient[7:]=[a,b]
   derived_grad.append(gradient)
 assert derived_grad==[[F(q)for q in row]for row in Q['all8_regular_function_first_gradients_in_FULL9_CA_F9']]
 nodes=Q['original_regular_function_nodes_unmodified']['old_CB']+Q['original_regular_function_nodes_unmodified']['new_CA']
 for node,gradient in zip(nodes,derived_grad):vals[node]=full_to_bar([F(0)]+gradient)
 assert [[i,encode(vals[i])]for i in nodes]==atom_saved['all8_g_quotient_overrides']
 tail=next(r for r in B['accepted_actual_regular_function_identities']if r['kind']=='accepted_actual_CB_f63_tail')
 raw_f63=tail['function'];plane_gB=tail['gB']
 assert raw_f63==J['both57_source2']['CB']['f63'] and raw_f63!=plane_gB
 assert B['old_multiplier_last_is_f63_not_new_gA']
 assert plane_gB not in C['actual_C0_basis8_scalar_nodes'];vals[plane_gB]={}
 dag=C['DAG_scalar_nodes'];memo=dict(vals);initial=set(vals)
 def original(i):
  if i in memo:return memo[i]
  n=dag[i]
  if n[0]=='c':p=const(F(n[1],n[2]))
  elif n[0]=='n':p={j:-q for j,q in original(n[1]).items()}
  elif n[0]=='a':p=plus(original(n[1]),original(n[2]))
  elif n[0]=='m':
   a=original(n[1]);p=times(a,original(n[2]))if a else{}
  else:raise AssertionError(('uncached original regular source function',i,n))
  memo[i]=p;assert 9*len(memo)<150000;return p
 actual_atoms={i:original(i)for i in atoms}
 assert [[i,encode(v)]for i,v in sorted(actual_atoms.items())]==atom_saved['all79_atom_values']
 constants={i:F(c)for i,c in X0['source_function_constant_refs']}
 assert all(v.get(0,F(0))==constants[i]for i,v in actual_atoms.items())
 # Explicit old f63 check in its own source then the actual YA adapter.
 source=J['both57_source2']['CB']['pure_function2jets'][str(raw_f63)]
 first=[F(source[words.index([j])])for j in range(9)]
 transformed=[sum((first[i]*adapter[i][j]for i in range(9)),F(0))for j in range(9)]
 assert original(raw_f63)==full_to_bar([F(source[0])]+transformed)
 typed={}
 for i in sorted(needed):
  n=ops[i]
  if n[0]=='c':v=const(F(n[1],n[2]))
  elif n[0]=='actual_C0_scalar_function':v=actual_atoms[n[1]]
  elif n[0]=='sum':
   v={}
   for j,a,b in n[1]:v=plus(v,typed[j],F(a,b))
  elif n[0]=='product':
   v=const(F(n[1],n[2]))
   for j in n[3]:v=times(v,typed[j])
  else:raise AssertionError(n)
  typed[i]=v
 output=[];nonzero=0;P_checked=0
 saved_P=Q['all640_exact_D0_firstjet81_columns']
 P_peer=read('reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63_REVIEW/report.json')
 P_source='reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json'
 assert P_peer['source4_quotient_firsts_and_D0_firstspan38_independently_accepted']
 assert P_peer['input_sha256'][P_source]==R['input_sha256'][P_source]
 for j,col in enumerate(B['B48_by640_sparse_columns']):
  entries={row:typed[node]for row,node in col if typed[node]}
  assert len(saved_P[j])==81 and len(C['honest_D0_relation_columns'][j])==9
  for row in range(1,9):
   own=full_to_bar([F(0)]+[F(q)for q in saved_P[j][9*row:9*(row+1)]])
   P_checked+=1
   if own:entries[47+row]=own
  assert all(entries.get(i,{}).get(0,F(0))==F(X0['full_X0_56_by768'][i][j])for i in range(56))
  output.append([[i,encode(v)]for i,v in sorted(entries.items())]);nonzero+=sum(map(len,entries.values()))
 assert output==expected['source_columns640_rows56_sparse_dual9'] and nonzero==409
 active=sum(map(len,typed.values()))+sum(map(len,memo.values()))+nonzero
 assert active<150000
 proof={'all_passed':True,'source_atoms79_and_typed_ancestors1302_independently_rebuilt':True,
  'all35_prefix_quotients_from_monomial_support_no_factor2':True,
  'all8_g_derivatives_from_source4_membership_and_same_full9_adapter':True,
  'source4_numerator_identity_coefficients_checked':g_identity_coefficients,
  'raw_CB_f63_node_distinct_from_gB_plane_preserved':{'raw_f63':raw_f63,'gB_plane':plane_gB},
  'source105_H20_first_coefficient_extraction_same_CA_F9':True,
  'all8x640_P_derivatives_from_accepted_exact_source_cache_and_s_equals_r':P_checked,
  'Pcache_exact_prior_independent_source_peer_hash_bound':True,
  'fresh_P_original_DAG_reconstruction_claimed':False,
  'all35840_constants_match_accepted_actual_X0':True,
  'all409_output_coefficients_of_actual_first640_X1_match':True,
  'sparse_live_coefficient_count':active,'dual_AD_not_source_implicit_solve':True,
  'no_native128_X1_newjet_point_prime_GB_Gaussian_or_T_length_depth':True}
 (out/'independent-prefix-source4-and-first640-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 for p in[str(Path(__file__).relative_to(ROOT)),'research/notes/'+STEM+'.md']:hashes[p]=sha(p)
 report={'id':STEM,'all_passed':True,'frozen':True,
  'actual_first640_X1_source_bound_independently':True,
  'all79_atoms_1302_nodes_35prefixes_8source4_gradients_checked':True,
  'all640_P8_derivatives_accepted_source_cache_hash_bound':True,
  'actual_nonzero_output_coefficients':409,
  'native128_X1_T_mu_length_K_torsion_or_depth_computed':False,
  'input_sha256':hashes,'evidence_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in out.iterdir()if p.is_file()},
  'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}),flush=True)

if __name__=='__main__':main()
