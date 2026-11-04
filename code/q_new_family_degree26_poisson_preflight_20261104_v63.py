#!/usr/bin/env python3
"""Metadata, constant majorants and graded candidate counts ONLY; no Inv26 job."""
import argparse,gzip,hashlib,itertools,json,math,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_DEGREE26_POISSON_PREFLIGHT_20261104_v63'
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mm(a,b):return [[sum(a[i][k]*b[k][j]for k in range(len(b)))for j in range(len(b[0]))]for i in range(len(a))]
def odd_support(d):
 cancel=d//2;ans=0
 for a in range(d+1):
  for Y in range(2*d+1):
   if not 0<=2*Y+1-a<=3*d:continue
   for W in range(3*d+1):
    for B in range(d-cancel+1):
     if 2*Y+1+W+B<=4*d-cancel and 2*(a+Y+W+B)+1<=9*d-2*cancel:ans+=1
 return ans
def main(out):
 started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 spinpath=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
 invpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
 bpath=ROOT/'reports/Q_NEW_FAMILY_HESSIAN_LOW_INVARIANT_EXPANSION_20261104_v63/native-input.txt.gz'
 basispath=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63_REVIEW/report.json'
 low14path=ROOT/'reports/Q_NEW_FAMILY_DEGREE12_14_INVARIANT_EXCLUSION_20261103_v63/certificate.json'
 charpath=ROOT/'work/old-plane-v62/invariant-character-through32.json'
 if not charpath.is_file():charpath=ROOT/'reports'/STEM/'character-records-through26.json'
 spin,inv,basis,low14,chars=map(read,[spinpath,invpath,basispath,low14path,charpath]);assert basis['all_passed']and low14['full_degree14_dimension']==3
 dimensions={r['degree']:r['invariant_dimension']for r in chars}
 assert [dimensions[d]for d in[12,14,18,24,25,26]]==[9,3,9,67,0,45]
 for r in chars[:27]:assert sum(t['sign']*t['multiplicity']for t in r['Weyl_terms'])==r['invariant_dimension']
 (out/'character-records-through26.json').write_text(json.dumps(chars[:27],indent=2)+'\n')
 O=s.Matrix(inv['symplectic_form16']);Oi=O.inv();q=math.lcm(*(int(s.denom(z))for z in Oi));assert q==5
 Oih=[[abs(int(q*Oi[i,j]))for j in range(16)]for i in range(16)];nnz=sum(z!=0 for r in Oih for z in r)
 tokens=iter(gzip.decompress(bpath.read_bytes()).decode().split());assert [int(next(tokens))for _ in range(4)]==[4,16,136,1080];norm=[]
 for _ in range(16):
  n=int(next(tokens));v=0
  for _ in range(n):v+=abs(int(next(tokens)));[next(tokens)for _ in range(4)]
  norm.append(v)
 mon=spin['quadratic_monomials'];norm_mats=[];derivatives=[]
 for key in['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136']:
  coefficients=[[4*Q(str(z))for z in row]for row in spin[key]];assert all(z.denominator==1 for row in coefficients for z in row)
  vals=[sum(abs(int(c))*norm[i]*norm[j]for c,(i,j)in zip(row,mon))for row in coefficients];norm_mats.append([vals[i:i+4]for i in range(0,16,4)])
  d=[]
  for k in range(16):
   vals=[]
   for row in coefficients:
    v=0
    for c,(i,j)in zip(row,mon):
     if i==k:v+=abs(int(c))*norm[j]
     if j==k:v+=abs(int(c))*norm[i]
    vals.append(v)
   d.append([vals[i:i+4]for i in range(0,16,4)])
  derivatives.append(d)
 w4=[(0,0,0,0),(0,0,0,1),(0,0,1,1),(0,1,0,1),(0,1,1,1)]
 w6=[(0,0,0,1,0,1),(0,0,0,1,1,1),(0,1,1,0,1,1)]
 prefixes=set();gradient_norms=[]
 for w in w4+w6:
  grad=[0]*16
  for i,letter in enumerate(w):
   rest=w[i+1:]+w[:i]
   for k in range(1,len(rest)+1):prefixes.add(rest[:k])
   mat=[[int(a==b)for b in range(4)]for a in range(4)]
   for v in rest:mat=mm(mat,norm_mats[v])
   for k in range(16):grad[k]+=sum(mat[a][b]*derivatives[letter][k][b][a]for a in range(4)for b in range(4))
  gradient_norms.append(grad)
 records=[]
 for degree,pairs in[(14,list(itertools.combinations(range(5),2))),(18,[(i,5+j)for i in range(5)for j in range(3)])]:
  nums=[sum(Oih[i][j]*gradient_norms[a][i]*gradient_norms[b][j]for i in range(16)for j in range(16))for a,b in pairs]
  support=odd_support(degree);digits=max(len(str(z))for z in nums)
  records.append(dict(beta_degree=degree,candidate_count=len(pairs),ambient_gradient_pair_indices=pairs,
   numerator_l1_bit_bounds=[z.bit_length()for z in nums],max_decimal_digits=digits,source_normalized_odd_support_bound=support,
   compact_JSON_bytes_per_output_upper=support*(digits+40),
   integer_denominator_before_content=f'5*4^{degree//2+1}*1080^{degree}',B_denominator_power=degree//2,H_denominator_power=degree,
   c_homogeneity_power=degree//2,source_common_X_removed_after_u_equals_aX=True))
 prefix_products=sum(len(w)>=2 for w in prefixes)*64
 gradient_products=sum(len(set(w))*16*16 for w in w4+w6)
 worst_new_maps=3*records[0]['compact_JSON_bytes_per_output_upper']+6*records[1]['compact_JSON_bytes_per_output_upper']
 report=dict(id=STEM,all_passed=True,preflight_only_no_invariant_or_gradient_polynomial_or_rank_job=True,
  character_dimensions={str(d):dimensions[d]for d in[12,14,18,24,25,26]},same_formal_ambient_plane_required=True,
  ambient_derivative_order='derive in ALL16 beta coordinates from136 quadratic C/D coefficients BEFORE substituting the2plane or source; never use2plane derivatives',
  octic_basis_CD_trace_words=w4,CD_degree12_trace_words=w6,
  N12_replacement_identity='N12=T011011-F4*T0011/20+F4*p4/5',
  F4_Poisson_central_on_invariants='Hamiltonian vector of F4 is a moment-dependent linear combination of the10 physical infinitesimal Sp4 actions',
  Inv14_known_complete_span='all10 Poisson brackets in any fullInv8 basis span the already-proved3D Inv14; require rank3 on SAMEplane',
  intended_rank_gates=[{'Inv14':3},{'F4_Inv14':3},{'Inv18':9},{'Inv26':45}],
  degree18_candidate_count=18,degree18_candidates='3F4*Inv14 +15 brackets of5CDtrace4 with3CDtrace6',
  degree26_candidate_count=57,degree26_candidates='3Inv14*full9Inv12=27 plus6primitive18*5Inv8=30',
  no_assumed_trace_generation_or_minimal_generator_number26=True,
  formal_coefficient_rows={str(d):(d+1)*(d+2)//2 for d in[14,18,26]},
  gradient_prefix_word_counts={str(n):sum(len(w)==n for w in prefixes)for n in range(1,6)},
  gradient_prefix_matrix_polynomial_products_upper=prefix_products,ambient_gradient_contraction_products_upper=gradient_products,
  inverse_Omega_nonzero_entries=nnz,Poisson_scalar_products_upper=25*nnz,
  formal_all8_ambient_gradients_storage_residues_upper=16*(5*36+3*78),formal_57_degree26_matrix_entries=57*378,
  future_source_bounds=records,future_source9_output_bytes_upper=worst_new_maps,
  future_source9_output_raw32MiB_completion_bound=worst_new_maps<=32*1024*1024,
  future_degree26_column_support_upper=odd_support(26),
  exact_degree26_source_denominator='B^13*H^26 times known integer; strip commonX ONLY afteru=aX, thenY=X^2',
  exact_source_parity='c^13 and finite tau imply odd(u,X); invariantTdegree even underT^2=B; afterdividingX resultingshapeparity even',
  degree14_and18_parity='c^7,c^9 are odd; normalizedmaps I/(c^(d/2)*X) are rational in(a=u/X,Y=X^2,W,B)',
  source_rank_failure_no_automatic_H16_gradient_fallback=True,
  H16_derivative_cost_if_later_needed='fulltheta5H ambientgradient would require16-variable differentiation ofHessian model/vector recurrence; not part of this proposal or budget',
  proposed_universal_job_limits=dict(cpu_and_wall_seconds=60,address_space_bytes=2_000_000_000,certificate_bytes=32*1024*1024),
  proposed_future_source_limits_not_yet_allocated=dict(cpu_and_wall_seconds=60,address_space_bytes=2_000_000_000,output_terms=100000,intermediate_terms=150000,total_raw_bytes=32*1024*1024),
  no_NEW_source_point_or_degree26_job_started=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),spinpath,invpath,bpath,basispath,low14path,charpath]},
  evidence_sha256={'character-records-through26.json':sha(out/'character-records-through26.json')},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:report[k]for k in['all_passed','gradient_prefix_word_counts','gradient_prefix_matrix_polynomial_products_upper','ambient_gradient_contraction_products_upper','inverse_Omega_nonzero_entries','future_source_bounds','future_source9_output_bytes_upper','future_degree26_column_support_upper','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
