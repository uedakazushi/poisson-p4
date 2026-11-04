#!/usr/bin/env python3
"""Constant-matrix/metadata bounds only: NO new invariant coefficient or rank job."""
import argparse,gzip,hashlib,itertools,json,math,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from weight1_structure import sp_basis
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_DEGREE24_THETA_H_REPAIR_PREFLIGHT_20261104_v63'
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def M(a):return [[Q(str(z))for z in r]for r in a]
def absM(a):return [[abs(x)for x in r]for r in a]
def T(a):return list(map(list,zip(*a)))
def mm(a,b):return [[sum((a[i][k]*b[k][j]for k in range(len(b))),Q())for j in range(len(b[0]))]for i in range(len(a))]
def add(a,b):return [[x+y for x,y in zip(r,t)]for r,t in zip(a,b)]
def denom(a):return math.lcm(*(x.denominator for r in a for x in r))
def necklaces(n):return [w for w in itertools.product(range(2),repeat=n)if w==min(w[i:]+w[:i]for i in range(n))]
def main(out):
 started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 invpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
 old12path=ROOT/'reports/Q_NEW_FAMILY_DEGREE12_14_INVARIANT_EXCLUSION_20261103_v63/certificate.json'
 old16path=ROOT/'reports/Q_NEW_FAMILY_DEGREE16_INVARIANT_PROBE_20261103_v63/report.json'
 prepath=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_PREFLIGHT_20261104_v63/report.json'
 failpath=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_UNIVERSAL_BASIS_20261104_v63/report.json'
 betapath=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63/native-input.txt.gz'
 data,old12,old16,pre,failed=map(read,[invpath,old12path,old16path,prepath,failpath])
 assert failed['status']=='degree12_rank_gate_failed'and failed['ranks']['12']==8
 assert old12['full_degree12_dimension']==9 and len(old12['degree12_basis_columns'])==9
 # Read only coefficient norms of ALREADY FROZEN integer beta numerators b=1080BH beta.
 tokens=iter(gzip.decompress(betapath.read_bytes()).decode().split());assert [int(next(tokens))for _ in range(4)]==[4,16,136,1080]
 norms=[]
 for _ in range(16):
  n=int(next(tokens));norm=0
  for k in range(n):norm+=abs(int(next(tokens)));[next(tokens)for _ in range(4)]
  norms.append(norm)
 O=M(data['symplectic_form16']);Oi=M(s.Matrix(O).inv().tolist());rep=[M(a)for a in data['all_ten_representation_matrices']]
 gs=sp_basis();G=M(s.Matrix([[s.trace(a*b)for b in gs]for a in gs]).inv().tolist())
 # Positive rational majorants for quadratic theta/H evaluated at coefficient l1 norms.
 V=[[mm(mm(absM(O),absM(r)),[[Q(x)]for x in norms])[i][0]for r in rep]for i in range(16)]
 mu=[[sum((norms[i]*V[i][j]/2 for i in range(16)),Q())]for j in range(10)]
 co=mm(absM(G),mu)
 theta=[[sum((co[k][0]*abs(rep[k][i][j])for k in range(10)),Q())for j in range(16)]for i in range(16)]
 H=add(theta,mm(mm(mm(absM(Oi),V),absM(G)),T(V)))
 # A sufficient common coefficient denominator, obtained from the unexpanded formulas.
 dr=math.lcm(*(denom(a)for a in rep));do,doi,dg=denom(O),denom(Oi),denom(G)
 q=math.lcm(2*dg*do*dr*dr,doi*do*do*dr*dr*dg)
 A=[[[int(q*x)for x in row]for row in a]for a in[theta,H]]
 assert all(q*x==int(q*x)for a in[theta,H]for row in a for x in row)
 bounds=[]
 for n in[6,8]:
  values=[]
  for word in necklaces(n):
   a=[[int(i==j)for j in range(16)]for i in range(16)]
   for letter in word:a=mm(a,A[letter])
   z=sum(a[i][i]for i in range(16));assert z.denominator==1
   values.append(int(z))
  maximum=max(values);support=pre['source_low_word_output_bounds'][0 if n==6 else 1]['rigorous_monomial_bound_after_guaranteed_B_cancel']
  bounds.append(dict(trace_length=n,beta_degree=2*n,necklace_count=len(values),integer_matrix_scale=q,
   integer_trace_denominator=f'{q}^{n}*(1080*B*H)^{2*n}',
   coefficient_l1_upper_bit_bounds=[v.bit_length()for v in values],coefficient_decimal_digits_upper=max(len(str(v))for v in values),
   output_support_bound=support,compact_JSON_bytes_upper=support*(len(str(maximum))+40),
   guaranteed_B_denominator_power_at_most=n,remaining_H_denominator_power=2*n))
 cd=pre['source_low_word_output_bounds'];mixed_cd16=2*cd[0]['compact_JSON_output_byte_bound']+bounds[0]['compact_JSON_bytes_upper']+cd[1]['compact_JSON_output_byte_bound']
 mixed_old16=2*cd[0]['compact_JSON_output_byte_bound']+bounds[0]['compact_JSON_bytes_upper']+bounds[1]['compact_JSON_bytes_upper']
 words16=[old16['invariant_definitions'][i]for i in old16['ambient_independent_columns']]
 assert len(words16)==20 and all(d['kind']=='trace8'for d in words16)
 report=dict(id=STEM,all_passed=True,preflight_only_no_new_invariant_or_rank_computation=True,
  precise_representation_repair='theta=rho_W(Gram^-1 mu), H=theta+Omega^-1 V Gram^-1 V^T are 16x16; H is not an sp4 four-dimensional-representation matrix',
  old12_full_basis_indices=old12['degree12_basis_columns'],old12_full_basis_words=[old12['degree12_trace_words'][i]for i in old12['degree12_basis_columns']],
  old16_full_basis_definitions=words16,same_formal_plane_from_failed_packet=True,
  reuse_failed_packet_CD_coefficients=True,no_new_actual_NEW25_point_or_alternate_plane=True,
  corrected_plan=['append old theta/H trace6 columns to saved rank8 CD space until rank9 (prefer exactly one old word)',
   'form F4*Inv12 + Sym2 Inv8 and require rank19',
   'first try saved36 CD trace8 coefficient columns for rank20',
   'only if preauthorized and CD trace8 fails, try frozen20 old theta/H trace8 basis words',
   'form70 products and require rank67; choose67 columns with coefficient minor'],
  formal_poly_supports_by_beta_degree={str(d):(d+1)*(d+2)//2 for d in[2,4,6,8,12,16,24]},
  formal_theta_H_build_polynomial_products_upper=2720,
  formal_all14_old_trace6_polynomial_products_upper=12*16**3+14*16**2,
  optional_all36_old_trace8_additional_products_upper=16*16**3+36*16**2,
  formal_theta_H_all_prefix_matrices_through4=30,
  formal30_matrices_full_coefficient_word8_storage_residues_upper=30*16**2*45,
  proposed_corrected_job_hard_limits=dict(cpu_and_wall_seconds=60,address_space_bytes=2_000_000_000,certificate_bytes=32*1024*1024),
  future_source_old_word_bounds=bounds,future_mixed_source2_CD_trace6_plus1_old_trace6_plus1_CD_trace8_bytes_upper=mixed_cd16,
  future_mixed_source2_CD_trace6_plus1_old_trace6_plus1_old_trace8_bytes_upper=mixed_old16,
  future_source_raw32MiB_bound_proven=mixed_old16<32*1024*1024,
  source_expansion_runtime_or_memory_success_not_guaranteed=True,
  future_source_word_at_a_time_liveness='at most two half-word matrices plus temporary product; do not build all16x16 prefix words on source',
  no_source_expansion_authorized_by_this_preflight=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),invpath,old12path,old16path,prepath,failpath,betapath]},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:report[k]for k in['all_passed','formal_all14_old_trace6_polynomial_products_upper','future_source_old_word_bounds','future_mixed_source2_CD_trace6_plus1_old_trace6_plus1_CD_trace8_bytes_upper','future_mixed_source2_CD_trace6_plus1_old_trace6_plus1_old_trace8_bytes_upper','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
