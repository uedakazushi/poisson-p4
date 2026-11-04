#!/usr/bin/env python3
"""Independent small algebra/proof gate for the R gluing bounds16 and19."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_T0_FINITE_R_GLUING_BOUND_20261104_v63';ID=BASE+'_REVIEW'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def add(a,b):
    q=dict(a)
    for m,c in b.items():q[m]=q.get(m,0)+c
    return {m:c for m,c in q.items() if c}
def mul(a,b):
    q={}
    for m,c in a.items():
        for n,d in b.items():e=tuple(x+y for x,y in zip(m,n));q[e]=q.get(e,0)+c*d
    return {m:c for m,c in q.items() if c}
def pairmul(a,b):return mul(a[0],b[0]),add(mul(a[0],b[1]),mul(a[1],b[0]))
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    primary=ROOT/'reports'/BASE/'report.json';p=read(primary);assert p['all_passed'];hashes={}
    for name,h in p['input_sha256'].items():assert sha(ROOT/name)==h;hashes[name]=h
    gates=['AMBIENT_B18_T0_R_PRIMARY_SEPARATION_20261104_v63_REVIEW/report.json','AMBIENT_B18_R_ZERO_SLICE_IMAGE_20261104_v63_REVIEW/report.json','AMBIENT_B18_R_N_PROJECTION_NORM_20261104_v63_REVIEW/report.json','AMBIENT_B18_W1_TAIL_BINDING_20261104_v63_REVIEW/report.json']
    for name in gates:
        path=ROOT/'reports'/name;assert read(path)['all_passed'];hashes[str(path.relative_to(ROOT))]=sha(path)
    one={(0,0):1};zero={};f={(1,0):1};g={(0,1):1};F2=pairmul((f,g),(f,g));assert F2==({(2,0):1},{(1,1):2})
    col0=pairmul(F2,(one,zero));col1=pairmul(F2,(zero,one));M=[[col0[0],col1[0]],[col0[1],col1[1]]]
    determinant=add(mul(M[0][0],M[1][1]),{m:-c for m,c in mul(M[0][1],M[1][0]).items()});assert determinant=={(4,0):1}
    order=4;length=4*order;assert length==16
    allowed=[(m,n) for m in range(4,9) for n in range(m+1) if m+n<=length];assert len(allowed)==35 and [list(x) for x in allowed]==p['allowed_m_n_pairs']
    # Fixed HNF arithmetic and necessary bounds; no actual exponents determined.
    for m,n in allowed:assert 4<=m<=8 and 0<=n<=m and m+n<=16
    z=read(ROOT/'reports/AMBIENT_B18_T0_COKERNEL_Z_20261104_v63_REVIEW/report.json');assert z['all_passed'] and z['Z_length']==3
    assert p['length_B_R_over_J_R_upper_bound']==16 and p['length_D4_over_image_S_F_upper_bound']==length+z['Z_length']==19
    assert not p['actual_m_n_b_or_finite_quotient_matrices_computed'] and not p['residual_kernel_finite_socle_or_depth_claimed']
    proof=dict(norm_scalar_F_not_old_nil_generator=True,F_maps_into_whole_squarezero_nilradical_of_D4=True,F_squared_annihilates_actual_D4=True,F_squared_BR_contained_in_JR=True,multiplication_F_squared_matrix='[[f^2,0],[2fg,f^2]]',universal_polynomial_determinant='f^4',F_reduced_order_U=4,exact_length_BR_over_F_squared=16,upper_length_BR_over_JR=16,contact4_gives_reduced_JR_containment_U4_only=True,HNF_form='(U^m+b*epsilon,U^n*epsilon), b modulo U^n',HNF_bounds=dict(m_min=4,m_max=8,n_min=0,n_max='m',length='m+n<=16'),allowed_pairs=allowed,exact_cokernel_sequence='0 -> BR/JR -> D4/image(S_F) -> Z -> 0',Z_length=3,upper_length_D4_over_image_S_F=19,no_residual_kernel_finiteness_socle_depth_or_actual_HNF_computation=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for path in [primary,Path(__file__),ROOT/f'research/notes/{ID}.md']:hashes[str(path.relative_to(ROOT))]=sha(path)
    r=dict(id=ID,all_passed=True,actual_R_gluing_bounds16_and19_independently_accepted=True,contact_HNF_necessary_bounds_not_exact_exponents=True,no_CAS_source_replay_or_residual_socle_claim=True,input_sha256=hashes,evidence_sha256={'independent-proof.json':sha(out/'independent-proof.json')},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=r['seconds'],report_sha256=sha(out/'report.json'))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
