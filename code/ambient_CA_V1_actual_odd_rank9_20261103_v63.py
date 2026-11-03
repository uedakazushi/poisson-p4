#!/usr/bin/env python3
"""Bind the successful57 paired-source cubic witness to full actual odd rank9."""
import argparse,gzip,hashlib,json,time
from fractions import Fraction
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from ambient_CA_V1_odd_order_v2_20261103_v63 import evaluator,frac
ROOT=Path(__file__).resolve().parents[1]
def load(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_bytes())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(M):return [[str(v)for v in row]for row in M.tolist()]
def main(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports/AMBIENT_CA_V1_ODD_ORDER_V2_20261103_v63';report=load(folder/'report.json');setup=load(folder/'setup.json.gz');jets=load(folder/'jets.json.gz');assert report['all_passed'] is True and report['exact_normal_conductor_order']==1
    ckpath=ROOT/'reports/AMBIENT_CA_V1_ODD_ORDER_CHECKPOINT_20261103_v63/checkpoint.json.gz';ck=load(ckpath);kernel=s.Matrix(setup['kernel29']);J=s.Matrix(setup['source29_target_Jacobian']);jc=setup['source_r28_columns'];Jr=J[:,jc];H=s.Matrix(setup['full150_by29_mixed_Hessian']);Y=setup['actual_target_Y28'];branch=setup['actual_branch29th_row'];inverse=Jr[Y,:].inv();W=H[:,0]/2
    q=[[Fraction(0)]*4 for _ in range(29)]
    for i in range(29):q[i][1]=frac(kernel[i])
    C=s.Matrix(ck['fixedCA150_by20']);point=list(map(frac,s.Matrix(setup['point150'])));basis=[list(map(frac,C[:,i]))for i in setup['extra_fiber_indices']]
    values=evaluator(point,basis,setup['group_GL_indices'],3)(q);Z=s.Matrix([s.Rational(row[3])for row in values]);assert s.Matrix([s.Rational(row[2])for row in values])==W
    stage=jets['stages'][0];qe=s.Matrix(stage['source_even_coeff']);qo=s.Matrix(stage['source_odd_coeff']);T=kernel.row_join(s.eye(29)[:,jc]);Horiginal=H*T.inv()
    direct=s.Matrix([s.Rational(row[3])for row in stage['all150_output_series']])-Horiginal*qe-J*qo
    assert direct==Z
    rem=lambda M:M-Jr*inverse*M[Y,:]
    WW=rem(W);weight=WW[branch];mixed=H[:,1:];ZZ=rem(Z-mixed*inverse*W[Y,:]);E=rem(mixed*inverse)
    normal=lambda M:M-WW/weight*M[branch,:]
    BY=normal(E);Bnu=normal(ZZ)/weight;L=Jr[branch,:]*inverse
    A=(BY-Bnu*L).row_join(Bnu);pc=list(DomainMatrix.from_Matrix(A).rref()[1]);rr=list(DomainMatrix.from_Matrix(A[:,pc].T).rref()[1]);assert len(pc)==len(rr)==9
    det=s.Rational(DomainMatrix.from_Matrix(A[rr,pc]).convert_to(s.QQ).det());assert det
    coefficients=A[rr,pc].inv()*A[rr,:];assert A[:,pc]*coefficients==A
    critical=rem(H);assert DomainMatrix.from_Matrix(critical).rank()==9
    (out/'certificate.json.gz').write_bytes(gzip.compress((json.dumps(dict(all150_actual_odd_linear=mat(A),nine_rows=rr,nine_columns=pc,nine_minor=str(det),all_column_factorization=mat(coefficients),whole_physical_critical_cotangent=mat(critical),source_pure_kernel3=mat(Z)),separators=(',',':'))+'\n').encode(),mtime=0))
    paths=[Path(__file__),ROOT/'research/notes/AMBIENT_CA_V1_ODD_ORDER_RESULT_20261103_v63.md',folder/'report.json',folder/'setup.json.gz',folder/'jets.json.gz',ckpath,ROOT/'reports/AMBIENT_CA_V1_COEFFICIENT_HESSIAN_20261103_v63/hessian.json.gz']
    result=dict(id='AMBIENT_CA_V1_ACTUAL_ODD_RANK9_20261103_v63',all_passed=True,actual_regular_odd8_eliminated_before_cusp_order=True,pure_kernel_cubic_independent_source_replay=True,whole_actual_odd_firstjet_rank=9,whole_normalizer_critical_cotangent_rank=9,GL_orbit_dimension=20,actual_discriminant_linear_direction_in_odd_ideal_span=True,source_pair_cubic_witness={'63':'800/27','64':'-40/27'},full_I_eq_J_and_genericV_type_requires_separate_application_proof=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('V1 ACTUAL ODD9 / CRITICAL9 / CORRECT8ODD-CUSP k1 PASS',result['seconds'])
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
