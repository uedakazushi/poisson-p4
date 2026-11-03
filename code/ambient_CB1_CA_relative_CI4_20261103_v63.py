#!/usr/bin/env python3
"""Publish/audit completed QQ(b) physical computation without repeating its98s run."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CB1_CA_RELATIVE_CI4_20261103_v63'
def read(path):
    with(gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as f:return json.load(f)
def sha(path):
    h=hashlib.sha256()
    with path.open('rb')as f:
        for q in iter(lambda:f.read(1<<20),b''):h.update(q)
    return h.hexdigest()
def main(out,recompute):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=True);assert not(out/'report.json').exists()
    data=read(out/'generic-quadratic.json.gz');mapping=read(out/'quotient-map.json')
    assert data['all_passed']and mapping['all_passed']and data['quadratic_obstruction_rank']==2
    workers=[ROOT/'code/ambient_CB1_CA_generic_quadratic_worker_20261103_v63.py',ROOT/'code/ambient_CB1_CA_quotient_map_worker_20261103_v63.py']
    if recompute:
        (ROOT/'work/v62-review').mkdir(parents=True,exist_ok=True)
        for p in workers:subprocess.run([str(ROOT/'.venv/bin/python'),str(p)],check=True,timeout=120)
        fresh=read(ROOT/'work/v62-review/CB1-CA-two-normal-quadratic-generic-probe.json')
        for key in['frame58_by58','implicit36_inverse','normal58_by2','raw150_by3','corrected150_by3']:
            assert fresh[key]==data[key],key
        freshmap=read(ROOT/'work/v62-review/CB1-CA-generic-quotient-map-probe.json')
        assert freshmap['normal2']==mapping['normal2']and freshmap['eta_square_zero_in_generic_CI4']
    b,r,z=s.symbols('b r z');field=s.QQ.frac_field(b)
    parse=lambda q:s.sympify(q,locals={'b':b})
    mat=lambda rows:s.Matrix([[parse(q)for q in row]for row in rows])
    R=mat(data['corrected150_by3']);assert R.shape==(150,3)
    assert len(DomainMatrix.from_Matrix(R).convert_to(field).rref()[1])==2
    rows=data['obstruction_basis_rows'];assert rows==[63,68]
    forms=[s.expand(R[i,0]*r*r+R[i,1]*r*z+R[i,2]*z*z)for i in rows]
    M=s.Matrix([[R[i,0],R[i,1],R[i,2],0]for i in rows]+[[0,R[i,0],R[i,1],R[i,2]]for i in rows])
    # This order differs by a row permutation from the saved witness.
    determinant=s.factor(M.det());assert determinant==-1024*b**4/(25*(b*b-3)**4)
    pointfolder=ROOT/'reports/AMBIENT_CB1_CA_TRANSVERSE_CI4_20261103_v63'
    point=read(pointfolder/'certificate.json');assert read(pointfolder/'report.json')['all_passed']
    for key in['frame58_by58','normal58_by2','implicit36_inverse','raw150_by3','corrected150_by3']:
        A=mat(data[key]).subs(b,1);B=s.Matrix(point[key]);assert(A-B).applyfunc(s.cancel)==s.zeros(*A.shape),key
    assert mat(data['frame58_by58']).subs(b,1).det()==-2
    c0,c1=map(parse,mapping['normal2']);assert s.expand(c0+b**3-b)==0 and s.expand(c1+2*b*b)==0
    eta=s.expand(c1*r-c0*z);v=b*r-z
    coeff1=9*b*b*(b*b-3)**2/8;coeff2=15*b*b*(b*b-3)**2/4
    assert s.cancel(eta*eta-coeff1*forms[0]-coeff2*forms[1])==0
    assert s.cancel(v*v-(b*b-3)**2*forms[0]/24)==0
    change=s.Matrix([[s.diff(eta,r),s.diff(eta,z)],[s.diff(v,r),s.diff(v,z)]]).det()
    assert s.factor(change)==-b*b*(b*b-3)
    assert s.factor(b*c0-c1)==-b*b*(b*b-3)
    assert mapping['eta_square_zero_in_generic_CI4']
    deps=[ROOT/'reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63/report.json',
          ROOT/'reports/AMBIENT_CB_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63_REVIEW/report.json',
          ROOT/'reports/AMBIENT_CB_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63_REVIEW/report.json']
    assert all(read(p)['all_passed']for p in deps)
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',out/'generic-quadratic.json.gz',out/'quotient-map.json',*workers,pointfolder/'report.json',pointfolder/'certificate.json',*deps]
    report=dict(id=STEM,all_passed=True,field='QQ(b), formal C base at b=1 and a nonempty open of S20',
        completed_full_generic_physical_worker_reused_without_repeating=True,
        completed_worker_seconds=data['seconds'],worker_full_recompute_requested=recompute,
        all150_generic_quadratic_obstruction_rank=2,generic_degree3_span_unit='1024*b^4/(25*(b^2-3)^4)',
        full_generic_frame_and_raw_corrected_data_specialize_exactly_to_closed_CI4=True,
        actual_relative_residuals_in_normal_square_by_constant_tangent22_along_smooth20_base=True,
        generic_rank2_descends_by_dominant_GL_times_b_and_constant_unit2minor=True,
        full_relative_ideal_is_two_quadrics_by_Nakayama_no_higher_jets_required=True,
        actual_CB1_intersection_CA_relative_free_rank=4,
        actual_completed_intersection='C[[s1..s20,u,v]]/(u^2,v^2)',
        CB_meeting_quotient='C[[s1..s20,v]]/(v^2)',
        actual_quotient_diagram_kernel='(u); square-zero free rank2 over C[[s20]]',
        explicit_isotropic_kernel_eta=str(eta),explicit_other_isotropic_coordinate=str(v),
        isotropic_change_determinant=str(s.factor(change)),
        eta_squared_zero_as_QQb_identity_not_just_point_test=True,
        full_actual_D_ring_Tor_basechange_nilmodule_or_second_embedded_prime_not_claimed=True,
        b_squared_equals3_exception_not_classified_or_declared_geometric=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in['all_passed','completed_worker_seconds','actual_completed_intersection','actual_quotient_diagram_kernel','seconds']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--recompute',action='store_true');a=p.parse_args();main(a.output.resolve(),a.recompute)
