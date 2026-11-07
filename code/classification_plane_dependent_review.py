#!/usr/bin/env python3
"""Independent six dependent global-vector sources; no producer/CAS imports.
The three unresolved whole sources are retained until the Ore/module proof.
"""
from fractions import Fraction as F
from pathlib import Path
import ast,hashlib,json,resource,signal,time
from classification_unimodular_logarithmic_review import RF
import classification_plane_m2_completion_review as M
q,H,B=M.q,M.H,M.B
ROOT=Path(__file__).resolve().parents[1]
SOURCE='reports/CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_DEPENDENT_20261107.json'
JOB='research/classification/jobs/plane_dependent_six_source_review.json'
OUT='reports/CLASSIFICATION_PLANE_DEPENDENT_SIX_SOURCE_REVIEW_20261107.json'
U=H.previous.RHO

def model(name):
    p=lambda t:H.parse(t,{'tau':U},RF)
    V=[p('a'),p('b'),p('-2*c'),{},{}] if '_H_' in name else [p('b'),{},{},{},{}]
    if name=='SS_N_LINE':coef=U
    elif 'TRACE_ONE'in name:coef=RF(1)
    else:coef=RF(0)
    T=[q.scale(v,coef)for v in V];En=[{},{},{},p('z'),p('w')]
    normal=[{},{},{},p('z'),p('-w')]if name.startswith('SS')else[{},{},{},{},p('z')]
    K=p('-5*z*w')if name=='SS_H_RESONANCE'else p('5*z**2/2')if name=='NIL_H_RESONANCE'else p('z**2')if name=='NIL_N_TRACE_ONE_K1'else{}
    return H.Padd(H.wedge(T,En),H.wedge(V,normal),q.antisymmetric({(3,4):K}))

def verify(Ms,c,env,base,ss):
    certs={};product=RF(1)
    for mat,k in zip(Ms,('D1','D2','actual_shear')):
        certs[k],v=H.certificate(mat,c[k],env,RF);product=product*v
    kernel=H.matrix(c['D1']['complete_kernel'],env,RF)
    C=kernel[27:30];D=kernel[18:21]
    assert C==H.matrix(c['full_C_projection'],env,RF)and D==H.matrix(c['full_D_projection'],env,RF)
    cz=not any(any(r)for r in C);dz=not any(any(r)for r in D)
    assert cz==c['all_C_zero']and dz==c['all_D_zero']and c['required_C_nonzero']and c['required_D_nonzero']==ss
    empty=cz or ss and dz;assert empty==c['G_empty_required_normal_hyperplane']
    d=c['D1']['nullity']+c['D2']['nullity'];bound=base+5+d
    assert d==c['higher_fiber_upper_bound']and bound==c['PGL_image_upper_bound']
    return dict(complete_kernel_certificates=certs,C_projection_zero=cz,D_projection_zero=dz,
                required_hyperplane_factor_forces_empty=empty,higher_fiber_upper=d,coarse_image_upper=bound,
                closed_by_linear_or_hyperplane=empty or bound<=24),product

def main():
    started=time.monotonic();signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    assert not(ROOT/OUT).exists();job=json.loads((ROOT/JOB).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    packet=json.loads((ROOT/SOURCE).read_text());assert not packet['all_passed']
    assert hashlib.sha256(packet['producer_bytes'].encode()).hexdigest()==packet['producer_sha256']
    first,second=M.bases();records=[]
    for c in packet['cases']:
        name=c['type'];P=model(name);assert not any(q.jacobi(P).values())and not H.divisor(P)
        expected=q.antisymmetric({ast.literal_eval(k):H.parse(v,{'tau':U},RF)for k,v in c['normal_form_P0'].items()});assert q.vec(P)==q.vec(expected)
        Ms=[H.previous.operator(P,B.convert_basis(basis,RF))for basis in (first,second)]+[H.shear(P,first,RF)]
        for mat,k in zip(Ms,('full_D1_350_by30','full_D2_350_by9','full_actual_shear_30_by6')):assert mat==H.matrix(c[k],{'tau':U},RF)
        r,prod=verify(Ms,c['analysis'],{'tau':U},c['affine_first_symbol_dimension_upper'],name.startswith('SS'));r['type']=name
        if name=='SS_N_LINE':
            a=c['analysis'];factors=M.factor_product(prod,a,{'tau':U},RF);roots=[]
            for text,f in factors:
                assert f.q==(F(1),)and len(f.p)==2;value=-f.p[0]/f.p[1]
                rows=[w for w in a['all_exceptional_roots']if F(w['parameter'])==value];assert len(rows)==1
                current=[H.convert_matrix(mat,lambda v:M.at_rf(RF(v),value))for mat in Ms]
                w,_=verify(current,rows[0],{'tau':RF(value)},6,True);w['parameter']=str(value);roots.append(w)
            assert len(roots)==6;r['all_six_artificial_minor_or_true_cohomology_roots']=roots
            assert set(w['parameter']for w in roots)=={'1/3','-1/3','1','-1','3','-3'}
        records.append(r);print(name,'full350/30/9/kernel PASS',r['coarse_image_upper'],'empty',r['required_hyperplane_factor_forces_empty'],flush=True)
    assert len(records)==6 and sum(r['closed_by_linear_or_hyperplane']for r in records)==3
    source=Path(__file__).read_text();result=dict(task_id='G',all_passed=True,all_six_full_source_arithmetic=records,
        whole_dependent_G_completion=False,remaining_whole_sources=['SS_N_LINE','NIL_N_TRACE_ZERO','NIL_N_TRACE_ONE_K0'],
        scope='Full six sources/whole tau line and all minor roots exactPASS; only Hss<=20,Hnil<=21,trace1K1empty closed here. Ore/module support is separate.',
        all_old_false_and_coarse_bounds_preserved=True,no_CAS_or_producer_imports=True,inputs=job['exact_inputs_and_hashes'],
        reviewer_sha256=hashlib.sha256(source.encode()).hexdigest(),reviewer_bytes=source,
        seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (ROOT/OUT).write_text(json.dumps(result,indent=2)+'\n');print('PASS',OUT,round(result['seconds'],3),flush=True)
if __name__=='__main__':main()
