#!/usr/bin/env python3
"""Independent two whole weak-orbit full ambient Fraction review.

Only the previously independent sparse Fraction engine is reused. No producer,
producer ambient engine, current runner, SymPy or native CAS is imported.
"""
import ast
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import re
import resource
import time
import classification_degenerate_quadric_q4_nil31_I_review as q

ROOT=Path(__file__).resolve().parents[1]
JOB='research/classification/jobs/triple_weak_wall_ambient_separate_review.json'
SOURCE='reports/CLASSIFICATION_PLANE_TRIPLE_WEAK_WALL_FULL_AMBIENT_20261108.json'
OUTPUT='reports/CLASSIFICATION_PLANE_TRIPLE_WEAK_WALL_AMBIENT_SEPARATE_REVIEW_20261108.json'


def parse(text):
    names={'a':'h','b':'x1','c':'x2','z':'x3','w':'x4'}
    return q.parse(re.sub(r'\b(a|b|c|z|w)\b',lambda m:names[m[0]],str(text)))


def full_source(t,integer_scale):
    Q={(0,1):f'-({3*t+2})*a*z+w**2',(0,4):'-3*z**2',
       (1,3):f'({6*t+1})*z**2',(1,4):f'({9*t})*z*w'}
    # These are the actual affine-module particular solutions, not the unsolved W_base.
    if t==F(4,21):W=['a**2+7*b*w','-2*a*b','4*a*z/7+w**2/2','a*w/7-3*b*z']
    else:
        assert t==F(1,12)
        W=['a**2-20*b*w','-2*a*b','a*z/4-w**2','-a*w/2-3*b*z']
    p={ij:q.scale(parse(text),integer_scale)for ij,text in Q.items()}
    for j,text in zip([0,1,3,4],W):
        p[2,j]=q.scale(parse(text),integer_scale)
        p[j,2]=q.scale(p[2,j],-1)
    return q.antisymmetric({ij:p.get(ij,{})for ij in q.PAIRS})


def rowsparse(a):return [{j:v for j,v in enumerate(r)if v}for r in a]


def zero_product(a,b):
    br=rowsparse(b)
    for r in rowsparse(a):
        out={}
        for k,c in r.items():
            for j,d in br[k].items():out[j]=out.get(j,F(0))+c*d
        assert not any(out.values())


def decode(r):
    a=[[F(0)for _ in range(r['shape'][1])]for _ in range(r['shape'][0])]
    for key,value in r['entries'].items():
        i,j=map(int,key.split(','));a[i][j]=F(value)
    return a


def check_complete_rank(a,cert):
    rank=cert['rank'];rr,cc=cert['minor_rows'],cert['minor_columns']
    assert len(rr)==len(cc)==rank
    d=q.determinant([[a[i][j]for j in cc]for i in rr]);assert d and d==F(cert['minor_determinant'])
    free=[j for j in range(len(a[0]))if j not in cc]
    k=[[F(v)for v in row]for row in cert['complete_kernel']]
    assert len(k)==len(a[0]) and all(len(row)==len(free)for row in k)
    assert [[k[i][j]for j in range(len(free))]for i in free]==[[F(i==j)for j in range(len(free))]for i in range(len(free))]
    zero_product(a,k)
    return {'rank':rank,'nullity':len(free),'saved_nonzero_maximal_minor':str(d),'full_kernel_original_rows_zero':True,'free_coordinate_identity_and_complete_kernel':True},k


def pfaffian(p):
    result={}
    for i in range(5):
        a,b,c,d=[j for j in range(5)if j!=i]
        pf=q.add(q.mul(p[a,b],p[c,d]),q.scale(q.mul(p[a,c],p[b,d]),-1),q.mul(p[a,d],p[b,c]))
        result=q.add(result,q.scale(q.mul(q.atom(1,i),pf),(-1)**i))
    return result


def main():
    resource.setrlimit(resource.RLIMIT_CPU,(90,91));resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,)*2);begun=time.monotonic()
    job=json.loads((ROOT/JOB).read_bytes())
    for r in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/r['path']).read_bytes()).hexdigest()==r['sha256'],r['path']
    saved=json.loads((ROOT/SOURCE).read_bytes())
    assert [(tuple(r['pair']),tuple(r['exponent']))for r in saved['full_original_coefficient_order']]==q.COLS
    assert [(tuple(r['quad']),tuple(r['exponent']))for r in saved['full_original_equation_order']]==q.ROWS
    records=[]
    for row in saved['cases']:
        t=F(row['tau']);integer_scale=row['integer_scale'];p=full_source(t,integer_scale)
        for ij in q.PAIRS:assert p.get(ij,{})==parse(row['complete_actual_brackets'][str(ij)])
        assert all(not r for r in q.jacobi(p).values())
        Dcols=[]
        for ij,e in q.COLS:
            variation={k:dict(v)for k,v in p.items()}
            variation[ij]=q.add(variation[ij],{e+(0,0,0):F(1)})
            variation[ij[::-1]]=q.scale(variation[ij],-1)
            Dcols.append(q.jvec(q.jacobi(variation)))
        D=q.transpose(Dcols)
        vertical=[];orbit=[]
        for a,b in product(range(5),repeat=2):
            V={};G={}
            for i,j in q.PAIRS:
                V[i,j]=q.add(q.atom(1,i,b)if j==a else{},q.atom(-1,j,b)if i==a else{})
                G[i,j]=q.add(p.get((b,j),{})if i==a else{},p.get((i,b),{})if j==a else{},q.scale(q.mul(q.atom(1,b),q.der(p.get((i,j),{}),a)),-1))
            vertical.append(q.vec(V));orbit.append(q.vec(G))
        V=q.transpose(vertical);M=q.transpose(vertical+orbit+[q.vec(p)])
        matrices={'full350_D':D,'Euler25':V,'Euler_GL5_scalar51':M}
        certificates={};fullkernel=None
        for name,a in matrices.items():
            assert a==decode(row['complete_sparse_matrices'][name]),name
            review,k=check_complete_rank(a,row['certificates'][name]);certificates[name]=review
            if name=='full350_D':fullkernel=k
        assert fullkernel==decode(row['full_D_kernel_sparse'])
        assert row['full_D_kernel_free_columns']==[i for i in range(150)if i not in row['certificates']['full350_D']['minor_columns']]
        zero_product(D,M)
        actualD=pfaffian(p);assert actualD==parse(row['actual_signed_Euler_Pfaffian']) and actualD
        weights=[3,2,0,5,4]
        wc=[sum(weights[i]for i in ij)-sum(a*b for a,b in zip(e,weights))+3 for ij,e in q.COLS]
        wr=[sum(weights[i]for i in ij)-sum(a*b for a,b in zip(e,weights))+6 for ij,e in q.ROWS]
        assert wc==row['all150_tangent_weights'] and wr==row['all350_equation_weights']
        assert all(wc[i]==0 for i,c in enumerate(q.vec(p))if c)
        assert all(wr[i]==wc[j]for i,r in enumerate(D)for j,c in enumerate(r)if c)
        free=row['full_D_kernel_free_columns']
        assert all(wc[i]==wc[free[j]]for i,r in enumerate(fullkernel)for j,c in enumerate(r)if c)
        wg=[]
        for col in q.transpose(M):
            support={wc[i]for i,c in enumerate(col)if c};assert len(support)<=1
            wg.append(next(iter(support))if support else None)
        assert wg==row['all51_group_column_weights']
        kc=Counter(wc[j]for j in free);gc=Counter(wg[j]for j in row['certificates']['Euler_GL5_scalar51']['minor_columns'])
        assert dict(sorted((str(k),v)for k,v in kc.items()))==row['kernel_weight_dimensions']
        assert dict(sorted((str(k),v)for k,v in gc.items()))==row['group_weight_dimensions']
        assert kc==gc and row['normal_quotient_weight_dimensions']=={}
        ranks={name:v['rank']for name,v in certificates.items()};assert ranks==row['ranks_QQ']=={'full350_D':102,'Euler25':24,'Euler_GL5_scalar51':48}
        records.append({'tau':str(t),'integer_scale':integer_scale,'all10_actual_brackets_from_full_Q_and_actual_particular_W':True,
          'all350_original_Jacobi_zero':True,'all_full_D350x150_Euler150x25_M150x51_reconstructed':True,'all_D_M_zero':True,
          'rank_certificates':certificates,'all48_kernel_directions_original_substitution_and_free_identity':True,
          'all_actual_Pfaffian_coefficients_exact':True,'all150_point_torus_identity_and350_D_weights':True,
          'all_kernel_and_group_weights_equal_no_normal_quotient':True,'projective_actual_orbit_dimension':23,'projective_tangent_dimension':23,
          'smooth_unique23_component_and_actual_orbitclosure_proved':True})
        print('PASS',str(t),ranks,flush=True)
    source=Path(__file__).read_bytes()
    report={'task_id':'E_reduced_nonnormal','all_independent_arithmetic_passed':True,'cases':records,
      'whole_affine_W6_two_actual_orbits_bridge':'Separately handchecked in final review note, using complete module kernels/D images and global triangular shear inverses. The unsolved W_base display is not the actual particular.',
      'two_smooth23_actual_orbitclosure_components_proved':True,'mutual_or_prior_distinctness_proved':False,'component_count_increment_claimed':0,
      'inputs':job['exact_inputs_and_hashes'],'reviewer_sha256':hashlib.sha256(source).hexdigest(),'reviewer_source_bytes':source.decode(),
      'no_producer_current_helper_or_CAS_import':True,'seconds':time.monotonic()-begun,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    with(ROOT/OUTPUT).open('x')as f:json.dump(report,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps({'all_independent_arithmetic_passed':True,'seconds':report['seconds'],'maximum_RSS_KiB':report['maximum_RSS_KiB']}))


if __name__=='__main__':main()
