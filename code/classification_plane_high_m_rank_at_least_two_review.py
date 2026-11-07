#!/usr/bin/env python3
"""Independent all-m k1/rankM>=2 G review; Fraction only, no producer/CAS.

Reconstructs every high-m30/9 operator and wall, all three full350 nonlinear
sources and Pfaffians, actual Laurent normalization and literal contradictions.
"""
import ast
from fractions import Fraction as F
from itertools import combinations
import hashlib,json,resource,signal,time
from pathlib import Path
import classification_plane_rank_zero_rank_two_review as r

ROOT=Path(__file__).resolve().parents[1]
JOB='research/classification/jobs/plane_high_m_rank_at_least_two_review.json'
OUTPUT='reports/CLASSIFICATION_PLANE_HIGH_M_RANK_AT_LEAST_TWO_REVIEW_20261107.json'
PREFIX='reports/CLASSIFICATION_PLANE_RANK_ZERO_'


def load(name):return json.loads((ROOT/(PREFIX+name+'_20261107.json')).read_text())
def constant(v):
    v=r.RF(v);assert len(v.p)<=1 and v.q==(F(1),)
    return v.p[0] if v.p else F()
def model(name):
    a,b,c,z,w=[r.atom(1,i)for i in range(5)];t=r.RHO
    En=[{},{},{},z,w];T=[{}for _ in range(5)]
    vn=[{},{},{},r.q.mul(c,w),r.q.mul(b,z)]
    vb=[{},{},{},r.q.add(r.q.mul(b,z),r.q.mul(c,w)),r.q.scale(r.q.mul(b,w),-1)]
    index=0;V=vn;K={}
    if name=='NDN0':index=1;K=r.q.scale(r.q.mul(z,z),-F(1,2))
    elif name=='BSS0_tau':
        index=1;V=vb;T=[r.q.scale(a,-2*t/3),r.q.scale(b,4*t/3),r.q.scale(c,-2*t/3),{},{}];K=r.q.scale(r.q.mul(z,w),1-2*t)
    elif name=='BN0_tau':index=2;V=vb;T[2]=r.q.scale(b,2*t);K=r.q.scale(r.q.mul(w,w),F(1,2)-t)
    elif name=='R0_ND_nonzero_row_rho':T[0]=r.q.add(b,r.q.scale(c,t))
    elif name=='R0_ND_zero_row':pass
    elif name=='R0_B_rc_nonzero':V=vb;T[0]=c
    elif name=='R0_B_rcline_rho':V=vb;T[0]=r.q.scale(b,t)
    else:raise ValueError(name)
    U=[r.atom(int(i==index))for i in range(5)]
    return r.padd(r.wedge(T,En),r.wedge(U,V),r.q.antisymmetric({(3,4):K}))


def kernel_tensors(basis,certificate):
    K=r.matrix(certificate['complete_kernel'])
    return [r.q.antisymmetric({ij:r.q.add(*(r.q.scale(B.get(ij,{}),c)for B,c in zip(basis,column)))for ij in r.PAIRS})for column in zip(*K)]
def common_U(P):return all(not P.get(ij)for ij in r.PAIRS if 0 not in ij)and all(not r.q.der(P.get(ij,{}),0)for ij in r.PAIRS)


def operators(packet,b1,b2):
    rows=[]
    for saved in packet['all_eight_whole_high_m_rank2_strata']:
        name=saved['type']
        if name=='NDSS0_reused':continue
        P=model(name);assert r.q.vec(P)==[r.scalar(v)for v in saved['full150_P0']]
        assert not r.divisor(P)
        A,B=r.operator(P,b1),r.operator(P,b2)
        assert A==r.matrix(saved['full350_by30_D1'])and B==r.matrix(saved['full350_by9_D2'])
        a=saved['analysis'];cert={k:r.certificate(M,a[k])for k,M in [('D1',A),('D2',B)]}
        walls=[]
        if saved['variable']:
            product=r.scalar(a['D1']['minor_determinant'])*r.scalar(a['D2']['minor_determinant'])*r.RHO
            assert product==r.scalar(a['exception_polynomial'])
            factors=r.RF(1);roots=set()
            for f,multiplicity in a['all_irreducible_factors']:
                f=r.scalar(f);assert len(f.p)==2 and f.q==(F(1),);roots.add(-f.p[0]/f.p[1])
                for _ in range(multiplicity):factors=factors*f
            ratio=product/factors;assert ratio and len(ratio.p)==len(ratio.q)==1
            assert roots=={F(w['parameter'])for w in a['all_exceptional_roots']}
            for w in a['all_exceptional_roots']:
                t=F(w['parameter']);checked={k:r.certificate(r.specialization(M,t),w[k])for k,M in [('D1',A),('D2',B)]}
                walls.append(dict(parameter=str(t),certificates=checked))
        if name.startswith('R0'):
            full=[(None,a)]+[(F(w['parameter']),w)for w in a.get('all_exceptional_roots',[])]
            for value,C in full:
                if name=='R0_B_rcline_rho'and value in {F(1),F(5,3),F(3)}:continue
                first=kernel_tensors(b1,C['D1']);second=kernel_tensors(b2,C['D2'])
                if name=='R0_B_rcline_rho'and value==F(-1):
                    assert all(not Q.get(ij)for Q in first for ij in r.PAIRS if 0 not in ij and ij!=(1,2))
                else:
                    assert all(common_U(Q)for Q in first+second)
                    assert common_U(P)
        else:
            Bdim=saved['affine_first_symbol_dimension']
            assert Bdim+cert['D1']['nullity']+cert['D2']['nullity']+5<=23
            assert all(Bdim-1+w['certificates']['D1']['nullity']+w['certificates']['D2']['nullity']+5<=23 for w in walls)
        rows.append(dict(type=name,generic=cert,all_walls=walls,
                         affine_symbol_dimension=saved['affine_first_symbol_dimension']))
    return rows


def atom(n,k=None,c=1):
    e=[0]*n
    if k is not None:e[k]=1
    return {tuple(e):F(c)}if c else{}
def variable(n,k):return atom(n,k)
def compose(poly,images):
    n=len(images);out={}
    for ex,c in poly.items():
        term=atom(n,c=c)
        for i,power in enumerate(ex):
            for _ in range(power):term=r.q.mul(term,images[i])
        out=r.q.add(out,term)
    return out
def decode(poly):return {tuple(e):F(c)for e,c in poly if F(c)}


def full_tensor(P0,bases):
    n=len(bases);out={ij:{}for ij in r.PAIRS}
    for j,vector in enumerate([P0]+bases):
        for (ij,ex),value in zip(r.q.COLS,vector):
            c=F(value)
            if c:
                e=list(ex)+( [0]*n)
                if j:e[5+j-1]=1
                out[ij]=r.q.add(out[ij],{tuple(e):c})
    return r.q.antisymmetric(out)
def jacobi_and_divisor(P,n):
    J={}
    for i,j,k in combinations(range(5),3):
        J[i,j,k]=r.q.add(*(r.q.mul(P.get((a,l),{}),r.q.der(P.get((b,c),{}),l))
                          for a,b,c in [(i,j,k),(j,k,i),(k,i,j)]for l in range(5)))
    E={quad:r.q.add(*(r.q.scale(r.q.mul(atom(5+n,quad[k]),J[quad[:k]+quad[k+1:]]),(-1)**k)for k in range(4)))for quad in r.q.QUADS}
    allrows=[{ex[5:]:c for ex,c in E[quad].items()if ex[:5]==mon}for quad,mon in r.q.ROWS]
    D=r.q.add(*(r.q.scale(r.q.mul(atom(5+n,i),r.pf4(P,[j for j in range(5)if j!=i])),(-1)**i)for i in range(5)))
    return allrows,D


def nonlinear(flag,old,new,b1):
    records=[]
    for case in flag['cases']:
        rho=F(case['rho']);n=case['parameters'];k1=case['whole_P1_kernel_dimension'];j_r=k1+6
        # Reconstruct the wall tensor without evaluating a generic sample.
        vector0=[constant(r.evaluate(c,rho))for c in r.q.vec(model('R0_B_rcline_rho'))]
        assert vector0==list(map(F,case['full150_P0']))
        kernels=case['full_P1_kernel_tensor_basis'];P2=old['all9_source_basis']
        rowline=next(c for c in old['all_eight_whole_high_m_rank2_strata']if c['type']=='R0_B_rcline_rho')
        wall=next(c for c in rowline['analysis']['all_exceptional_roots']if F(c['parameter'])==rho)
        tensors=kernel_tensors(b1,wall['D1'])
        assert [list(map(F,v))for v in kernels]==[[constant(v)for v in r.q.vec(P)]for P in tensors]
        P=full_tensor(vector0,kernels+P2);J,D=jacobi_and_divisor(P,n)
        monomials=[]
        for m in case['parameter_monomials']:
            e=[0]*n
            if 'linear'in m:e[m['linear']]=1
            else:
                for j in m['quadratic_pair']:e[j]+=1
            monomials.append(tuple(e))
        saved=[{ex:F(c)for ex,c in zip(monomials,row)if F(c)}for row in case['full350_Jacobi_coefficient_matrix']]
        assert J==saved
        targets=[]
        for spatial,row in zip(case['all_Pf_mod_w_spatial_monomials'],case['full_Pf_mod_w_parameter_coefficient_matrix']):
            actual={ex[5:]:c for ex,c in D.items()if ex[:5]==tuple(spatial)}
            assert actual=={ex:F(c)for ex,c in zip(monomials,row)if F(c)}
            if actual:targets.append(tuple(spatial))
        assert set(targets)=={tuple(v['spatial_monomial'])for v in case['all_nonzero_actual_Pf_mod_w_membership_tests']}
        d,rr,l,m=[variable(n,j)for j in (8,j_r,11,13)]
        assert {ex[5:]:c for ex,c in P[(0,4)].items()if ex[:5]==(0,0,0,2,0)}==d
        assert {ex[5:]:c for ex,c in P[(1,2)].items()if ex[:5]==(0,0,0,2,0)}==rr
        for test in case['all_nonzero_actual_Pf_mod_w_membership_tests']:
            spatial=tuple(test['spatial_monomial']);target={ex[5:]:c for ex,c in D.items()if ex[:5]==spatial}
            combination=r.q.add(*(r.q.scale(J[int(i)],F(c))for i,c in test['original350_row_functional'].items()))
            residual={monomials[int(i)]:F(c)for i,c in test['nonmembership_remainder'].items()}
            assert target==r.q.add(combination,residual)
            assert test['member']==(not residual)
            if residual:
                assert spatial==(0,0,0,5,0)and residual==r.q.scale(r.q.mul(d,rr),-1)
        p=[variable(n,j)for j in range(n)]
        if rho==1:
            expected=[r.q.add(r.q.scale(rr,-6),r.q.mul(r.q.add(r.q.scale(p[0],2),r.q.scale(p[3],-1)),l),r.q.scale(r.q.mul(d,m),4)),
                      r.q.scale(r.q.mul(d,l),-3),
                      r.q.add(r.q.scale(rr,-4),r.q.mul(r.q.add(r.q.scale(p[0],2),r.q.scale(p[10],-2)),l),r.q.scale(r.q.mul(d,m),-4))]
        elif rho==F(5,3):
            expected=[r.q.add(r.q.scale(rr,-8),r.q.mul(r.q.add(r.q.scale(p[0],F(8,3)),r.q.scale(p[3],-F(4,3))),l),r.q.scale(r.q.mul(d,m),F(8,3))),
                      r.q.scale(r.q.mul(d,l),-F(11,3)),
                      r.q.add(r.q.scale(rr,-6),r.q.mul(r.q.add(r.q.scale(p[0],F(8,3)),r.q.scale(p[3],-F(1,3)),r.q.scale(p[10],-F(8,3))),l),r.q.scale(r.q.mul(d,m),-F(14,3)))]
        else:
            assert rho==3;expected=[r.q.add(r.q.scale(rr,-12),r.q.scale(r.q.mul(d,m),-4)),{},r.q.add(r.q.scale(rr,-10),r.q.scale(r.q.mul(d,m),6))]
        assert [J[i]for i in (51,121,122)]==expected
        # Entire coefficient Laurent normalization, not a point substitution.
        weights=[(-3,1),(-1,0),(-1,0),(0,0),(0,0)];allweights=[]
        for vector in [vector0]+kernels+P2:
            wset={tuple((2,-1)[k]+sum(weights[i][k]for i in ij)-sum(ex[i]*weights[i][k]for i in range(5))for k in range(2))
                  for value,(ij,ex)in zip(vector,r.q.COLS)if F(value)}
            assert len(wset)==1;allweights.append(list(next(iter(wset))))
        record=next(v for v in new['cases']if F(v['rho'])==rho)
        assert allweights[0]==[0,0]and allweights[1:]==record['complete_parameter_normalization_weights']
        assert allweights[9]==[-1,0]and allweights[1+j_r]==[0,-1]
        images=[variable(n,j)for j in range(n)];images[8]=images[j_r]=atom(n)
        equations=[compose(j,images)for j in J]
        for step in record['all_affine_elimination_steps']:
            if step.get('nonzero_constant_consequence'):
                assert r.q.add(*(r.q.scale(equations[int(i)],F(c))for i,c in step['current350_row_functional'].items()))==atom(n)
                continue
            for constraint in step['all_affine_consequences']:
                assert r.q.add(*(r.q.scale(equations[int(i)],F(c))for i,c in constraint['current350_row_functional'].items()))==decode(constraint['full_affine_equation'])
                j=constraint['solved_parameter'];assert decode(constraint['solution'])==r.q.scale(r.q.add(decode(constraint['full_affine_equation']),r.q.scale(variable(n,j),-1)),-1)
            substitutions=[variable(n,j)for j in range(n)]
            for constraint in step['all_affine_consequences']:substitutions[constraint['solved_parameter']]=decode(constraint['solution'])
            equations=[compose(j,substitutions)for j in equations];images=[compose(j,substitutions)for j in images]
        assert record['source_empty']
        assert images==[decode(j)for j in record['whole_source_original_parameter_images']]
        assert equations==[decode(j)for j in record['full350_residual_Jacobi_equations']]
        assert all(compose(j,images)==v for j,v in zip(J,equations))
        records.append(dict(rho=str(rho),full350_nonlinear_rows_reconstructed=True,
                            actual_Pf_mod_w_all_coefficients_and_functionals=True,
                            kappa_equals_minus_d_r_on_full_solution=True,
                            all150_Laurent_weights_and_inverse=True,
                            literal_rows_51_121_122_checked=True,normalized_original350_gives_one=True))
    return records


def main():
    resource.setrlimit(resource.RLIMIT_CPU,(175,176));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2);signal.alarm(180)
    start=time.monotonic();job=json.loads((ROOT/JOB).read_text())
    for row in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
    for manifest in ['plane_rank_zero_rank_two_high_m_manifest_20261107.json','plane_rank_zero_kappa_open_manifest_20261107.json']:
        for row in json.loads((ROOT/'research/classification/inputs'/manifest).read_text())['files']:
            assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256']
    rank3=load('FULLRANK_WEIGHTS');b1,b2,audit=r.source_bases(rank3)
    old,flag,new=load('RANK_TWO_HIGH_M'),load('RANK_TWO_HIGH_M_FLAG'),load('KAPPA_OPEN')
    assert not old['all_G24_closed']and not flag['all_three_whole_flag_exclusions_pass']
    assert old['all30_source_basis']==[[str(constant(v))for v in r.q.vec(B)]for B in b1]
    assert old['all9_source_basis']==[[str(constant(v))for v in r.q.vec(B)]for B in b2]
    rows=operators(old,b1,b2)
    # The fourth nonkernel source is exactly the old NDSS tau0 whole symbol.
    ndss=next(c for c in load('RANK_TWO_NONKERNEL_WEIGHTS')['cases']if c['type']=='NDSS')
    P=r.q.antisymmetric({ast.literal_eval(ij):{ex:r.evaluate(c,F(0))for ex,c in r.parse(v).items()if r.evaluate(c,F(0))}for ij,v in ndss['normal_form_P0'].items()})
    w=next(v for v in ndss['analysis']['all_exceptional_roots']if v['parameter']=='0')
    ndcert={k:r.certificate(r.operator(P,b),w[k])for k,b in [('D1',b1),('D2',b2)]}
    assert ndcert['D1']['nullity']==8 and ndcert['D2']['nullity']==0
    rank3cert=[]
    for c in rank3['cases']:
        x=[r.atom(1,i)for i in range(5)];a,b,cc,z,w=x
        Mn=[{},{},{},r.q.add(r.q.mul(b,z),r.q.mul(cc,w)),r.q.add(r.q.mul(a,z),r.q.scale(r.q.mul(b,w),-1))]
        ss=c['type']=='semisimple';U=[r.atom(int(i==(1 if ss else 0)))for i in range(5)]
        P=r.padd(r.wedge(U,Mn),r.q.antisymmetric({(3,4):r.q.mul(z,w)if ss else r.q.scale(r.q.mul(z,z),-F(1,2))}))
        wall=next(v for v in c['all_exceptional_roots']if v['tau']=='0')
        cert={k:r.certificate(r.operator(P,basis),wall[k])for k,basis in [('D1',b1),('D2',b2)]}
        assert cert['D1']['nullity']==8 and cert['D2']['nullity']==0
        rank3cert.append(dict(type=c['type'],certificates=cert,restored_symbol_dimension=11 if ss else 10,image_upper=24 if ss else 23))
    full=nonlinear(flag,old,new,b1)
    out=dict(task_id='G',all_passed=True,producer_or_CAS_imports=False,source_quotient=audit,
             exact_input_bindings=job['exact_inputs_and_hashes'],all_high_m_rank2_operators=rows,
             reused_NDSS0_full_reconstruction=ndcert,three_whole_kappa_sources=full,
             restored_rank3_tauzero=rank3cert,whole_high_m_rank2_image_upper=23,
             all_m_k1_rank_at_least_two_image_upper=24,
             complements=['k0','rankM<=1','rank2 plane restriction','fullE membership and fullG exhaustion'],
             old_partial_false_flags_preserved=True,elapsed_seconds=time.monotonic()-start,
             peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (ROOT/OUTPUT).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k]for k in ['all_passed','all_m_k1_rank_at_least_two_image_upper','elapsed_seconds','peak_rss_kib']}))


if __name__=='__main__':main()
