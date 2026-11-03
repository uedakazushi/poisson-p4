#!/usr/bin/env python3
"""Independent Fraction physical source/polarization/kernel3 and W1 tangent review."""
import argparse
from fractions import Fraction as F
import itertools
import json
from math import factorial
from pathlib import Path
import resource
import signal
import time

import ambient_B18_zero_linear_contact4_review_20261103_v63 as own
from ambient_B18_actual_D_tangent_cone_support_review_20261103_v63 import inverse, matrix
from ambient_B18_Noether_fiber_even_quadratic_review_20261103_v63 import canonical_and_swap, exact_row_rank

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_W1_ACTUAL_TANGENT_INTERSECTION_20261103_v63_REVIEW'


def source_basis(cb):
    def raw(f=None,g=None,H=None):
        P=[F(0)]*150
        def put(i,j,mon,c):
            P[15*own.PI[i,j]+own.MI[tuple(sorted(mon))]]+=c
        if f is not None:
            put(1,2,f,F(1))
            if cb:
                for i,k in enumerate(f):
                    if k==1:
                        other=f[1-i]
                        put(2,3,(other,3),-F(1,2));put(2,4,(other,4),-F(1,2))
        if g is not None: put(3,4,g,F(1))
        if H is not None:
            for k,sign,pair in [(4,-1,(2,3)),(3,1,(2,4))]:
                count=H.count(k)
                if count:
                    mon=list(H);mon.remove(k)
                    put(*pair,mon,F(sign*count))
        return canonical_and_swap(P)
    cubics=[tuple([0,1,3,4][i] for i in m) for m in itertools.combinations_with_replacement(range(4),3) if any(i>=2 for i in m)]
    C=[raw(f=m) for m in [(0,0),(0,1),(1,1)]]+[raw(g=(0,0))]+[raw(H=m) for m in cubics]
    CG=[raw(g=m) for m in [(0,0),(0,1),(1,1)]] if cb else []
    return C,CG


def physical(M,point,delta,Q,order):
    n=order+1
    zero=[F(0)]*n
    def add(a,b):return [x+y for x,y in zip(a,b)]
    def mul(a,b):return [sum((a[j]*b[i-j] for j in range(i+1)),F(0)) for i in range(n)]
    def total(values):
        out=zero[:]
        for v in values:out=add(out,v)
        return out
    power=[[[F(i==j) for j in range(5)] for i in range(5)]]
    for _ in range(order):power.append(own.matrix_product(power[-1],M))
    G=[[[power[k][i][j]/factorial(k) for k in range(n)] for j in range(5)] for i in range(5)]
    Gi=[[[(-1)**k*power[k][i][j]/factorial(k) for k in range(n)] for j in range(5)] for i in range(5)]
    quadratic=[]
    for a,b in own.MONS:
        row=[]
        for i,j in own.MONS:
            v=mul(Gi[a][i],Gi[b][j])
            if i!=j:v=add(v,mul(Gi[a][j],Gi[b][i]))
            row.append(v)
        quadratic.append(row)
    P=[[a,b,c]+[F(0)]*(n-3) for a,b,c in zip(point,delta,Q)]
    changed=[[total(mul(P[15*k+l],quadratic[l][m]) for l in range(15)) for m in range(15)] for k in range(10)]
    wedge=[[add(mul(G[i][a],G[j][b]),[-v for v in mul(G[i][b],G[j][a])]) for a,b in own.PAIRS] for i,j in own.PAIRS]
    return [total(mul(wedge[i][k],changed[k][m]) for k in range(10)) for i in range(10) for m in range(15)]


def join(*matrices):
    return [sum((row for row in rows),[]) for rows in zip(*matrices)]


def main(output):
    start=time.monotonic()
    signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000))
    assert not output.exists()
    ap=ROOT/'reports/AMBIENT_W1_ACTUAL_TANGENT_INTERSECTION_20261103_v63/report.json'
    gates=[ROOT/'reports/AMBIENT_B18_GENERIC_W_FULL_COMPLETIONS_20261103_v63_REVIEW.json',
        ROOT/'reports/AMBIENT_B18_W_UNIFORM_DARBOUX_20261103_v63_REVIEW/report.json',
        ROOT/'reports/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63_REVIEW.json',
        ROOT/'reports/AMBIENT_B18_ALL_NONZERO_W_FULL_COMPLETIONS_20261103_v63/report.json']
    author,generic,uniform,globalgate,allw=[own.accepted(p) for p in [ap,*gates]]
    assert generic['CA_CB_cone_dim_edim_depth_multiplicity'][1]==38
    assert generic['N_cone_dim_edim_depth_multiplicity'][1]==36
    cp=ap.parent/'certificate.json.gz';c=own.read(cp)
    point=[F(0)]*150
    for pair,mon,val in [((0,2),(1,1),-1),((2,3),(0,1),1),((2,4),(0,0),1),((2,4),(3,3),1),((3,4),(1,1),1)]:
        point[15*own.PI[pair]+own.MI[mon]]=F(val)
    assert point==list(map(F,c['physical_point150']))
    value=[[F(0)]*4 for _ in range(4)];rows4=[0,2,3,4]
    for i,a in enumerate(rows4):
        for j,b in enumerate(rows4):
            if a<b:value[i][j]=point[15*own.PI[a,b]+own.MI[1,1]]
            elif b<a:value[i][j]=-point[15*own.PI[b,a]+own.MI[1,1]]
    _,rank4minor=inverse(value);assert rank4minor==1
    bases,dimensions,paths={}, {}, []
    for name,expected in [('CA',38),('CB',38),('N',36)]:
        rp=ap.parent/f'{name}-tangent.json.gz';paths.append(rp);r=own.read(rp)
        C,CG=source_basis(name=='CB')
        alpha=[F(0)]*20;alpha[0]=alpha[7]=-1;alpha[3]=alpha[11]=1;alpha[16]=F(1,3)
        assert alpha==list(map(F,r['source_point20']))
        indices=r['source_fixed20_indices']
        if name=='CB':
            p=[sum((alpha[i]*C[i][j] for i in [0,1,2,*range(4,20)]),F(0))-alpha[3]*sum((alpha[i]*CG[i][j] for i in range(3)),F(0)) for j in range(150)]
            fixed=[[C[i][j]-alpha[3]*CG[i][j] for j in range(150)] for i in range(3)]
            fixed+=[[-sum((alpha[i]*CG[i][j] for i in range(3)),F(0)) for j in range(150)]]+C[4:]
        else:
            p=[sum((a*row[j] for a,row in zip(alpha,C)),F(0)) for j in range(150)]
            fixed=[C[i] for i in indices]
        assert p==point
        assert indices==([i for i in range(20) if i!=2] if name=='N' else list(range(20)))
        critical=ROOT/f'reports/AMBIENT_B18_{name}_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';paths.append(critical)
        flags=own.read(critical)['flag_positions'];assert flags==r['source_flag_positions']
        generators=[]
        for i,j in flags:
            E=[[F(0)]*5 for _ in range(5)];E[i][j]=1;generators.append(E)
        nf=len(fixed)
        cols=fixed+[own.lie(E,point) for E in generators]
        J=list(map(list,zip(*cols)))
        assert J==matrix(r['source_Jacobian150_by_n']) and exact_row_rank(J)==nf+8
        k=list(map(F,r['source_kernel']))
        assert own.matvec(J,k)==[0]*150 and any(k)
        def direction(v):
            M=[[sum((v[nf+a]*E[i][j] for a,E in enumerate(generators)),F(0)) for j in range(5)] for i in range(5)]
            delta=[sum((v[a]*row[j] for a,row in enumerate(fixed)),F(0)) for j in range(150)]
            fv=[F(0)]*20
            for a,i in zip(v[:nf],indices):fv[i]=a
            Q=[-fv[3]*sum((fv[i]*CG[i][j] for i in range(3)),F(0)) for j in range(150)] if name=='CB' else [F(0)]*150
            return M,delta,Q,fv
        M,delta,Q,fv=direction(k)
        assert M==matrix(r['source_kernel_GL_matrix']) and fv==list(map(F,r['source_kernel_fixed_velocity20']))
        assert Q==list(map(F,r['full_CB_source_quadratic_velocity150']))
        jet=physical(M,point,delta,Q,3)
        assert [row[0] for row in jet]==point and all(row[1]==0 for row in jet)
        square,cube=[[row[d] for row in jet] for d in [2,3]]
        assert square==list(map(F,r['full150_kernel_square'])) and cube==list(map(F,r['full150_kernel_cube']))
        L2=own.lie(M,own.lie(M,point))
        assert cube==[a/6+b/2+cc for a,b,cc in zip(own.lie(M,L2),own.lie(M,own.lie(M,delta)),own.lie(M,Q))]
        mixed=[]
        for j in range(nf+9):
            ej=[F(i==j) for i in range(nf+9)]
            v=[a+b for a,b in zip(k,ej)]
            Mt,dt,Qt,_=direction(v);Me,de,Qe,_=direction(ej)
            total=physical(Mt,point,dt,Qt,2);pure=physical(Me,point,de,Qe,2)
            mixed.append([a[2]-b-cp[2] for a,b,cp in zip(total,square,pure)])
        H=list(map(list,zip(*mixed)))
        assert H==matrix(r['full150_by_n_mixed2'])
        second=join(J,[[a] for a in square],H)
        assert exact_row_rank(second)==r['source2_span_rank']==expected-1
        full=join(second,[[a] for a in cube])
        assert exact_row_rank(full)==expected
        B=[[row[j] for j in r['independent_column_indices']] for row in full]
        assert B==matrix(r['actual_tangent150_by_dimension'])
        _,minor=inverse([B[i] for i in r['unit_minor_rows']])
        assert minor==F(r['unit_minor_determinant'])
        bases[name]=B;dimensions[name]=dict(source_rank=nf+8,second_rank=expected-1,full_rank=expected,minor=str(minor))
    A,B,N=[bases[n] for n in ['CA','CB','N']]
    assert exact_row_rank(join(A,N))==38 and exact_row_rank(join(B,N))==38
    assert exact_row_rank(join(A,B))==39
    TD=matrix(c['full150_actual_D_tangent_basis'])
    assert exact_row_rank(TD)==37 and exact_row_rank(join(A,TD))==38 and exact_row_rank(join(B,TD))==38
    assert exact_row_rank(join(N,TD))==37
    inputs=[Path(__file__),Path(own.__file__),ROOT/'code/ambient_B18_Noether_fiber_even_quadratic_review_20261103_v63.py',
        ROOT/'code/ambient_B18_actual_D_tangent_cone_support_review_20261103_v63.py',ROOT/f'research/notes/{STEM}.md',ap,cp,*gates,*paths]
    report=dict(id=STEM,all_passed=True,
        independent_raw_fixed20_CA_CB_N_and_physical_W1_point_reconstructed=True,
        independent_actual_source_Jacobian_and_unique_kernel_all150=True,
        independent_inverse_coordinate_polarization_all_kernel_mixed2=True,
        independent_one_true_kernel_cube_per_family_and_Lie_comparison=True,
        CB_full_ratio_times_f_quadratic_source_term_retained=True,
        actual_Ac8_Ac7_upper_bound_and_proper_source_completeness_used=True,
        physical_tangent_records=dimensions,actual_CA_CB_N_tangent_dimensions=[38,38,36],
        actual_CA_plus_CB_span_dimension=39,actual_D_tangent_dimension=37,actual_Dred_tangent_dimension=36,
        rank4_value_minor='1',actual_Dred_germ_equals_N=True,nilradical_cotangent_image_dimension=1,
        every_nonzero_W_by_scheme_preserving_GL_scalar_transport=True,
        full_D_type_depth_nil_index_nil_support_or_embedded_components_claimed=False,
        verification_kind='Own raw coefficient source basis, Fraction physical polarization/kernel3 and exact ranks/minors; producer not rerun, no higher source jets.',
        input_sha256={str(p.relative_to(ROOT)):own.sha(p) for p in inputs},seconds=time.monotonic()-start)
    output.write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k] for k in ['all_passed','actual_CA_CB_N_tangent_dimensions','actual_D_tangent_dimension','nilradical_cotangent_image_dimension','seconds']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path)
    main(parser.parse_args().output)
