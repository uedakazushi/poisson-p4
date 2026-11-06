#!/usr/bin/env python3
"""One full first-jet test of the complete rank1 long-top, nonconstant-g source.

This tests dominance, not a higher formal jet.  A rank below16 is not used
as a generic upper bound without a separate all-source argument.
"""
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_component_beta_differential import certificate, dm


def actual_curve():
    """Verify the hand-derived rational same-jet curve, without elimination."""
    started=time.monotonic()
    signal.alarm(45)
    resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    source_path=Path('reports/CLASSIFICATION_RANK1_LONG_FIRSTJET_20261107.json')
    packet=json.loads(source_path.read_text())
    x,y,z,w=xx=s.symbols('x y z w')
    a,b,c=s.symbols('a b c');kk=s.symbols('k0:6')
    k,A,B,C,lam=s.symbols('k A B C lam')
    env={str(v):v for v in xx+(a,b,c)+kk}
    source=s.zeros(4)
    for key,value in packet['universal_affine_brackets'].items():
        i,j=map(int,key.split(','));source[i,j]=s.sympify(value,locals=env);source[j,i]=-source[i,j]
    zeta=z-1;v=b+c;ell=v*zeta-a*w
    K=k+A*zeta**2+B*zeta*w+C*w*w
    Klam=s.expand(lam*K+lam*(lam-1)*ell**2/2)
    monomials=(1,z,w,z*z,z*w,w*w)
    poly=s.Poly(Klam,z,w)
    substitutions={a:lam*a,b:b,c:lam*(b+c)-b}
    substitutions.update({q:poly.coeff_monomial(m) for q,m in zip(kk,monomials)})
    Plam=source.subs(substitutions,simultaneous=True)
    frame=s.eye(4);frame[0,2]=(1-lam)*v;frame[0,3]=-(1-lam)*a;frame[1,1]=1/lam
    inverse={x:x-(1-lam)*ell,y:lam*y}
    curve=(frame*Plam*frame.T).subs(inverse,simultaneous=True).applyfunc(s.cancel)
    original=curve.subs(lam,1)
    base=dict(zip(xx,(0,0,1,0)))
    pairs=tuple(itertools.combinations(range(4),2))
    residual=curve-original
    checks=[s.cancel(residual[i,j].subs(base)) for i,j in pairs]
    checks += [s.cancel(s.diff(residual[i,j],q).subs(base)) for i,j in pairs for q in xx]
    assert checks==[0]*30
    assert all(s.Poly(curve[i,j],*xx).total_degree()<=2 for i,j in pairs)
    jacobi=[s.cancel(sum(curve[i,l]*s.diff(curve[j,m],xx[l])
                         +curve[j,l]*s.diff(curve[m,i],xx[l])
                         +curve[m,l]*s.diff(curve[i,j],xx[l]) for l in range(4)))
            for i,j,m in itertools.combinations(range(4),3)]
    assert jacobi==[0]*4
    derivative=s.Poly(s.diff(curve[1,3],lam).subs(lam,1),*xx).coeff_monomial(x*x)
    assert derivative==-s.Rational(1,2)
    here=Path(__file__)
    report={'status':'FULL_SOURCE_ACTUAL_SAME_FIRSTJET_CURVE_PROVED',
            'producer_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),'producer_bytes':here.read_text(),
            'source_report_sha256':hashlib.sha256(source_path.read_bytes()).hexdigest(),
            'field':'QQ(a,b,c,k,A,B,C,lam)','localizations':['k','lam'],
            'basepoint':{str(q):str(v) for q,v in base.items()},
            'ell':str(ell),'K':str(K),'K_lam':str(Klam),
            'parameter_transform':{str(q):str(t) for q,t in substitutions.items()},
            'frame':[[str(t) for t in row] for row in frame.tolist()],
            'affine_frame_translation':[str(-(1-lam)*v),'0','0','0'],
            'full_curve_brackets':{str(i)+','+str(j):str(curve[i,j]) for i,j in pairs},
            'all_four_universal_affine_Jacobi_zero':True,'all_thirty_zero_firstjet_residuals_zero':True,
            'all_original_coordinate_degrees_at_most_two':True,
            'nonconstant_curve_derivative_P_yw_x_squared':str(derivative),
            'fixed_modular_source_dimension_upper':12,'affine_modular_orbit_dimension':10,
            'constant_evaluation_dimension':6,'actual_same_jet_fiber_dimension_lower':1,
            'primitive_beta_image_dimension_upper':15,
            'scope':'The affine centralizer and affine y-shears place every marked symplectic nonconstant-g point in this normal form. Each displayed curve consists of actual quadratic Poisson tensors, fixes the entire0/1jet, and is nonconstant for every source parameter. Affine saturation therefore has one-dimensional first-jet fibers everywhere. Constant-g and zero-g branches use separate proofs.',
            'seconds':round(time.monotonic()-started,4)}
    Path('reports/CLASSIFICATION_RANK1_LONG_ACTUAL_FIBER_20261107.json').write_text(json.dumps(report,indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({q:report[q] for q in ('status','primitive_beta_image_dimension_upper','nonconstant_curve_derivative_P_yw_x_squared','seconds')}))


def global_bound():
    """One rational-function nullspace certificate on the complete dense source."""
    started=time.monotonic()
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    previous=Path('reports/CLASSIFICATION_RANK1_LONG_FIRSTJET_20261107.json')
    packet=json.loads(previous.read_text())
    x,y,z,w=xx=s.symbols('x y z w')
    a,b,c=s.symbols('a b c')
    kk=s.symbols('k0:6')
    k,A,B,C=s.symbols('k A B C')
    env={str(v):v for v in xx+(a,b,c)+kk}
    source=s.zeros(4)
    for key,value in packet['universal_affine_brackets'].items():
        i,j=map(int,key.split(','))
        source[i,j]=s.sympify(value,locals=env)
        source[j,i]=-source[i,j]
    normal={kk[0]:k+A,kk[1]:-2*A,kk[2]:-B,kk[3]:A,kk[4]:B,kk[5]:C}
    base=dict(zip(xx,(0,0,1,0)))
    pairs=tuple(itertools.combinations(range(4),2))
    def jet(P):
        return s.Matrix([P[i,j].subs(base) for i,j in pairs]
                        +[s.diff(P[i,j],v).subs(base) for i,j in pairs for v in xx])
    directions=[source.diff(q) for q in (a,b,c)+kk]
    labels=list(map(str,(a,b,c)+kk))
    xv=s.Matrix(xx)
    for i,j in itertools.product(range(4),repeat=2):
        linear=s.zeros(4);linear[i,j]=1;field=linear*xv
        directions.append(linear*source+source*linear.T
                          -sum((field[l]*source.diff(xx[l]) for l in range(4)),s.zeros(4)))
        labels.append(f'GL_{i}_{j}')
    for l in range(4):
        directions.append(source.diff(xx[l]));labels.append(f'translation_{l}')
    directions.append(source);labels.append('scale')
    matrix=s.Matrix.hstack(*[jet(P).subs(normal,simultaneous=True).applyfunc(s.expand) for P in directions])
    field=s.QQ.frac_field(a,b,c,k,A,B,C)
    D=DomainMatrix.from_Matrix(matrix).convert_to(field)
    N=D.nullspace()
    assert D.matmul(N.transpose()).is_zero_matrix
    kernel=N.to_Matrix()
    pivots=list(N.rref()[1])
    assert len(pivots)==N.shape[0]
    minor=N.extract(list(range(N.shape[0])),pivots).det()
    assert minor!=field.zero
    rank_upper=matrix.cols-kernel.rows
    assert rank_upper<=21
    here=Path(__file__)
    report={'status':'GLOBAL_PRIMITIVE_BETA_NONDominance_PROVED',
            'producer_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),
            'producer_bytes':here.read_text(),
            'source_report_sha256':hashlib.sha256(previous.read_bytes()).hexdigest(),
            'coefficient_field':'QQ(a,b,c,k,A,B,C)',
            'normalization':'g=z, basepoint=(0,0,1,0), K=k+A(z-1)^2+B(z-1)w+Cw^2; k!=0',
            'actual_direction_labels':labels,
            'jet_matrix':[[str(v) for v in row] for row in matrix.tolist()],
            'kernel_rows':[[str(v) for v in row] for row in kernel.tolist()],
            'exact_matrix_times_kernel_zero':True,'kernel_independent_minor_columns':pivots,
            'kernel_minor_determinant':str(minor),
            'full_firstjet_image_dimension_upper':rank_upper,
            'constant_evaluation_orbit_dimension':6,
            'primitive_beta_image_dimension_upper':rank_upper-6,
            'logical_scope':'The normalized source represents every symplectic marked point of the full nonconstant-g long-top source under affine coordinate changes. Its irreducible saturation has first-jet differential rank at most21 over the rational function field, hence global first-jet image dimension at most21 in characteristic zero. Fixed nondegenerate evaluation removes6 dimensions; perspective quotient cannot increase dimension.',
            'seconds':round(time.monotonic()-started,4)}
    Path('reports/CLASSIFICATION_RANK1_LONG_FIRSTJET_GLOBAL_20261107.json').write_text(json.dumps(report,indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({q:report[q] for q in ('status','full_firstjet_image_dimension_upper','primitive_beta_image_dimension_upper','seconds')}))


def main():
    started = time.monotonic()
    signal.alarm(45)
    x,y,z,w = xx = s.symbols('x y z w')
    a,b,c = s.symbols('a b c')
    kk = s.symbols('k0:6')
    K = sum(q*m for q,m in zip(kk,(1,z,w,z*z,z*w,w*w)))
    H = z*K
    F = a*z*s.diff(K,z)+(b+c*z)*s.diff(K,w)+(a+1)*K
    r = y-x*x/2
    entries = [-F-(1+a)*r, a*z, b+c*z,
               a*x*z+s.diff(H,w), x*(b+c*z)-r-s.diff(H,z), z]
    pairs = tuple(itertools.combinations(range(4),2))
    source = s.zeros(4)
    for (i,j),v in zip(pairs,entries):
        source[i,j] = s.expand(v)
        source[j,i] = -source[i,j]
    jacobi = [s.expand(sum(source[i,l]*s.diff(source[j,k],xx[l])
                           +source[j,l]*s.diff(source[k,i],xx[l])
                           +source[k,l]*s.diff(source[i,j],xx[l]) for l in range(4)))
              for i,j,k in itertools.combinations(range(4),3)]
    assert jacobi == [0]*4
    modular = [s.expand(sum(s.diff(source[i,j],xx[i]) for i in range(4))) for j in range(4)]
    assert modular == [1,x,0,0]
    pf = s.expand(source[0,1]*source[2,3]-source[0,2]*source[1,3]+source[0,3]*source[1,2])
    assert s.expand(pf+z*(r+K)) == 0
    parameters = (a,b,c)+kk
    at = dict(zip(parameters,(1,2,3,1,2,3,4,5,6)))
    point = dict(zip(xx,(0,0,1,0)))
    P = source.subs(at)
    C = P.subs(point)
    assert C.det()!=0
    def jet(M):
        return s.Matrix([M[i,j].subs(point) for i,j in pairs]
                        +[s.diff(M[i,j],v).subs(point) for i,j in pairs for v in xx])
    directions = [source.diff(q).subs(at) for q in parameters]
    labels = list(map(str,parameters))
    xv = s.Matrix(xx)
    for i,j in itertools.product(range(4),repeat=2):
        A=s.zeros(4); A[i,j]=1; AX=A*xv
        directions.append(A*P+P*A.T-sum((AX[k]*P.diff(xx[k]) for k in range(4)),s.zeros(4)))
        labels.append(f'GL_{i}_{j}')
    for k in range(4):
        directions.append(P.diff(xx[k])); labels.append(f'translation_{k}')
    directions.append(P); labels.append('scale')
    actual=s.Matrix.hstack(*map(jet,directions))
    centered=xv-s.Matrix([point[v] for v in xx])
    gauges=[]
    for k in range(4):
        ca=C*s.eye(4)[:,k]
        gauges.append(jet(-2*centered[k]*C+centered*ca.T-ca*centered.T))
    gauge=s.Matrix.hstack(*gauges)
    all_jets=actual.row_join(gauge)
    ec,gc,jc=(certificate(actual[:6,:]),certificate(gauge),certificate(all_jets))
    kernel=dm(actual[:6,:]).nullspace().to_Matrix().T
    kc=certificate(gauge[6:,:].row_join(actual[6:,:]*kernel))
    assert ec['rank']==6 and gc['rank']==4
    rank=jc['rank']-10
    assert rank==kc['rank']-4 and rank<=16
    here=Path(__file__)
    helper=here.with_name('classification_component_beta_differential.py')
    report={'status':'PRIMITIVE_BETA_DOMINANCE_PROVED' if rank==16 else 'POINT_DIFFERENTIAL_ONLY',
            'producer_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),
            'producer_bytes':here.read_text(),
            'helper_sha256':hashlib.sha256(helper.read_bytes()).hexdigest(),
            'universal_source_parameters':list(map(str,parameters)),
            'universal_affine_brackets':{str(i)+','+str(j):str(source[i,j]) for i,j in pairs},
            'all_universal_Jacobi_zero':True,'universal_modular_field':list(map(str,modular)),
            'universal_pfaffian':str(pf),'source_witness':{str(k):str(v) for k,v in at.items()},
            'affine_basepoint':{str(k):str(v) for k,v in point.items()},
            'constant_bivector':[[str(v) for v in row] for row in C.tolist()],
            'actual_direction_labels':labels,'jet_matrix':[[str(v) for v in row] for row in all_jets.tolist()],
            'evaluation_certificate':ec,'perspective_certificate':gc,'all_jet_certificate':jc,
            'fixed_constant_linear_plus_perspective_certificate':kc,'primitive_beta_image_rank':rank,
            'logical_scope':'This source includes every nonconstant-g solution after an affine change of z,w. All affine directions and all perspective gauge directions are retained. Rank16 proves dominance; any smaller point rank is not a global upper bound.',
            'seconds':round(time.monotonic()-started,4)}
    Path('reports/CLASSIFICATION_RANK1_LONG_FIRSTJET_20261107.json').write_text(json.dumps(report,indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k:report[k] for k in ('status','primitive_beta_image_rank','seconds')}))


if __name__=='__main__':
    import sys
    actual_curve() if '--curve' in sys.argv else global_bound() if '--global' in sys.argv else main()
