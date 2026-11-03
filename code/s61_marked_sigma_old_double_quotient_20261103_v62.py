#!/usr/bin/env python3
"""Universal t=0 square-zero thickening of the old distinct-root incidence."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.ddm import DDM
from s61_marked_sigma_point_lift_20261103_v62 import load_model,coerce
from toric_closed_jet_coordinates_20261102 import coords,pairs,m2
from s61_model import X
from weight1_structure import sp_basis

ROOT=Path(__file__).resolve().parents[1]


def dm(values,K):return DomainMatrix.from_rep(DDM(values,(len(values),len(values[0])),K))


def verify(model_folder,output):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    L=load_model(model_folder);K=s.QQ.algebraic_field(s.sqrt(-2),alias='h');h=K.unit
    zero,one=K.zero,K.one;normals=[zero]*13;normals[9]=h;normals[10]=one
    point=L.lift(h,normals,domain=K);assert point['all135_zero'] and not point['t']
    q=point['q'];bi,qi=L.one,L.two;fixed=25;chart=[i for i in range(45) if i!=fixed]
    highest=[[zero]*16 for _ in range(135)]
    pure=[[zero]*44 for _ in range(135)]
    highest_terms=[]
    for row,form in enumerate(L.original):
        terms=[]
        for mon,n,d in form:
            qs=[qi.index(i) for i in mon if i in qi]
            if len(qs)!=2:continue
            bs=[bi.index(i) for i in mon if i in bi]
            scalar=coerce(K,s.Rational(n,d));terms.append((qs,bs,scalar))
            value=scalar*q[qs[0]]*q[qs[1]]
            if len(bs)==1:highest[row][bs[0]]+=value
            if not bs:
                for k in range(2):
                    index=qs[k]
                    if index!=fixed:pure[row][chart.index(index)]+=scalar*q[qs[1-k]]
        highest_terms.append(terms)
    M=dm(highest,K);H=M.nullspace().transpose();hrows=list(H.transpose().rref()[1])
    H=H*H.extract(hrows,range(H.shape[1])).inv();assert H.shape==(16,9)
    C=dm(pure,K);ker=C.nullspace().transpose()
    physical=coords()['E']*s.Matrix([K.to_sympy(value) for value in q]);y=s.Matrix(X[1:]);Q=s.zeros(4)
    for k,(i,j) in enumerate(pairs):Q[i,j]=sum(physical[k*10+a]*y[z]*y[w] for a,(z,w) in enumerate(m2))
    Q-=Q.T;orbit=[]
    for A in sp_basis():
        dQ=A*Q+Q*A.T-sum(((A*y)[i]*Q.diff(y[i]) for i in range(4)),s.zeros(4))
        coefficients=s.Matrix([s.Poly(dQ[i,j],*y).coeff_monomial(y[z]*y[w]) for i,j in pairs for z,w in m2])
        col=[K.from_sympy(value) for value in coords()['Pi']*coefficients]
        orbit.append([col[i]+q[i]*col[fixed] for i in chart])
    O=dm([list(row) for row in zip(*orbit)],K);ocols=list(O.rref()[1]);assert len(ocols)==8
    basis=O.extract(range(44),ocols);candidates=[]
    for j in range(ker.shape[1]):
        column=ker.extract(range(44),[j])
        if basis.hstack(column).rank()>basis.shape[1]:
            values=column.to_list();scale=next(row[0] for row in values if row[0])
            column=column.mul(one/scale);basis=basis.hstack(column)
            eta=[zero]*45
            for i,index in enumerate(chart):eta[index]=column.to_list()[i][0]
            candidates.append(eta)
    print('Q45',C.rank(),ker.shape[1],'EXTRA',len(candidates),flush=True)
    comp=[i for i in range(16) if i not in hrows]
    selected=list(M.extract(range(135),comp).transpose().rref()[1])
    block=M.extract(selected,comp);assert block.shape==(7,7);inverse=block.inv().to_list()
    z=s.symbols('z0:9');R=K.poly_ring(*z);F=K.frac_field(*z)
    Hr=H.convert_to(R).to_list()
    beta=[sum((c*R.ring.gens[j] for j,c in enumerate(row) if c),R.zero) for row in Hr]
    derivatives=[];beta_directions=[]
    for number,eta in enumerate(candidates):
        raw=[]
        for terms in highest_terms:
            value=R.zero
            for qs,bs,c in terms:
                bprod=R.one
                for b in bs:bprod*=beta[b]
                value+=c*(eta[qs[0]]*q[qs[1]]+q[qs[0]]*eta[qs[1]])*bprod
            raw.append(value)
        db=[R.zero]*16
        for i,index in enumerate(comp):db[index]=-sum((c*raw[row] for c,row in zip(inverse[i],selected) if c),R.zero)
        total=[]
        for i,terms in enumerate(highest_terms):
            value=raw[i]
            for qs,bs,c in terms:
                for k,index in enumerate(bs):
                    term=c*q[qs[0]]*q[qs[1]]*db[index]
                    for j,b in enumerate(bs):
                        if j!=k:term*=beta[b]
                    value+=term
            total.append(value)
        assert not any(total[row] for row in selected)
        derivatives.append(total);beta_directions.append(db)
        print('CANDIDATE',number,'TERMS',sum(len(p) for p in total),flush=True)
    for degree in [0,1]:
        shifts=[(0,)*9] if degree==0 else [tuple(int(i==j) for i in range(9)) for j in range(9)]
        columns=len(candidates)*len(shifts);coefficients={}
        for j,rows in enumerate(derivatives):
            for k,shift in enumerate(shifts):
                for i,poly in enumerate(rows):
                    for mon,c in poly.items():
                        exponent=tuple(a+b for a,b in zip(mon,shift))
                        coefficients.setdefault((i,exponent),{})[j*len(shifts)+k]=c
        conditions=DomainMatrix(dict(enumerate(coefficients.values())),(len(coefficients),columns),K)
        kernel=conditions.nullspace().transpose()
        print('DEGREE',degree,'CONDITIONS',conditions.shape,'KERNEL',kernel.shape[1],flush=True)
        if kernel.shape[1]:break
    assert kernel.shape[1]>=1
    scalars=[row[0] for row in kernel.to_list()]
    scale=next(value for value in scalars if value);scalars=[value/scale for value in scalars]
    combination=[]
    for j in range(len(candidates)):
        combination.append(sum((R.ring.from_dict({shift:scalars[j*len(shifts)+k]})
                                for k,shift in enumerate(shifts)),R.zero))
    eta=[sum((c*v[i] for c,v in zip(combination,candidates)),R.zero) for i in range(45)]
    db=[sum((c*v[i] for c,v in zip(combination,beta_directions)),R.zero) for i in range(16)]
    total=[sum((c*v[i] for c,v in zip(combination,derivatives)),R.zero) for i in range(135)]
    assert not any(total)
    assert any(combination) and basis.shape[1]==O.rank()+len(candidates)
    # beta+epsilon*db, q+epsilon*eta, epsilon^2=0, t=0 is exact full135.
    expression=lambda value:str(L.expression(value,K,s.Symbol('h')))
    fexpression=lambda value:str(L.expression(F.convert(value,R),F,s.Symbol('h')))
    report=dict(id='S61_MARKED_SIGMA_OLD_DOUBLE_QUOTIENT_20261103_v62',all_passed=True,
                field='QQ(h), h^2=-2',normal_origin='old distinct roots u=1,v=0',
                Q45_projective_jacobian_rank=C.rank(),Q45_tangent_dimension=ker.shape[1],
                projective_orbit_dimension=8,extra_q_candidates=len(candidates),
                constraint_shape=list(conditions.shape),constraint_kernel_dimension=kernel.shape[1],
                H=[[expression(v) for v in row] for row in H.to_list()],
                eta_beta_parameter_degree=degree,
                eta_q=list(map(fexpression,eta)),delta_beta=list(map(fexpression,db)),
                normal_not_in_projective_orbit=True,all135_highest_first_derivatives_zero=True,
                exact_dual_number_thickening_with_t_zero=True,
                same_incidence_generic_length_lower=2,
                limitations=['The local closed quotient is transported over the smooth orbit/kernel incidence.',
                             'An upper2 certificate is separately needed for ambient equality.',
                             'This does not by itself settle the double-root junction.'],
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in [Path(__file__),ROOT/'code/s61_marked_sigma_point_lift_20261103_v62.py',
                                        ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']},
                seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model-folder',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();r=verify(args.model_folder,args.output)
    print('PASS',r['id'],r['seconds'])
