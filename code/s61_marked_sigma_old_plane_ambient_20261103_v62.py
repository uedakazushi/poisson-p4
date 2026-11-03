#!/usr/bin/env python3
"""Full135 ambient Jacobian and quadratic incidence cuts on old-plane types."""
import argparse
from collections import Counter
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
P=32003


def dm(values,K):return DomainMatrix.from_rep(DDM(values,(len(values),len(values[0])),K))


def verify(model_folder,output):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    L=load_model(model_folder);K=s.QQ.algebraic_field(s.sqrt(-2),alias='h');h=K.unit
    zero,one=K.zero,K.one
    expr=lambda value:str(L.expression(value,K,s.Symbol('h')))
    def integral(value):return all(int(c.denominator)%P for c in value.to_list())
    def unit(value):
        if not integral(value):return False
        coefficients=value.to_list()
        a=coefficients[-1] if coefficients else 0
        b=coefficients[-2] if len(coefficients)>1 else 0
        norm=s.Rational(a)**2+2*s.Rational(b)**2
        return bool(norm and int(norm.p)%P and int(norm.q)%P)
    cases=[]
    source=L.original;betaindices,qindices=L.one,L.two;fixed=25
    for u,v,name in [(one,zero,'distinct-u'),(zero,one,'distinct-v'),
                     (h,coerce(K,s.Rational(2,3))*h,'double')]:
        normal=[zero]*13;normal[9]=h*u;normal[10]=u;normal[12]=v
        point=L.lift(h,normal,domain=K);assert point['all135_zero'] and not point['t']
        q=point['q'];actual_beta=point['beta']
        highest=[[zero]*16 for _ in range(135)];quad=[[zero]*16 for _ in range(16)]
        for row,form in enumerate(source):
            for mon,n,d in form:
                if sum(i in qindices for i in mon)!=2:continue
                beta=[betaindices.index(i) for i in mon if i in betaindices]
                value=coerce(K,s.Rational(n,d))
                for i in mon:
                    if i in qindices:value*=q[qindices.index(i)]
                if len(beta)==1:highest[row][beta[0]]+=value
                elif len(beta)==2:quad[beta[0]][beta[1]]+=value
        HD=dm(highest,K);H=HD.nullspace().transpose()
        hrows=list(H.transpose().rref()[1]);H=H*H.extract(hrows,range(H.shape[1])).inv()
        assert (HD*H).is_zero_matrix
        QD=dm(quad,K);assert (H.transpose()*(QD+QD.transpose())*H).is_zero_matrix
        qexpr=s.Matrix([K.to_sympy(value) for value in q]);physical=coords()['E']*qexpr
        spatial=s.Matrix(X[1:]);Q=s.zeros(4)
        for k,(i,j) in enumerate(pairs):
            Q[i,j]=sum(physical[k*10+a]*spatial[z]*spatial[w] for a,(z,w) in enumerate(m2))
        Q-=Q.T
        orbit=[]
        for A in sp_basis():
            dQ=A*Q+Q*A.T-sum(((A*spatial)[i]*Q.diff(spatial[i]) for i in range(4)),s.zeros(4))
            coefficients=s.Matrix([s.Poly(dQ[i,j],*spatial).coeff_monomial(spatial[z]*spatial[w])
                                   for i,j in pairs for z,w in m2])
            column=coords()['Pi']*coefficients
            column=[K.from_sympy(value) for value in column]
            orbit.append([value+q[i]*column[fixed] for i,value in enumerate(column)])
        O=dm([list(row) for row in zip(*orbit)],K)
        qcuts=list(O.transpose().rref()[1]);ocols=list(O.extract(qcuts,range(10)).rref()[1])
        orbit_minor=O.extract(qcuts,ocols).det();assert unit(orbit_minor)
        frozen={qindices[fixed],*[qindices[i] for i in qcuts],*[betaindices[i] for i in hrows]}
        variables=[i for i in range(62) if i not in frozen];location=dict(zip(variables,range(len(variables))))
        N=len(variables)
        generic=(H*dm([[coerce(K,i+2)] for i in range(H.shape[1])],K)).to_list()
        generic=[row[0] for row in generic]
        for kind,beta in [('actual',actual_beta),('generic-H',generic)]:
            values=dict(zip(betaindices,beta))|dict(zip(qindices,q))|{61:zero}
            encoded=[];linear=[[zero]*N for _ in range(135)];full=[[zero]*61 for _ in range(135)]
            fullvars=[i for i in range(62) if i!=qindices[fixed]];fullloc=dict(zip(fullvars,range(61)))
            for row,equation in enumerate(source):
                terms=[]
                for mon,n,d in equation:
                    ids=tuple(mon)+(61,)*(2-sum(i in qindices for i in mon))
                    coefficient=coerce(K,s.Rational(n,d));terms.append((ids,coefficient))
                    for variable,exponent in Counter(ids).items():
                        if variable==qindices[fixed]:continue
                        rest=list(ids);rest.remove(variable);value=coefficient*exponent
                        for i in rest:value*=values[i]
                        full[row][fullloc[variable]]+=value
                        if variable in location:linear[row][location[variable]]+=value
                encoded.append(terms)
            LD=dm(linear,K);pivots=list(LD.rref()[1]);rank=len(pivots);n=N-rank
            fullrank=dm(full,K).rank()
            print(name,kind,'FULLRANK',fullrank,'N',N,'RANK',rank,'NORMAL',n,flush=True)
            record=dict(case=name,beta=kind,H_dimension=H.shape[1],orbit_dimension=len(qcuts),
                        full_ambient_Jacobian_rank=fullrank,cut_variables=N,cut_Jacobian_rank=rank,
                        normal_dimension=n,incidence_dimension=H.shape[1]+len(qcuts),
                        beta_cut_indices=hrows,q_cut_indices=qcuts,orbit_minor=expr(orbit_minor))
            if n<=3:
                selected=list(LD.extract(range(135),pivots).transpose().rref()[1])
                block=LD.extract(selected,pivots);jdet=block.det();assert unit(jdet)
                T=LD.nullspace().transpose()
                free=[i for i in range(N) if i not in pivots]
                if n:T=T*T.extract(free,range(n)).inv()
                directions=[[row[j] for row in T.to_list()] for j in range(n)]
                def raw2(direction):
                    result=[]
                    for equation in encoded:
                        coefficient=zero
                        for ids,scalar in equation:
                            polynomial=[scalar,zero,zero]
                            for i in ids:
                                delta=direction[location[i]] if i in location else zero
                                polynomial=[polynomial[0]*values[i],
                                            polynomial[1]*values[i]+polynomial[0]*delta,
                                            polynomial[2]*values[i]+polynomial[1]*delta]
                            coefficient+=polynomial[2]
                        result.append(coefficient)
                    return result
                coefficients=[]
                squares=[raw2(direction) for direction in directions]
                for i in range(n):
                    for j in range(i,n):
                        if i==j:coefficients.append(squares[i])
                        else:
                            total=raw2([a+b for a,b in zip(directions[i],directions[j])])
                            coefficients.append([z-a-b for z,a,b in zip(total,squares[i],squares[j])])
                if coefficients:
                    raw=dm([list(row) for row in zip(*coefficients)],K)
                    residual=raw-LD.extract(range(135),pivots)*block.inv()*raw.extract(selected,range(len(coefficients)))
                    columns=list(residual.rref()[1]);rows=list(residual.extract(range(135),columns).transpose().rref()[1])
                    rdet=residual.extract(rows,columns).det()
                    record.update(quadratic_rank=len(columns),quadratic_monomials=len(coefficients),
                                  quadratic_minor=expr(rdet),quadratic_minor_norm_prime_unit=unit(rdet))
                    if len(columns)==n*(n+1)//2 and unit(rdet):record['transverse_length_upper']=1+n
                else:record['transverse_length_upper']=1
                record['jacobian_minor']=expr(jdet)
            cases.append(record)
            (output/'partial.json').write_text(json.dumps(cases,indent=2)+'\n')
    report=dict(id='S61_MARKED_SIGMA_OLD_PLANE_AMBIENT_20261103_v62',all_passed=True,cases=cases,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in [Path(__file__),ROOT/'code/s61_marked_sigma_point_lift_20261103_v62.py']},
                seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model-folder',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=verify(args.model_folder,args.output)
    print('PASS',report['id'],report['seconds'])
