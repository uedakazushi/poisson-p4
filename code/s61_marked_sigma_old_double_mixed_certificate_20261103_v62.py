#!/usr/bin/env python3
"""Exact mixed-characteristic length-five certificate at the old double junction.

The 17 affine cuts replace delta_u4=0 by delta_u4+delta_u20=0.
All Taylor coefficients below are evaluated after that exact substitution.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.ddm import DDM
from s61_marked_sigma_point_lift_20261103_v62 import load_model,coerce

ROOT=Path(__file__).resolve().parents[1]


def dm(rows,K):return DomainMatrix.from_rep(DDM(rows,(len(rows),len(rows[0])),K))


def rational_pair(value):
    coefficients=list(value.to_list())
    assert len(coefficients)<=2
    return ([s.Rational(0)]*(2-len(coefficients))+[s.Rational(x) for x in coefficients])


def integral(value,prime=32003):
    return all(int(c.q)%prime for c in rational_pair(value))


def unit(value,prime=32003):
    b,a=rational_pair(value);norm=a*a+2*b*b
    return bool(norm and int(norm.p)%prime and int(norm.q)%prime)


def verify(folder,output):
    start=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    L=load_model(folder.resolve());K=s.QQ.algebraic_field(s.sqrt(-2),alias='h');h=K.unit;zero=K.zero
    v=[zero]*13;v[9]=h*h;v[10]=h;v[12]=coerce(K,s.Rational(2,3))*h
    point=L.lift(h,v,domain=K);assert point['all135_zero'] and not point['t']
    reference=ROOT/'reports/S61_MARKED_SIGMA_OLD_PLANE_AMBIENT_20261103_v62/report.json'
    old=next(c for c in json.loads(reference.read_text())['cases'] if c['case']=='double' and c['beta']=='actual')
    bi,qi=L.one,L.two;beta_cuts=[bi[i] for i in old['beta_cut_indices']]
    frozen={qi[25],*[qi[i] for i in old['q_cut_indices']],*beta_cuts}
    variables=[i for i in range(62) if i not in frozen];loc=dict(zip(variables,range(44)))
    values=dict(zip(bi,point['beta']))|dict(zip(qi,point['q']))|{61:zero}
    full=[[zero]*62 for _ in range(135)];encoded=[]
    for row,form in enumerate(L.original):
        terms=[]
        for mon,n,d in form:
            ids=tuple(mon)+(61,)*(2-sum(i in qi for i in mon));c=coerce(K,s.Rational(n,d))
            zeros=sum(not values[i] for i in ids)
            if zeros<=2:terms.append((ids,c))
            if zeros>=2:continue
            for index,count in Counter(ids).items():
                rest=list(ids);rest.remove(index);value=c*count
                for i in rest:value*=values[i]
                full[row][index]+=value
        encoded.append(terms)
    default=dm([[row[i] for i in variables] for row in full],K)
    piv=list(default.rref()[1]);assert len(piv)==41
    T=default.nullspace().transpose();free=[i for i in range(44) if i not in piv];T=T*T.extract(free,range(3)).inv()
    # The default quadratic GB has reduced x1-axis in precisely this free basis.
    axis=[row[1] for row in T.to_list()]
    qlambda=next(i for i in variables if i in qi and axis[loc[i]])
    m2=shutil.which('M2');assert m2;cases=[]
    monomials=[(2,0,0),(1,1,0),(1,0,1),(0,2,0),(0,1,1),(0,0,2)]
    for beta_cut in beta_cuts[:1]:
        assert beta_cut==3 and qlambda==19
        linear=[[row[i] for i in variables] for row in full]
        for i,row in enumerate(linear):row[loc[qlambda]]-=full[i][beta_cut]
        J=dm(linear,K);piv=list(J.rref()[1]);assert len(piv)==41
        selected=list(J.extract(range(135),piv).transpose().rref()[1]);block=J.extract(selected,piv)
        determinant=block.det();assert unit(determinant)
        T=J.nullspace().transpose();free=[i for i in range(44) if i not in piv];T=T*T.extract(free,range(3)).inv()
        directions=[]
        for j in range(3):
            direction=[zero]*62
            for i,index in enumerate(variables):direction[index]=T.to_list()[i][j]
            direction[beta_cut]=-direction[qlambda];directions.append(direction)
            assert direction[beta_cut]+direction[qlambda]==zero
            assert all(not sum(a*b for a,b in zip(row,direction)) for row in full)
        def square(direction):
            result=[]
            for terms in encoded:
                total=zero
                for ids,c in terms:
                    p=[c,zero,zero]
                    for i in ids:p=[p[0]*values[i],p[1]*values[i]+p[0]*direction[i],p[2]*values[i]+p[1]*direction[i]]
                    total+=p[2]
                result.append(total)
            return result
        squares=[square(direction) for direction in directions];coeff=[]
        for i in range(3):
            for j in range(i,3):
                if i==j:coeff.append(squares[i])
                else:
                    total=square([a+b for a,b in zip(directions[i],directions[j])])
                    coeff.append([z-a-b for z,a,b in zip(total,squares[i],squares[j])])
        raw=dm([list(row) for row in zip(*coeff)],K)
        residual=raw-J.extract(range(135),piv)*block.inv()*raw.extract(selected,range(6))
        rows=list(residual.transpose().rref()[1]);basis=residual.extract(rows,range(6)).to_list()
        assert len(rows)==5
        columns=list(dm(basis,K).rref()[1]);quadratic_minor=dm(basis,K).extract(range(5),columns).det()
        assert unit(quadratic_minor)
        all_scalars=[x for matrix in [full,J.to_list(),T.to_list(),raw.to_list(),residual.to_list(),directions]
                     for row in matrix for x in row]
        assert all(integral(x) for x in all_scalars)
        expression=lambda value:str(L.expression(value,K,s.Symbol('a'))).replace('**','^')
        polynomials=['+'.join('('+expression(c)+')'+''.join(f'*x{k}^{e}' for k,e in enumerate(mon) if e)
                              for mon,c in zip(monomials,row) if c) or '0' for row in basis]
        outputs={}
        for field,prefix in [('QQ','A=QQ[h];K=toField(A/ideal(h^2+2));a=promote(h,K);'),
                             ('root3083','K=ZZ/32003;a=promote(3083,K);'),
                             ('root28920','K=ZZ/32003;a=promote(28920,K);')]:
            script=output/f'{field}.m2';log=output/f'{field}.log'
            script.write_text(prefix+'R=K[x0,x1,x2];I=ideal('
                              +','.join(polynomials)+');assert(dim(R/I)==0);assert(numColumns basis(R/I)==5);'
                              +'print("DIM="|toString dim(R/I));'
                              +'print("HS="|toString hilbertSeries(R/I,Reduce=>true));'
                              +'print("LENGTH="|toString numColumns basis(R/I));'
                              +'print("GB="|toString gens gb I);print "MIXED CUT PASS";exit 0;\n')
            with log.open('w') as stream:subprocess.run([m2,'--script',str(script.resolve())],stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=60)
            outputs[field]=log.read_text();assert 'MIXED CUT PASS' in outputs[field]
        text=outputs['QQ']
        # These saved directions include the mixed beta coordinate.  They permit
        # an independent evaluation in the unmodified homogeneous135 equations.
        details=dict(global_variables=variables,linear_pivot_rows=selected,linear_pivot_columns=piv,
                     free_columns=free,quadratic_basis_rows=rows,quadratic_minor_columns=columns,
                     center=[expression(values[i]) for i in range(62)],
                     directions=[[expression(x) for x in row] for row in directions],
                     raw_quadratic_matrix=[[expression(x) for x in row] for row in raw.to_list()],
                     corrected_quadratic_matrix=[[expression(x) for x in row] for row in residual.to_list()])
        (output/'certificate.json').write_text(json.dumps(details,indent=2)+'\n')
        record=dict(beta_cut_global_index=beta_cut,q_link_global_index=qlambda,
                    cut=f'delta_u{beta_cut+1}+delta_u{qlambda+1}=0',quadratic_rank=len(rows),
                    beta_cuts_global=beta_cuts,q_cuts_local=old['q_cut_indices'],
                    chart_q25_global=qi[25],
                    ambient_Jacobian_rank=41,normal_dimension=3,cut_count=17,
                    mixed_substitution_verified_in_all135_tangent_directions=True,
                    jacobian_minor=expression(determinant),jacobian_minor_norm_prime_unit=unit(determinant),
                    quadratic_minor=expression(quadratic_minor),quadratic_minor_norm_prime_unit=unit(quadratic_minor),
                    scalar_denominators_prime_units=True,checked_scalar_count=len(all_scalars),
                    quadratic_basis=[[expression(value) for value in row] for row in basis],outputs=outputs)
        cases.append(record);(output/'partial.json').write_text(json.dumps(cases,indent=2)+'\n')
        print(record['cut'],text,'SECONDS',time.monotonic()-start,flush=True)
        if 'LENGTH=5' in text:break
    report=dict(id='S61_MARKED_SIGMA_OLD_DOUBLE_MIXED_CERTIFICATE_20261103_v62',all_passed=True,cases=cases,
                theorem='The actual mixed17-cut transverse fiber is Artinian of length at most5 at both residue roots.',
                limits=['Ambient t=0 additionally uses the two distinct17-dimensional incidence closures of generic lengths2 and3.'],
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in [Path(__file__),ROOT/'code/s61_marked_sigma_point_lift_20261103_v62.py',reference,
                                        L.input_path,ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']},
                seconds=time.monotonic()-start)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n');return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--model-folder',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    r=verify(args.model_folder,args.output);print('PASS',r['id'],r['seconds'])
