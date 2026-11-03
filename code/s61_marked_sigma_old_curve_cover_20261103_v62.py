#!/usr/bin/env python3
"""Univariate full-incidence quadratic rank cover on u=1,v=k, h^2=-2.

Uses BOTH residue roots of h^2+2 modulo32003. A row-module determinant
certifies an entire algebraic k-line, not just rational parameter samples.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.ddm import DDM
from s61_marked_sigma_point_lift_20261103_v62 import load_model,coerce
from s61_binary_zero_pair_stratum_removal_20261102 import compile_span

ROOT=Path(__file__).resolve().parents[1]
P=32003


def dm(rows,K):return DomainMatrix.from_rep(DDM(rows,(len(rows),len(rows[0])),K))


def matrix_input(path,rows,F):
    polynomials=[]
    for row in rows:
        out=[]
        for value in row:
            assert value.denom.degree()==0
            out.append(value.numer/value.denom[(0,)])
        polynomials.append(out)
    with path.open('w') as stream:
        print(P,len(rows),len(rows[0]),0,file=stream);print(0,1,file=stream)
        for row in polynomials:
            nonzero=[(j,value) for j,value in enumerate(row) if value]
            print(len(nonzero),file=stream)
            for j,value in nonzero:
                degree=value.degree();print(j,degree+1,*[int(value.get((i,),F.dom.zero))%P for i in range(degree+1)],file=stream)


def verify(folder,output):
    start=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    L=load_model(folder);reference=ROOT/'reports/S61_MARKED_SIGMA_OLD_PLANE_AMBIENT_20261103_v62/report.json'
    old=next(c for c in json.loads(reference.read_text())['cases'] if c['case']=='distinct-u' and c['beta']=='actual')
    native=ROOT/'code/s61_marked_sigma_polynomial_rowmodule_20261103.cpp'
    executable=output/'row-module';compile_span(native,executable)
    k=s.Symbol('k');results=[]
    for root in [3083,28920]:
        base=s.GF(P);F=base.frac_field(k);kk=F.gens[0];h=coerce(F,root)
        values=[F.zero]*13;values[9]=h;values[10]=F.one;values[12]=kk
        point=L.lift(h,values,domain=F);assert point['all135_zero'] and not point['t']
        q,beta=point['q'],point['beta'];bi,qi=L.one,L.two
        frozen={qi[25],*[qi[i] for i in old['q_cut_indices']],*[bi[i] for i in old['beta_cut_indices']]}
        variables=[i for i in range(62) if i not in frozen];location=dict(zip(variables,range(44)))
        position=dict(zip(bi,beta))|dict(zip(qi,q))|{61:F.zero}
        linear=[[F.zero]*44 for _ in range(135)];encoded=[]
        scalar_cache={}
        for row,equation in enumerate(L.original):
            terms=[]
            for mon,n,d in equation:
                ids=tuple(mon)+(61,)*(2-sum(i in qi for i in mon))
                rational=s.Rational(n,d)
                if rational not in scalar_cache:scalar_cache[rational]=coerce(F,rational)
                coefficient=scalar_cache[rational]
                if any(not position[i] and i not in location for i in ids):continue
                zeros=sum(not position[i] for i in ids)
                if zeros<=2:terms.append((ids,coefficient))
                if zeros>=2:continue
                for variable,exponent in Counter(ids).items():
                    if variable not in location:continue
                    rest=list(ids);rest.remove(variable);value=coefficient*exponent
                    for index in rest:value*=position[index]
                    linear[row][location[variable]]+=value
            encoded.append(terms)
        LD=dm(linear,F);_,pivots=LD.rref();assert len(pivots)==43
        kernel=LD.nullspace().transpose();assert kernel.shape==(44,1)
        direction=[row[0] for row in kernel.to_list()]
        ring=F.field.ring;denominator=ring.one
        for value in direction:denominator=denominator.lcm(value.denom)
        polynomials=[value.numer*denominator.exquo(value.denom) for value in direction]
        gcd=ring.zero
        for value in polynomials:gcd=gcd.gcd(value)
        polynomials=[value.exquo(gcd) for value in polynomials]
        direction=[F.new(value) for value in polynomials]
        assert (LD*dm([[value] for value in direction],F)).is_zero_matrix
        raw2=[]
        for row in encoded:
            total=F.zero
            for ids,coefficient in row:
                values=[coefficient,F.zero,F.zero]
                for index in ids:
                    delta=direction[location[index]] if index in location else F.zero
                    values=[values[0]*position[index],values[1]*position[index]+values[0]*delta,
                            values[2]*position[index]+values[1]*delta]
                total+=values[2]
            raw2.append(total)
        discriminant=s.Poly(9*k*k+2*root,k,modulus=P)
        determinant_gcd=s.Poly(0,k,modulus=P);records=[]
        for j,value in enumerate(polynomials):
            if not value:continue
            columns=[i for i in range(44) if i!=j]
            rows=[[linear[i][col] for col in columns]+[raw2[i]] for i in range(135)]
            path=output/f'h{root}-omit{j}.in';matrix_input(path,rows,F)
            answer=output/f'h{root}-omit{j}.out';log=output/f'h{root}-omit{j}.log'
            with log.open('w') as stream:
                subprocess.run([str(executable.resolve()),str(path.resolve()),str(answer.resolve())],
                               check=True,timeout=60,stdout=stream,stderr=subprocess.STDOUT)
            line=next(line for line in answer.read_text().splitlines() if line.startswith('DETERMINANT_COEFFICIENTS '))
            data=list(map(int,line.split()[1:]));assert data[0]==len(data)-1
            determinant=s.Poly.from_list(list(reversed(data[1:])),k,modulus=P)
            determinant_gcd=s.gcd(determinant_gcd,determinant)
            remaining=determinant_gcd
            while remaining.degree()>0:
                common=s.gcd(remaining,discriminant)
                if common.degree()==0:break
                remaining=remaining.exquo(common)
            records.append(dict(omitted_column=j,determinant=str(determinant.as_expr()),
                                factors=str(determinant.factor_list()),gcd_degree=determinant_gcd.degree()))
            print('ROOT',root,'OMIT',j,'GCD',determinant_gcd.monic().factor_list(),'SECONDS',time.monotonic()-start,flush=True)
            if remaining.degree()==0:break
        results.append(dict(h_mod_P=root,discriminant=str(discriminant.as_expr()),
                            normal_direction_polynomials=list(map(str,polynomials)),
                            determinant_gcd=str(determinant_gcd.monic().as_expr()),
                            determinant_gcd_factors=str(determinant_gcd.monic().factor_list()),
                            remaining_extra_factors=str(remaining.monic().factor_list()),
                            all_distinct_root_finite_k_upper2_proved=remaining.degree()==0,records=records))
        (output/'partial.json').write_text(json.dumps(results,indent=2)+'\n')
    executable.unlink()
    report=dict(id='S61_MARKED_SIGMA_OLD_CURVE_COVER_20261103_v62',all_passed=True,
                prime=P,results=results,scope='Full135 actual fixed-beta univariate chart; singular discriminant points and infinity require separate certificates.',
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in [Path(__file__),ROOT/'code/s61_marked_sigma_point_lift_20261103_v62.py',reference,native]},
                seconds=time.monotonic()-start)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n');return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--model-folder',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    r=verify(args.model_folder,args.output);print('PASS',r['id'],r['seconds'])
