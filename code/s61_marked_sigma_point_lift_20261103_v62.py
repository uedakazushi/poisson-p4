#!/usr/bin/env python3
"""Portable inverse of the frozen all-h marked-square Smith model.

The chart has q25=-1. Thirteen normal variables v0,...,v12 and rho=u5
are retained. h and all values may lie in QQ, a number field, a rational
function field, or GF(p). No native binary or historical output is needed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import sympy as s
from s61_model import S61
from s61_marked_sigma_pencil_model_20261102 import prepare

ROOT=Path(__file__).resolve().parents[1]


def coerce(domain,value):
    # SymPy GF.convert(Rational(1,2)) can silently truncate to zero. Map
    # numerator and denominator separately and divide inside the target field.
    if isinstance(value,s.Rational):
        return domain.convert(int(value.p))/domain.convert(int(value.q))
    try:
        return domain.convert(value)
    except Exception:
        if hasattr(domain,'dom'):
            try:
                return domain.convert(value,domain.dom)
            except Exception:
                pass
        return domain.from_sympy(value)


class PointLifter:
    """Read prepare()'s checked J/T/U and evaluate its exact weighted inverse."""
    def __init__(self,input_path):
        self.input_path=Path(input_path)
        tokens=iter(self.input_path.read_text().split())
        integer=lambda:int(next(tokens))
        n,rank=integer(),integer()
        assert (n,rank)==(14,28)
        self.weights=[integer() for _ in range(45)]
        self.normal_weights=[integer() for _ in range(n)]
        self.selected=[integer() for _ in range(rank)]
        def core():
            return [(s.Rational(next(tokens)),integer()) for _ in range(integer())]
        def matrix(rows,columns):
            return [[core() for _ in range(columns)] for _ in range(rows)]
        self.J=matrix(135,45);self.T=matrix(45,13);self.U=matrix(45,28)
        self.equations=[]
        for _ in range(135):
            row=[]
            for _ in range(integer()):
                i,j,rho_power=integer(),integer(),integer()
                row.append((i,j,rho_power,core()))
            self.equations.append(row)
        self.row_weights=[integer() for _ in range(135)]
        assert next(tokens,None) is None
        assert self.normal_weights==[1]*7+[2]*4+[3]*2+[1]
        self.model=S61()
        self.one=[i for i,w in enumerate(self.model.reduced['parameter_weights']) if w==1]
        self.two=[i for i,w in enumerate(self.model.reduced['parameter_weights']) if w==2]
        self.original=self.model.reduced['equations']+[json.loads(
            (ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']]

    def lift(self,h,normals,rho=0,domain=s.QQ,verify=True):
        """Return domain-valued beta16/q45/t and all raw135 residuals.

        Use F.convert(K.unit,K) for h when F=K.frac_field(...), K quadratic.
        Input normals are in v0,...,v12 order, never including rho.
        """
        assert len(normals)==13
        h=coerce(domain,h);rho=coerce(domain,rho)
        normals=[coerce(domain,value) for value in normals]
        zero,one=domain.zero,domain.one
        powers={0:one}
        def evaluate_core(poly):
            total=zero
            for coefficient,power in poly:
                if power not in powers:powers[power]=h**power
                total+=coerce(domain,coefficient)*powers[power]
            return total
        J=[[evaluate_core(x) for x in row] for row in self.J]
        T=[[evaluate_core(x) for x in row] for row in self.T]
        U=[[evaluate_core(x) for x in row] for row in self.U]
        equations=[[(i,j,evaluate_core(c)*(rho**power if power else one)) for i,j,power,c in row]
                   for row in self.equations]
        coordinates=[sum((value*normals[j] for j,value in enumerate(row) if value),zero)
                     for row in T]
        def nonlinear(row):
            return sum((coefficient*coordinates[i]*(coordinates[j] if j>=0 else one)
                        for i,j,coefficient in equations[row] if coefficient and coordinates[i]
                        and (j<0 or coordinates[j])),zero)
        for weight in range(1,max(self.weights)+1):
            residual={j:nonlinear(row) for j,row in enumerate(self.selected)
                      if self.row_weights[row]==weight}
            for i in range(45):
                if self.weights[i]==weight:
                    coordinates[i]-=sum((U[i][j]*value for j,value in residual.items()
                                         if U[i][j] and value),zero)
        smith=[nonlinear(i)+sum((c*x for c,x in zip(J[i],coordinates) if c and x),zero)
               for i in range(135)]
        assert not any(smith[i] for i in self.selected)
        chart=[i for i in range(45) if i!=25]
        q=[zero]*45;q[25]=-one
        for i,j in enumerate(chart):q[j]=coordinates[i]
        beta=[zero]*16
        for source,value in {4:rho,11:one,25:one,26:h,50:2*one}.items():
            beta[self.one.index(source)]=value
        parameters=dict(zip(self.one,beta))|dict(zip(self.two,q))
        t=coordinates[44]
        residuals=[]
        if verify:
            for row in self.original:
                value=zero
                for mon,numerator,denominator in row:
                    term=coerce(domain,s.Rational(numerator,denominator))
                    qdegree=0
                    for index in mon:
                        term*=parameters[index]
                        qdegree+=index in self.two
                    term*=t**(2-qdegree) if qdegree!=2 else one
                    value+=term
                residuals.append(value)
            # A supplied point need not satisfy all residual equations. Their
            # values are returned, allowing callers to impose parameter ideals.
            # Vanishing here is the actual full135 condition in original coords.
        return dict(beta=beta,q=q,t=t,chart_coordinates=coordinates,
                    smith_residuals=smith,original_residuals=residuals,
                    all135_zero=verify and not any(residuals),domain=domain)

    def expression(self,value,domain,alias=None):
        ground=domain if hasattr(domain,'ext') else getattr(domain,'dom',domain)
        if alias is not None and hasattr(ground,'ext'):
            def coefficient(element):
                values=ground.convert(element).to_list()
                return sum((s.Rational(c)*alias**(len(values)-1-i)
                            for i,c in enumerate(values)),s.S.Zero)
            if domain is ground:
                return s.factor(coefficient(value))
            def polynomial(poly):
                return sum((coefficient(c)*s.prod(z**e for z,e in zip(domain.symbols,mon))
                            for mon,c in poly.items()),s.S.Zero)
            return s.factor(polynomial(value.numer)/polynomial(value.denom))
        return s.factor(domain.to_sympy(value))


def load_model(folder):
    """Create a NEW folder with frozen prepare(), or reuse its existing input."""
    folder=Path(folder)
    if not (folder/'input').exists():
        folder.mkdir(parents=True,exist_ok=True)
        prepare(folder)
    return PointLifter(folder/'input')


def verify(output):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    lifter=load_model(output/'model')
    H=s.Symbol('h');results=[]
    for name,square in [('new',-s.Rational(1127,1000)),('old-axis',-2),('old-end',-2)]:
        K=s.QQ.algebraic_field(s.sqrt(square),alias='h');h=K.unit
        v=[K.zero]*13
        if name=='new':
            v[7]=coerce(K,s.Rational(100,161))*h;v[8]=-K.one
            v[9]=coerce(K,s.Rational(10,7))*h;v[10]=K.one
        elif name=='old-axis':
            v[9]=h;v[10]=K.one
        else:
            v[12]=K.one
        point=lifter.lift(h,v,domain=K)
        assert point['all135_zero'] and not point['t']
        convert=lambda values:list(map(str,[lifter.expression(value,K,H) for value in values]))
        results.append(dict(case=name,h_square=str(square),beta=convert(point['beta']),
                            q=convert(point['q']),t=str(lifter.expression(point['t'],K,H)),
                            all135_original_equations_zero=True))
        (output/f'{name}.json').write_text(json.dumps(results[-1],indent=2)+'\n')
    # A full chart of the old weighted projective plane, before imposing k^2.
    K=s.QQ.algebraic_field(s.sqrt(-2),alias='h');sigma,k=s.symbols('s k')
    F=K.frac_field(sigma,k);h=F.convert(K.unit,K)
    sv=F.gens[0];kv=F.gens[1];v=[F.zero]*13
    v[9]=h*sv**2;v[10]=sv**2;v[12]=sv**3*kv
    point=lifter.lift(h,v,domain=F)
    assert point['all135_zero'] and not point['t']
    convert=lambda values:list(map(str,[lifter.expression(value,F,H) for value in values]))
    results.append(dict(case='old-full-plane',h_square='-2',parameters=['s','k'],
                        beta=convert(point['beta']),q=convert(point['q']),t='0',
                        all135_original_equations_zero=True))
    for name,domain,h,rho,v in [
            ('rho-retained',s.QQ,3,1,[0]*13),
            ('old-modp',s.GF(32003),3083,0,[0]*9+[3083,1,0,0])]:
        point=lifter.lift(h,v,rho=rho,domain=domain)
        assert point['all135_zero'] and not point['t']
        results.append(dict(case=name,coefficient_field=str(domain),h=str(h),rho=str(rho),
                            all135_original_equations_zero=True,t='0'))
    paths=[Path(__file__),ROOT/'code/s61_marked_sigma_pencil_model_20261102.py',
           ROOT/'code/s61_model.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',
           ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json']
    report=dict(id='S61_MARKED_SIGMA_POINT_LIFT_20261103_v62',all_passed=True,
                convention='q physical45 with q25=-1, v0..v12 from frozen Smith order, rho retained',
                cases=results,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                                            for p in paths},
                prepared_input_sha256=hashlib.sha256(lifter.input_path.read_bytes()).hexdigest(),
                seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=verify(args.output)
    print('PASS',result['id'],result['seconds'])
