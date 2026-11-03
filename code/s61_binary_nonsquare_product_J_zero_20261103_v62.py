#!/usr/bin/env python3
"""Exact QQ: product-rank defects on J=5ad-6e-4=0 have dimension at most zero."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import random
import shutil
import subprocess
import tarfile
import time

import numpy as np
import sympy as s
from sympy.polys.rings import ring

from s61_binary_nonsquare_mod_ell_coordinates_20261103_v62 import load_chart
from s61_binary_deeper_uniform_quadratic_20261101 import elimination_rank

ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(output,reuse):
    begin=time.time();output.mkdir(parents=True,exist_ok=False)
    c=load_chart();R=c['R'];a,e,d=R.gens
    R2,aa,ee=ring('a,e',s.QQ);K2=R2.to_field();ak,ek=K2(aa),K2(ee)
    dk=(6*ek+4)/(5*ak)
    def cv(poly):
        return sum((K2(co)*ak**ex[0]*ek**ex[1]*dk**ex[2] for ex,co in poly.items()),K2.zero)
    P,*x=ring('x0,x1,x2,x3,x4,x5,x6',K2)
    images=[x[0],3*x[0]+x[3]-x[5],x[1],x[6],x[2],P.zero,x[3],x[4]]
    def move(poly):
        result=P.zero
        for mon,co in poly.items():
            term=P(cv(co))
            for n,img in zip(mon,images):
                if n:term*=img**n
            result+=term
        return result
    ell=move(c['ell']);quad=[move(c['old_quadrics'][i])for i in [0,1,3,5]]
    coeff=[ell.get(tuple(int(i==j)for j in range(7)),K2.zero)for i in range(7)]
    pivot=2;g=coeff[pivot]
    assert g==-s.Rational(32,27)*(3*ek+2)
    bs=300*aa*ee+54*ee**3+9*ee**2-24*ee-4
    units=[aa,500*aa+90*ee**2+3*ee-24,bs,100*aa+18*ee**2-9*ee-2,3*ee+2]
    assert cv(c['B'])==2*K2(bs)/(5*ak)
    assert dk==2*(3*ek+2)/(5*ak)
    def is_unit(poly):
        if not poly:return False
        for f in units:
            while not poly.rem(f):poly=poly.exquo(f)
        return poly.is_ground
    tail=sum((P(coeff[j])*x[j]for j in range(7)if j!=pivot),P.zero)
    cleared=[]
    for q in quad:
        result=P.zero
        for mon,co in q.items():
            n=mon[pivot];rest=list(mon);rest[pivot]=0
            result+=P.from_dict({tuple(rest):co})*((-tail/g)**n)
        den=R2.one
        for co in result.values():den=den.lcm(co.denom)
        assert is_unit(den)
        row={mon:(co*den).numer.exquo((co*den).denom)for mon,co in result.items()if co}
        for f in units:
            while all(not co.rem(f)for co in row.values()):row={mon:co.exquo(f)for mon,co in row.items()}
        cleared.append(row)
    positions=[j for j in range(7)if j!=pivot]
    monomials=sorted(tuple(triple.count(j)for j in range(7))for triple in itertools.combinations_with_replacement(positions,3))
    rows=[{tuple(n+(k==j)for k,n in enumerate(mon)):co for mon,co in q.items()}for q in cleared for j in positions]
    assert(len(rows),len(monomials))==(24,56)
    matrix=output/'matrix.txt';unitpath=output/'units.txt'
    with matrix.open('w')as f:
        f.write('24 56 0\n')
        for row in rows:
            for mon in monomials:
                p=row.get(mon,R2.zero);f.write(str(len(p))+''.join(' '+str(co)+' '+' '.join(map(str,ex))+' 0'for ex,co in p.items())+'\n')
    unitpath.write_text('\n'.join(str(f).replace('**','^')for f in units)+'\n')
    nativepath=ROOT/'code/polynomial_unit_span_pivots_20261103_v62.cpp'
    flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
    executable=output/'unit-span'
    subprocess.run(['c++','-w','-O3','-std=c++17',str(nativepath),*flags,'-lgmp','-o',str(executable)],check=True,capture_output=True,text=True)
    residual=output/'residual.txt'
    cp=subprocess.run([str(executable),str(matrix),str(residual),str(unitpath),'0'],check=True,capture_output=True,text=True,timeout=60)
    (output/'unit-span.log').write_text(cp.stderr)
    with residual.open()as f:
        nr,nc,npiv=map(int,f.readline().split());f.readline();f.readline();data=[]
        for _ in range(nr*nc):
            z=f.readline().split();nt=int(z[0]);data.append(R2.from_dict({tuple(map(int,z[2+4*t:4+4*t])):s.QQ(z[1+4*t])for t in range(nt)}))
    assert(nr,nc,npiv)==(12,44,12)
    prime=32003
    def mod(co):
        co=s.QQ(co);return int(co.numerator)*pow(int(co.denominator),-1,prime)%prime
    ev=lambda p:sum(mod(co)*pow(2,ex[1],prime)for ex,co in p.items())%prime
    A=np.array([ev(p)for p in data],dtype=np.int64).reshape(nr,nc)
    rng=random.Random(62711);choices=[]
    for _ in range(4):
        perm=list(range(nc));rng.shuffle(perm);rank,pivots=elimination_rank(A[:,perm],prime);assert rank==12
        choices.append([perm[j]for j in pivots])
    inp=output/'minor-input.txt'
    with inp.open('w')as f:
        f.write('44 12\n')
        for j in range(nc):
            for i in range(nr):
                p=data[i*nc+j];f.write(str(len(p))+''.join(' '+str(co)+' '+' '.join(map(str,ex))+' 0'for ex,co in p.items())+'\n')
        f.write('4\n')
        for choice in choices:f.write(' '.join(map(str,choice))+'\n')
    determinant_source=ROOT/'code/s61_binary_hyperplane_heavy_minors_20261102.cpp'
    source=determinant_source.read_text().replace('assert(m==7)','assert(m==12)').replace('vector<int>ids(7)','vector<int>ids(m)')
    (output/'minors.cpp').write_text(source)
    if reuse:
        assert inp.read_bytes()==(reuse/'input').read_bytes()
        raw=(reuse/'output').read_bytes();log=(reuse.parent/'jzero-minors.log').read_text()
        assert 'DONE' in log
        provenance=str(reuse)
    else:
        exe=output/'minors'
        subprocess.run(['c++','-w','-O3','-std=c++17',str(output/'minors.cpp'),*flags,'-lgmp','-o',str(exe)],check=True,capture_output=True,text=True)
        cp=subprocess.run([str(exe),str(inp),str(output/'minor-output.txt')],check=True,capture_output=True,text=True,timeout=100)
        raw=(output/'minor-output.txt').read_bytes();(output/'minor-output.txt').unlink();log=cp.stderr;exe.unlink();provenance='executed by this invocation'
    (output/'minors.log').write_text(log)
    with gzip.open(output/'minor-output.txt.gz','wb')as f:f.write(raw)
    lines=raw.decode().splitlines()
    def readpoly(line):
        z=line.split();nt=int(z[0]);assert len(z)==1+4*nt
        assert all(int(z[4+4*t])==0 for t in range(nt))
        return R2.from_dict({tuple(map(int,z[2+4*t:4+4*t])):s.QQ(z[1+4*t])for t in range(nt)})
    polynomials=[readpoly(lines[4*k+1])for k in range(4)]
    gcd=readpoly(lines[15]);orders=[]
    for f in units:
        n=0
        while not gcd.rem(f):gcd=gcd.exquo(f);n+=1
        orders.append(n)
    assert gcd and gcd.is_ground
    # Independent exact numerical determinant checks of all selected QQ minors.
    for point in [(1,2),(2,1)]:
        evqq=lambda p:s.Rational(p.evaluate(0,point[0]).evaluate(0,point[1]))
        for choice,p in zip(choices,polynomials):
            assert s.Matrix([[evqq(data[i*nc+j])for j in choice]for i in range(nr)]).det()==evqq(p)
    executable.unlink()
    paths=[Path(__file__),nativepath,determinant_source,ROOT/'code/s61_binary_nonsquare_mod_ell_coordinates_20261103_v62.py',*c['input_paths']]
    report=dict(all_passed=True,field='QQ',locus='J=5ad-6e-4=0, a*F*B*C*d!=0',
                d_substitution='(6e+4)/(5a)',ell_pivot=str(g),known_units=list(map(str,units)),
                original_product_matrix_shape=[24,56],known_unit_pivots=12,residual_shape=[12,44],
                selected_residual_columns=choices,exact_maximal_minor_gcd_known_unit_orders=orders,
                exact_gcd_remaining_constant=str(gcd),independent_QQ_determinant_points=[[1,2],[2,1]],
                conclusion='On the stated J=0 locus, the rank24 defect set has dimension<=0. By the pair lemma the rank41 defect set of all six quadrics light multiples also has dimension<=0.',
                all_core_product_rank41_proved=False,ambient_nilpotence_claimed=False,
                determinant_calculation_reused=provenance,
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
                generated_sha256={p.name:sha(p)for p in output.iterdir()},seconds=time.time()-begin)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k]for k in ['all_passed','conclusion','seconds']}))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--reuse-minors',type=Path)
    args=ap.parse_args();run(args.output,args.reuse_minors)
