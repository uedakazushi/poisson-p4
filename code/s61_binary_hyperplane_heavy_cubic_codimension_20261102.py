#!/usr/bin/env python3
"""Direct heavy-light cubic block; rank7 drop has codimension at least two."""
import argparse,hashlib,json,time,os,subprocess,tempfile
from pathlib import Path
from itertools import combinations,combinations_with_replacement
import sympy as s
from sympy.polys.rings import ring
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
root=Path(__file__).resolve().parents[1]
NOTE=root/'research/notes/S61_BINARY_HYPERPLANE_HEAVY_CUBIC_CODIMENSION_20261102.md'
def verify():
    begin=time.time()
    data=json.loads((root/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json').read_text())
    old=json.loads((root/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json').read_text())['linear']
    R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field();pol=lambda x:R.from_expr(s.sympify(x));frac=lambda x:K.from_expr(s.sympify(x))
    T=[[frac(row[j]) for j in data['normal34']] for row in data['T34']]
    U=[[frac(v) for v in row] for row in data['U34']]
    td=[]
    for j in range(8):
        den=R.one
        for row in T:den=den.lcm(row[j].denom)
        td.append(den)
    ud=R.one
    for row in U:
        for v in row:ud=ud.lcm(v.denom)
    Tn=[[v.numer*td[j].exquo(v.denom) for j,v in enumerate(row)] for row in T]
    Un=[[v.numer*ud.exquo(v.denom) for v in row] for row in U]
    m=S61();beta={int(k[1:])-1:s.sympify(v) for k,v in data['beta'].items()};two,F=specialize(m,beta,True)
    fixed=old['fixed'];cols=old['cols'];J=[[R.zero for j in range(45)] for i in range(135)];H=[]
    for i,f in enumerate(F):
        row=[]
        for mon,c in f.items():
            cc=pol(c)
            if len(mon)==2:
                x,y=mon
                if x!=fixed and y==fixed:J[i][cols.index(x)]+=cc
                if y!=fixed and x==fixed:J[i][cols.index(y)]+=cc
            elif len(mon)==1 and mon[0]==fixed:J[i][44]+=cc
            if fixed not in mon:
                ij=[cols.index(j) for j in mon]+[44]*(2-len(mon))
                if cc:row.append((ij,cc))
        H.append(row)


    template_path=root/'code/data/s61_binary_hyperplane_heavy_template_20261102.json'
    template=json.loads(template_path.read_text())
    assert template['variables']==['a','e','d']
    mons=template['monomials'];assert len(mons)==7
    light=[]
    for mon in mons:
        assert len(mon)==8 and mon[5]==1 and sum(mon)==2
        light.append(next(i for i in range(8) if i!=5 and mon[i]==1))
    assert sorted(light)==[0,1,2,3,4,6,7]
    known=[a,500*a+90*e**2+3*e-24,
           50*a*d+18*d*e**2-9*d*e-2*d+60*e-40,
           100*a+18*e**2-9*e-2]
    def unit(poly):
        assert poly
        exponents=[]
        for f in known:
            n=0
            while poly.rem(f)==0:poly=poly.exquo(f);n+=1
            exponents.append(n)
        assert max(map(sum,poly))==0,poly
        return str(poly[(0,0,0)]),exponents
    denominator_units=[unit(v) for v in [*td,ud]]
    zero=R.zero
    # Recheck the frozen linear section against the newly rebuilt original J.
    assert all(sum((v*Tn[j][k] for j,v in enumerate(row) if v),zero)==0
               for row in J for k in range(8))
    selected=data['selected34_rows']
    assert all(sum((v*Un[j][k] for j,v in enumerate(J[i]) if v),zero)==(ud if r==k else zero)
               for r,i in enumerate(selected) for k in range(34))
    cache={}
    def prod(i,j):
        key=tuple(sorted((i,j)))
        if key not in cache:
            cache[key]=[Tn[i][5]*Tn[j][ell]+Tn[i][ell]*Tn[j][5] for ell in light]
        return cache[key]
    Q=[]
    for terms in H:
        row=[zero]*7
        for (i,j),c in terms:
            for k,p in enumerate(prod(i,j)):
                if p:row[k]+=c*p
        Q.append(row)
    W=[[-sum((v*Q[selected[j]][k] for j,v in enumerate(row) if v),zero)
        for k in range(7)] for row in Un]
    A=[[ud*Q[i][k]+sum((v*W[j][k] for j,v in enumerate(row) if v),zero)
        for k in range(7)] for i,row in enumerate(J)]
    ids=[i for i,row in enumerate(A) if any(row)]
    assert ids==template['row_ids'] and len(ids)==16
    assert all(data['equation_weights'][i]==3 for i in ids)
    readpoly=lambda p:R.from_dict({tuple(mon):s.QQ.convert(s.Rational(c)) for mon,c in p})
    M=[[readpoly(p) for p in row] for row in template['matrix']]
    assert len(M)==16 and all(len(row)==7 for row in M)
    template_contents=[]
    for row in M:
        g=zero
        for entry in row:g=g.gcd(entry)
        assert g and max(map(sum,g))==0
        template_contents.append(str(g[(0,0,0)]))
    scales=[]
    for i,row in zip(ids,M):
        j=next(j for j,p in enumerate(row) if p)
        ratio=K(A[i][j])/K(row[j])
        assert all(K(A[i][k])==ratio*K(row[k]) for k in range(7))
        nc,ne=unit(ratio.numer);dc,de=unit(ratio.denom)
        scales.append(dict(row=i,constant=str(s.Rational(nc)/s.Rational(dc)),
                           known_factor_exponents=[x-y for x,y in zip(ne,de)]))
    change=[[s.QQ.convert(s.Rational(v)) for v in row] for row in template['constant_row_change']]
    assert len(change)==15 and all(len(row)==16 for row in change)
    assert s.Matrix(change).rank()==15
    N=[[sum((row[i]*M[i][j] for i in range(16)),zero) for j in range(7)] for row in change]
    minor_rows=template['selected_minor_rows'];assert len(minor_rows)==3
    print('DIRECT_HEAVY_BLOCK',len(ids),'terms',sum(len(p) for row in A for p in row),
          'template terms',sum(len(p) for row in M for p in row),'seconds',time.time()-begin,flush=True)
    native=root/'code/s61_binary_hyperplane_heavy_minors_20261102.cpp'
    library=root/'work/m2-runtime/lib';env=dict(os.environ,DYLD_LIBRARY_PATH=str(library)+':/opt/homebrew/lib')
    with tempfile.TemporaryDirectory(prefix='s61-hyperplane-heavy-') as temp:
        wd=Path(temp);inp=wd/'input';result=wd/'output';exe=wd/'native'
        with inp.open('x') as stream:
            stream.write('15 7\n')
            for row in N:
                for p in row:
                    stream.write(str(len(p))+'\n')
                    for mon,c in p.items():stream.write(str(c)+' '+' '.join(map(str,mon))+'\n')
            stream.write('3\n')
            for row in minor_rows:stream.write(' '.join(map(str,row))+'\n')
        flags=['clang++','-O3','-std=c++17','-I'+str(root/'work/flint-headers/flint/3.5.0/include'),
               '-I/opt/homebrew/include',str(native),str(library/'libflint.24.0.dylib'),
               '-L/opt/homebrew/lib','-lgmp','-o',str(exe)]
        compiled=subprocess.run(flags,capture_output=True,text=True)
        if compiled.returncode:raise RuntimeError(compiled.stderr)
        run=subprocess.run([str(exe),str(inp),str(result)],env=env,check=True,
                           capture_output=True,text=True,timeout=30)
        lines=result.read_text().splitlines();minors=[];gcds=[]
        def parsepoly(line):
            values=line.split();n=int(values[0]);assert len(values)==1+4*n
            return R.from_dict({tuple(map(int,values[2+4*k:5+4*k])):
                                s.QQ.convert(s.Rational(values[1+4*k])) for k in range(n)})
        for i,rows in enumerate(minor_rows):
            assert lines[4*i].split()==['MINOR',str(i),*map(str,rows)]
            minors.append(parsepoly(lines[4*i+1]));assert lines[4*i+2]==f'GCD {i}'
            gcds.append(parsepoly(lines[4*i+3]))
        expected=a**4*known[3]/100
        assert gcds[-1]==expected
        assert all(poly.rem(expected)==0 for poly in minors)
        # Small independent exact numerical determinant verification.
        for point in [(1,2,3),(2,1,1)]:
            ev=lambda p:s.Rational(p.evaluate(0,point[0]).evaluate(0,point[1]).evaluate(0,point[2]))
            for rows,poly in zip(minor_rows,minors):
                assert s.Matrix([[ev(N[i][j]) for j in range(7)] for i in rows]).det()==ev(poly)
    encode=lambda p:[[list(mon),str(c)] for mon,c in p.items()]
    paths=[Path(__file__),NOTE,native,template_path,
           root/'reports/S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102.json',
           root/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json',
           root/'code/s61_model.py',root/'code/s61_flag_transverse_obstruction_20261101.py']
    return dict(id='S61_BINARY_HYPERPLANE_HEAVY_CUBIC_CODIMENSION_20261102',all_passed=True,
        core_variables=['a','e','d'],raw_four_jet_cache_used=False,
        direct_original_equations=135,original_weight3_equation_count=sum(w==3 for w in data['equation_weights']),heavy_variable=5,light_variables_in_column_order=light,
        heavy_light_monomials=mons,nonzero_heavy_row_ids=ids,nonzero_heavy_row_count=16,
        polynomial_kernel_and_selected_right_inverse_checked=True,
        column_and_inverse_denominator_units=denominator_units,
        direct_row_to_template_unit_scales=scales,
        template_full_polynomial_proportionality_checked=True,
        template_row_contents_are_nonzero_constants=template_contents,
        constant_change_rank=15,selected_minor_rows=minor_rows,
        selected_minors_sparse=[encode(p) for p in minors],
        maximum_minor_gcd_up_to_unit=str(expected.as_expr()),
        bad_locus_dimension_upper_on_a_F_B_C_open=1,
        native_log=run.stderr,independent_scalar_points=[[1,2,3],[2,1,1]],
        limitations=['Nuisance-zero heavy-light block only; no all-point rank7 assertion.',
                    'Pure cubic rank, quartic surjectivity, and ambient saturation remain separate.'],
        proof_note=str(NOTE.relative_to(root)),
        input_sha256={str(path.relative_to(root)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
        seconds=time.time()-begin)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();assert not args.output.exists()
    result=verify()
    with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({k:result[k] for k in ['all_passed','maximum_minor_gcd_up_to_unit',
          'bad_locus_dimension_upper_on_a_F_B_C_open','seconds']}))
