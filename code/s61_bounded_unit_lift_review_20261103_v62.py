#!/usr/bin/env python3
"""Independent Fraction scalar lift and selected coefficient-matrix audit."""
import argparse,hashlib,itertools,json,shlex,subprocess,tarfile,sys
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1]
STEM='S61_BINARY_NONSQUARE_BOUNDED_UNIT_LIFT_20261103_v62'

def polynomial(line):
    t=line.split();n=int(t[0]);assert len(t)==1+4*n
    return {tuple(map(int,t[2+4*k:5+4*k])):F(t[1+4*k])for k in range(n)}
def accumulate(a,b,c=F(1)):
    for m,v in b.items():
        z=a.get(m,0)+c*v
        if z:a[m]=z
        elif m in a:del a[m]
def product(a,b):
    out={}
    for m,c in a.items():
        for n,d in b.items():accumulate(out,{tuple(x+y for x,y in zip(m,n)):c*d})
    return out

def main(out):
    out.mkdir(parents=True,exist_ok=False);folder=ROOT/'reports'/STEM;d=json.loads((folder/'report.json').read_text())
    with tarfile.open(folder/'evidence.tar.gz')as t:
        saved={m.name:t.extractfile(m).read()for m in t.getmembers()}
    for n,v in d['evidence'].items():assert len(saved[n])==v['bytes']and hashlib.sha256(saved[n]).hexdigest()==v['sha256']
    lines=saved['source14.txt'].decode().splitlines();assert lines[0]=='14 46 10'
    assert len(lines[1].split())==10 and len(lines[2].split())==14
    A=[[polynomial(lines[3+46*i+j])for j in range(46)]for i in range(14)]
    z=saved['scalar.witness'].decode().splitlines();assert z[0]=='46';co=[polynomial(s)for s in z[1:]];assert len(co)==46
    total={}
    for j in range(46):accumulate(total,product(A[1][j],co[j]))
    target={(2,0,1):F(5),(1,1,0):F(-6),(1,0,0):F(-4)}
    assert total==target and sum(map(len,co))==154 and all(sum(m)<=3 for p in co for m in p)
    assert set(co[22])=={(0,0,0)}and co[22][(0,0,0)]!=0
    replacement_determinant=co[22][(0,0,0)]
    vectors=[[A[j][i]for j in range(14)]for i in range(46)]
    vectors.append([({m:-c for m,c in target.items()}if j==1 else{})for j in range(14)])
    mus=[(a,b,n-a-b)for n in range(6)for a in range(n+1)for b in range(n-a+1)]
    assert len(mus)==56;bounds=[mus]*46+[[(0,0,0)]];offset=[0]
    for r in bounds:offset.append(offset[-1]+len(r))
    assert offset[-1]==2577
    equations={}
    for j in range(14):
        for i in range(47):
            for m in sorted(vectors[i][j],reverse=True):
                for u in bounds[i]:
                    key=(j,*(x+y for x,y in zip(m,u)))
                    if key not in equations:equations[key]=len(equations)
    assert len(equations)==11302
    selection=list(map(int,saved['isolated.kernel.selected-equations'].decode().split()))
    assert len(selection)==len(set(selection))==2577
    selected={v:i for i,v in enumerate(selection)};entries=[];prime=32003
    for j in range(14):
        for i in range(47):
            for m,c in vectors[i][j].items():
                assert c.denominator%prime
                value=c.numerator*pow(c.denominator,-1,prime)%prime
                for k,u in enumerate(bounds[i]):
                    key=(j,*(x+y for x,y in zip(m,u)));row=selected.get(equations[key])
                    if row is not None and value:entries.append((row,offset[i]+k,value))
    inp=out/'selected-matrix.in';inp.write_text(f'2577 {prime} {len(entries)}\n'+''.join(f'{i} {j} {c}\n'for i,j,c in entries))
    cpp=ROOT/'code/nmod_selected_rank_review_20261103_v62.cpp';binary=out/'rank'
    flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True))
    subprocess.run(['c++','-O2',str(cpp),'-o',str(binary),*flags],check=True)
    r=subprocess.run([str(binary),str(inp)],capture_output=True,text=True,timeout=90)
    (out/'rank.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stdout+r.stderr
    fresh=ROOT/'work/old-plane-v62/unit-lift-review-fresh/report.json'
    f=json.loads(fresh.read_text());assert f['all_passed']and f['residual_sha256']==d['residual_sha256']
    paths=[Path(__file__),cpp,folder/'report.json',folder/'evidence.tar.gz',
        ROOT/'research/notes'/f'{STEM}.md',ROOT/'research/notes'/f'{STEM}_REVIEW.md']
    for n,want in d['input_sha256'].items():
        q=ROOT/n;assert hashlib.sha256(q.read_bytes()).hexdigest()==want;paths.append(q)
    report=dict(id=STEM+'_REVIEW',all_passed=True,independent_all_archive_hashes=True,
        independent_Fraction_degree3_scalar_identity=True,independent_constant_unit_column_replacement22=True,
        replacement_determinant=str(replacement_determinant),
        independent_full_finite_coefficient_matrix=[11302,2577],independent_selected_square_rank=2577,
        independent_rank_prime=prime,QQ_no_isolated_column_lift_degree_at_most5=True,
        fresh_full_residual13_sha_matches=True,all_point_rank_or_ambient_exhaustion_proved=False,
        input_sha256={str(q.relative_to(ROOT)):hashlib.sha256(q.read_bytes()).hexdigest()for q in paths},
        evidence_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in [inp,out/'rank.log']})
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('Bounded unit lift independent coefficient/rank review PASS')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
