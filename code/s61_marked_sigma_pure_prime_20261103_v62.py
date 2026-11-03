#!/usr/bin/env python3
"""Prime-unit spread of the known rho=1 pure-endpoint extra points at h=4."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from s61_model import S61

ROOT=Path(__file__).resolve().parents[1]
PRIME=32003


def verify():
    started=time.monotonic();model=S61();one=[i for i,w in enumerate(model.reduced['parameter_weights'])if w==1]
    two=[i for i,w in enumerate(model.reduced['parameter_weights'])if w==2]
    point_path=ROOT/'code/certificates/s61_marked_sigma_pure_endpoint_points_20261102.json'
    prior_path=ROOT/'reports/S61_MARKED_SIGMA_PURE_ENDPOINT_REMOVAL_20261102.json'
    prior=json.loads(prior_path.read_text());assert prior['all_passed']
    slice_path=ROOT/'reports/S61_MARKED_SIGMA_RHO_ONE_SPECIAL_20261103_v62/manifest.json'
    slice_evidence=json.loads(slice_path.read_text());assert slice_evidence['rho_one_lengths']['4']==6
    for name,expected in prior['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
    extra_path=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'
    source=model.reduced['equations']+[json.loads(extra_path.read_text())['universal_new_equation']]
    qs=[s.Matrix(q).applyfunc(s.Rational)for q in json.loads(point_path.read_text())['q_points']]
    beta=s.zeros(16,1)
    for variable,value in {4:1,11:1,25:1,26:4,50:2}.items():beta[one.index(variable)]=value
    def unit(value):value=s.Rational(value);return int(value.p)%PRIME and int(value.q)%PRIME
    def reduce(value):value=s.Rational(value);assert int(value.q)%PRIME;return int(value.p)*pow(int(value.q),-1,PRIME)%PRIME
    reduced_points=[];results=[]
    for index,q in enumerate(qs):
        assert q[25]==-1 and all(q[i]==0 for i in [0,19,27,32])
        assert all(int(value.q)%PRIME for value in q)
        values=dict(zip(one,beta))|dict(zip(two,q))
        highest=s.zeros(135,16);quadratic=s.zeros(16)
        for row,equation in enumerate(source):
            total=s.S.Zero
            for mon,n,d in equation:
                assert d%PRIME
                if sum(i in two for i in mon)!=2:continue
                total+=s.Rational(n,d)*s.prod(values[i]for i in mon)
                bs=[one.index(i)for i in mon if i in one]
                coefficient=s.Rational(n,d)*s.prod(q[two.index(i)]for i in mon if i in two)
                if len(bs)==1:highest[row,bs[0]]+=coefficient
                if len(bs)==2:quadratic[bs[0],bs[1]]+=coefficient
            assert total==0
        reduced=highest.applyfunc(reduce).to_DM().convert_to(s.GF(PRIME))
        pivots=list(reduced.rref()[1]);assert len(pivots)==7
        rows=list(reduced.extract(list(range(135)),pivots).transpose().rref()[1]);block=highest[rows,pivots]
        determinant=block.det(method='domain-ge');assert unit(determinant)
        free=[i for i in range(16)if i not in pivots];H=s.zeros(16,9)
        for j,i in enumerate(free):H[i,j]=1
        tail=-block.inv()*highest[rows,free]
        for j,i in enumerate(pivots):H[i,:]=tail[j,:]
        assert highest*H==s.zeros(135,9) and H.T*(quadratic+quadratic.T)*H==s.zeros(9)
        assert all(int(value.q)%PRIME for value in H)
        geometry=prior['finite_fiber_point_geometry'][index+1]
        orbit=s.Rational(geometry['orbit_minor']['determinant']);ambient=s.Rational(geometry['total_jacobian_minor']['determinant'])
        assert unit(orbit) and unit(ambient)
        assert len(geometry['orbit_minor']['columns'])==8 and len(geometry['total_jacobian_minor']['columns'])==44
        reduced_points.append(tuple(map(reduce,q)))
        results.append(dict(point_index=index,highest_rank=7,H_dimension=9,
                            highest_minor_rows=rows,highest_minor_columns=pivots,highest_minor=str(determinant),
                            kernel=[[str(c)for c in row]for row in H.tolist()],
                            all135_highest_identities_and_quadratic_kernel_verified=True,
                            orbit_minor=str(orbit),ambient_44_minor=str(ambient),
                            ambient_smooth_relative_dimension=17,ambient_t_zero_over_V=True))
    flag=tuple(int(i==25)*(-1)%PRIME for i in range(45))
    assert len({flag,*reduced_points})==3
    files=[Path(__file__),point_path,prior_path,extra_path,slice_path,ROOT/'code/s61_model.py',
           ROOT/'reports/S61_MARKED_SIGMA_STANDARD_FLAG_SPREAD_20261103_v62.json',
           ROOT/'research/notes/S61_MARKED_SIGMA_PURE_PRIME_20261103_v62.md']
    return dict(id='S61_MARKED_SIGMA_PURE_PRIME_20261103_v62',all_passed=True,prime=PRIME,
                h=4,rho=1,three_points_remain_distinct_mod_prime=True,known_lower_lengths=[4,1,1],
                supporting_slice_length_mod_prime=6,cases=results,
                scope='Two explicit h=4 extra points spread to full mixed-characteristic t=0 neighborhoods; standard flag separate.',
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in files},
                seconds=time.monotonic()-started)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=verify()
    with args.output.open('x')as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print('PASS',report['id'],report['seconds'])
