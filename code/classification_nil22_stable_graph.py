#!/usr/bin/env python3
"""One exact T39/T40 graph test for both entire stable nil22 sources."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

import sympy as s
from classification_nil22_image_source_review import build_P,coords,mag,x,y
from s61_model import S61,X,pushforward

ROOT=Path(__file__).resolve().parents[1]


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();begun=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(170,175));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    model=S61();k0,k1,k2=s.symbols('k0 k1 k2');specs=[('stable_line_conic',s.Matrix([[3,-15],[-15,3]]),s.Matrix([s.Rational(3,2),-s.Rational(3,2)]),[0,0,0,0,k0,k1],[k0,k1]),
        ('stable_three_lines',s.Matrix([[-6,3],[3,-6]]),s.zeros(2,1),[0,0,0,k0,k1,k2],[k0,k1,k2])];records=[]
    for name,S,V0,coefficients,parameters in specs:
        affine,_,_=build_P(x*y*(x-y),S,V0,-s.eye(2)/3);affine=affine.subs(dict(zip(mag,coefficients)),simultaneous=True)
        P=s.zeros(5)
        for i,j in itertools.combinations(range(4),2):
            entry=sum(c*X[0]**(2-sum(e))*s.prod(z**a for z,a in zip(X[1:],e))for e,c in s.Poly(affine[i,j],*coords).terms());P[i+1,j+1]=entry;P[j+1,i+1]=-entry
        K=P.subs(dict(zip(X,[1,0,0,0,0])))[1:,1:];pf=s.factor(K[0,1]*K[2,3]-K[0,2]*K[1,3]+K[0,3]*K[1,2]);assert pf!=0 and not pf.free_symbols
        unit=[s.eye(4)[:,i]for i in range(4)];r0=unit[0];r1=unit[2]/K[0,2];other=[]
        for i in (1,3):
            r=unit[i];other.append(r-(r.T*K*r1)[0]*r0+(r.T*K*r0)[0]*r1)
        r2,r3=other;r3=r3/s.cancel((r2.T*K*r3)[0]);frame=s.eye(5);frame[1:,1:]=s.Matrix.vstack(*(r.T for r in [r0,r1,r2,r3]))
        J=s.zeros(4);J[0,1]=J[2,3]=1;J[1,0]=J[3,2]=-1;assert(frame[1:,1:]*K*frame[1:,1:].T-J).applyfunc(s.cancel)==s.zeros(4)
        assert all(not value.free_symbols for value in frame)
        normalized=model.normalize(pushforward(P,frame));verification=model.verify_point(normalized,check_equations=True);point=normalized['parameters'];us=s.symbols('u1:62')
        assert all(not s.denom(s.cancel(value)).free_symbols for value in [*point,*normalized['perspective']])
        graphs={}
        for target,path in [('T39','reports/Q51_UNIVERSAL_SECTION_20261101.json'),('T40','reports/Q52_UNIVERSAL_SECTION_20261101.json')]:
            packet=json.loads((ROOT/path).read_text());indices=packet['weight2_u_indices_zero_based']
            residuals=[s.factor(s.cancel(point[i]-s.sympify(expression).subs(dict(zip(us,point)),simultaneous=True)))for i,expression in zip(indices,packet['quadratic_section_q'])]
            constant=[dict(u_index=i+1,residual=str(v))for i,v in zip(indices,residuals)if v!=0 and not v.free_symbols]
            graphs[target]=dict(all45_residuals=list(map(str,residuals)),whole_source_in_graph=all(v==0 for v in residuals),
                entire_source_excluded_by_nonzero_constant=bool(constant),nonzero_constant_residuals=constant,
                scope='A nonzero constant excludes every residual parameter, including its zero loci. Without that constant, only noncontainment of the whole source is concluded from a nonzero polynomial.')
        records.append(dict(name=name,all_residual_parameters=list(map(str,parameters)),canonical_F0=str(affine[0,1].subs(dict.fromkeys(coords[:2],0))),
            fixed_Darboux_frame=frame.tolist(),origin_Pf=str(pf),no_parameter_localization=True,all150_full134_and350_checks=verification,
            perspective=list(map(str,normalized['perspective'])),all61_parameters=list(map(str,point)),exact_graph_membership=graphs,
            whole_original_six_source_covered=any(value['whole_source_in_graph']for value in graphs.values())))
    paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/classification_nil22_image_source_review.py',ROOT/'reports/CLASSIFICATION_NIL22_STABLE_PROJECTIVE_ACTION_V2_20261107.json',ROOT/'reports/Q51_UNIVERSAL_SECTION_20261101.json',ROOT/'reports/Q52_UNIVERSAL_SECTION_20261101.json',*[ROOT/('legacy/v02/data/'+name)for name in ['symplectic_slice_certificate.json','symplectic_slice_equations.json','symplectic_slice_reduced_equations.json']]]
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,cases=records,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_source_bytes=Path(__file__).read_text(),
        scope='One exact whole-source T39/T40 graph recognition; no other MPS exclusion or componenthood follows.',seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as handle:json.dump(result,handle,indent=2,default=str);handle.write('\n')
    print(json.dumps(dict(all_passed=True,covered=[r['whole_original_six_source_covered']for r in records],seconds=result['seconds'])))


if __name__=='__main__':main()
