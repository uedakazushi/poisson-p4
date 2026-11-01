#!/usr/bin/env python3
"""Exact B-rank-seven locus on H_plus and three non-dominant components."""
import argparse,hashlib,json,os,subprocess,tempfile,time
from itertools import combinations
from pathlib import Path
import sympy as s
from weight1_structure import sp_basis
ROOT=Path(__file__).resolve().parents[1]


def verify(m2_binary):
    start=time.time();cover_path=ROOT/'reports/S61_BINARY_RANK_DROP_COVER_20261101.json';cover=json.loads(cover_path.read_text());assert cover['all_passed']
    B=s.Matrix(cover['block_matrices'][1]).applyfunc(s.sympify);u=s.symbols('u1:62');a,b,c,d,e,f=u[3],u[13],u[34],u[35],u[48],u[49];kk=5*c+f;ll=5*b-e
    qq=600*b*b-255*b*e-63*e*e+100*a*f
    ideals=[[kk,ll],[kk,qq],
        [kk,25*b*d+25*d*e+8*f*f,5*a*d-4*b*f+2*e*f,10*b*b+5*b*e-5*e*e+4*a*f],
        [kk,50*b*d-25*d*e-4*f*f,5*a*d+6*b*f,30*b*b-15*b*e+2*a*f]]
    expr=lambda v:str(v).replace('**','^')
    mat='matrix{'+','.join('{'+','.join(expr(v) for v in row)+'}' for row in B.tolist())+'}'
    script='R=QQ[u4,u14,u35,u36,u49,u50];\nB='+mat+';\nkk=5*u35+u50;\nI=ideal(kk)+minors(8,sub(B,{u50=>-5*u35}));\n'
    for i,gens in enumerate(ideals):script+='P'+str(i)+'=ideal('+','.join(map(expr,gens))+');\n'
    script+='''J=intersect{P0,P1,P2,P3};assert isSubset(I,J);
assert all(flatten entries gens J,g->((g^2)%I==0));
assert all({P0,P1,P2,P3},isPrime);
assert (apply({P0,P1,P2,P3},P->dim(R/P))=={4,4,3,3});
assert (dim(R/(P2+ideal(u50*(u14+u49))))==2);
assert (dim(R/(P3+ideal(u50*(2*u14-u49))))==2);
print("ALL_POWERS=2;PRIME_DIMENSIONS=4,4,3,3;BOUNDARY_DIMENSIONS=2,2");
print("PASS");
'''
    env=os.environ.copy();env['DYLD_LIBRARY_PATH']=':'.join([str(ROOT/'work/m2-runtime/lib'),'/opt/homebrew/lib','/opt/homebrew/opt/givaro/lib']+([env['DYLD_LIBRARY_PATH']] if env.get('DYLD_LIBRARY_PATH') else []))
    with tempfile.TemporaryDirectory(prefix='binary-deeper-reduction-',dir=ROOT/'work') as temp:
        path=Path(temp)/'certificate.m2';path.write_text(script)
        run=subprocess.run([str(m2_binary),'--script',str(path)],env=env,capture_output=True,text=True,timeout=60)
    assert run.returncode==0,(run.stdout,run.stderr);assert 'PASS' in run.stdout
    family_path=ROOT/'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json';rank_path=ROOT/'reports/S61_BINARY_UNIFORM_RANK_20261101.json'
    family=json.loads(family_path.read_text());rank=json.loads(rank_path.read_text());uu=list(map(s.Symbol,family['sample'].keys()))
    action=s.Matrix(family['primitive_action']).applyfunc(s.sympify)
    # These five fixed group generators complement the Q-stabilizer. Thus
    # vanishing maximal minors after selecting them control ALL ten directions.
    sp=sp_basis();flat=s.Matrix.hstack(*[H.reshape(16,1) for H in sp]);rows=list(flat.T.to_DM().rref()[1]);left=flat[rows,:].inv()
    stabilizer=s.Matrix.hstack(*[left*s.Matrix(H).applyfunc(s.sympify).reshape(16,1)[rows,:] for H in rank['stabilizer_generators']])
    complement=[2,4,5,6,7];det_complement=stabilizer.row_join(s.eye(10)[:,complement]).det();assert det_complement
    parameterizations=[{u[11]:0,u[12]:0,c:-f/5,e:5*b},
        {u[11]:0,u[12]:0,c:-f/5,a:-5*(2*b-e)*(b+e)/(4*f),d:-8*f*f/(25*(b+e))},
        {u[11]:0,u[12]:0,c:-f/5,a:-15*b*(2*b-e)/(2*f),d:4*f*f/(25*(2*b-e))}]
    results=[]
    for index,param in zip([0,2,3],parameterizations):
        cons=s.Matrix([u[11],u[12]]+ideals[index]);assert cons.subs(param,simultaneous=True).applyfunc(s.cancel)==s.zeros(cons.rows,1)
        N=(cons.jacobian(uu)*action).subs(param,simultaneous=True).applyfunc(s.cancel)
        assert (N*stabilizer).applyfunc(s.cancel)==s.zeros(cons.rows,5)
        selected=N[:,complement];den=s.lcm([s.denom(v) for v in selected]);selected=selected.applyfunc(lambda v:s.cancel(v*den))
        size=4 if index==0 else 5;count=0
        for rr in combinations(range(selected.rows),size):
            for cc in combinations(range(5),size):
                M=selected[list(rr),list(cc)].to_DM();assert M.det()==M.domain.zero;count+=1
        source_dimension=12 if index==0 else 11;normal_rank_bound=size-1
        assert source_dimension+normal_rank_bound==15
        results.append(dict(ideal_index=index,dense_parameterization={str(k):str(v) for k,v in param.items()},
            source_dimension=source_dimension,group_normal_rank_upper_bound=normal_rank_bound,beta_image_dimension_upper_bound=15,
            clearing_factor=str(den),selected_normal_matrix=[[str(v) for v in row] for row in selected.tolist()],
            universally_vanishing_maximal_minor_count=count,whole_image_avoided_by_general_beta=True))
    paths=[Path(__file__),cover_path,family_path,rank_path,ROOT/'code/weight1_structure.py']
    return dict(id='S61_BINARY_DEEPER_REDUCTION_20261101',all_passed=True,field='QQ',B_rank_upper_bound=7,
        reduced_prime_ideals=[list(map(str,gens)) for gens in ideals],m2_certificate=script,m2_output=run.stdout,
        m2_binary=str(m2_binary),m2_binary_sha256=hashlib.sha256(m2_binary.read_bytes()).hexdigest(),
        complement_columns=complement,group_complement_determinant=str(det_complement),excluded_components=results,
        only_remaining_component_over_general_beta=1,
        proof_note='research/notes/S61_BINARY_DEEPER_REDUCTION_20261101.md',
        limitations=['Classifies B-rank<=7 on K=0; lower full Jacobian rank can also occur through the C block.',
            'In particular K=0,a=0,25*b-14*e=0 is a separate dominant16D candidate, not excluded here.',
            'The remaining17D component is only locally removed on its established free65 open; no endpoint exhaustion.'],
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},seconds=time.time()-start)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);p.add_argument('--m2-binary',type=Path,default=Path('/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary'));a=p.parse_args();r=verify(a.m2_binary)
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,non_dominant_components=len(r['excluded_components']),remaining_component=1,seconds=r['seconds'])))

if __name__=='__main__':main()
