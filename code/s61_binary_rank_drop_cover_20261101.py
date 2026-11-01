#!/usr/bin/env python3
"""Exact reduced full-Jacobian rank-drop locus on the nonsquare binary chart.

Uses bounded QQ ideal memberships in six variables, not sampling of rank loci.
The three orbit families are not asserted to be removed by saturation.
"""
import argparse,hashlib,json,os,subprocess,tempfile,time
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
ROOT=Path(__file__).resolve().parents[1]


def verify(m2_binary):
    start=time.time();m=S61();u=s.symbols('u1:62');weights=m.reduced['parameter_weights'];one=[i for i,w in enumerate(weights) if w==1]
    beta={i:(s.S.Zero if i in [11,12] else u[i]) for i in one};two,F=specialize(m,beta,True);fixed=two.index(29);q0=s.zeros(45,1);q0[fixed]=1
    D=s.zeros(135,45);tc=s.zeros(135,1)
    for row,f in enumerate(F):
        value=0
        for mon,c in f.items():
            if len(mon)==2:
                i,j=mon;value+=c*q0[i]*q0[j];D[row,i]+=c*q0[j];D[row,j]+=c*q0[i]
            elif len(mon)==1:tc[row]+=c*q0[mon[0]]
        assert value==0
    K0=s.Matrix.hstack(*D[:70,:].nullspace());assert K0.shape==(45,29);E=(D[70:,:]*K0).applyfunc(s.expand)
    block_path=ROOT/'reports/S61_BINARY_HYPERSURFACE_RANK_DROP_20261101.json';block=json.loads(block_path.read_text());assert block['all_passed']
    rr=block['residual_block_rows'];cc=block['residual_block_columns']
    assert K0==s.Matrix(block['constant_first70_kernel']).applyfunc(s.sympify)
    for i in range(4):
        for j in range(4):
            if i>j or i==3 or j==3 or (i,j)==(0,2):assert E[rr[i],cc[j]]==s.zeros(len(rr[i]),len(cc[j]))
    A,B,C=[E[rr[i],cc[i]] for i in range(3)]
    for i,M in enumerate([A,B,C]):assert M==s.Matrix(block['diagonal_blocks'][i]).applyfunc(s.sympify)
    # Twelve independent LINEAR rows suffice for C. A full minor here gives a
    # full minor of the actual17x5 C, including the extra quadratic G6 row.
    crow=[0,1,2,3,4,5,6,7,8,10,12,14];Clinear=C[crow,:]
    a,b,c,d,e,f=u[3],u[13],u[34],u[35],u[48],u[49];core=[a,b,c,d,e,f]
    kk=5*c+f;ll=5*b-e
    p=600*c*c+25*b*d-155*d*e+230*c*f-18*f*f
    q=600*b*b-25*a*c-255*b*e-63*e*e+95*a*f
    prime_generators=[[kk],[ll],[d,15*c-f],[5*b-3*e,a],[p,q]]
    def m2expr(x):return str(x).replace('**','^')
    def mat(M):return 'matrix{'+','.join('{'+','.join(m2expr(x) for x in row)+'}' for row in M.tolist())+'}'
    script='R=QQ['+','.join(map(str,core))+'];\nA='+mat(A)+';\nB='+mat(B)+';\nC='+mat(Clinear)+';\n'
    script+='IA=minors(4,A);IB=minors(9,B);IC=minors(5,C);\nkk='+m2expr(kk)+';ll='+m2expr(ll)+';\n'
    for i,gens in enumerate(prime_generators,1):script+='P'+str(i)+'=ideal('+','.join(map(m2expr,gens))+');\n'
    script+='''JJ=intersect{P1,P2,P3,P4,P5};
assert isSubset(IB,JJ);
assert (((kk*ll)^2)%IA==0);
assert (((kk*ll)^3)%IC==0);
assert all(flatten entries gens JJ,g->((g^2)%IB==0));
assert isPrime P5;
assert (apply({P1,P2,P3,P4,P5},P->dim(R/P))=={5,5,4,4,4});
print("A_POWER=2;C_POWER=3;B_GENERATOR_POWERS=2");
print("PRIME_DIMENSIONS=5,5,4,4,4");
print("B_RADICAL_GENERATORS="|toExternalString gens JJ);
print("PASS");
'''
    env=os.environ.copy();libs=[str(ROOT/'work/m2-runtime/lib'),'/opt/homebrew/lib','/opt/homebrew/opt/givaro/lib']
    env['DYLD_LIBRARY_PATH']=':'.join(libs+([env['DYLD_LIBRARY_PATH']] if env.get('DYLD_LIBRARY_PATH') else []))
    with tempfile.TemporaryDirectory(prefix='binary-rank-cover-',dir=ROOT/'work') as temp:
        path=Path(temp)/'certificate.m2';path.write_text(script)
        proc=subprocess.run([str(m2_binary),'--script',str(path)],env=env,capture_output=True,text=True,timeout=60)
    assert proc.returncode==0,(proc.stdout,proc.stderr)
    assert 'PASS' in proc.stdout and 'A_POWER=2;C_POWER=3;B_GENERATOR_POWERS=2' in proc.stdout
    assert 'PRIME_DIMENSIONS=5,5,4,4,4' in proc.stdout

    # Exact first70 compatibility: its radical is (kk,ll). In particular
    # B full rank implies tc is not in the leading-column span.
    left=s.Matrix.hstack(*D[:70,:].T.nullspace()).T;compat=(left*tc[:70,:]).applyfunc(s.expand)
    assert s.expand(10*compat[8]-20*compat[28]-ll**2)==0
    assert s.expand(10*compat[13]+50*compat[29]-kk**2)==0
    assert all(s.expand(v.subs({f:-5*c,e:5*b}))==0 for v in compat)

    family_path=ROOT/'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json';rank_path=ROOT/'reports/S61_BINARY_UNIFORM_RANK_20261101.json'
    family=json.loads(family_path.read_text());rank=json.loads(rank_path.read_text());uu=[u[i] for i in one]
    action=s.Matrix(family['primitive_action']).applyfunc(s.sympify);stab=s.Matrix.hstack(*[s.Matrix(v).applyfunc(s.sympify) for v in rank['primitive_actions']])
    for gens in prime_generators:
        polys=[u[11],u[12]]+gens;gb=s.groebner(polys,*uu,domain=s.QQ)
        deriv=s.Matrix(polys).jacobian(uu)*stab
        assert all(gb.reduce(s.expand(v))[1]==0 for v in deriv)
    base={s.Symbol(k):s.Rational(v) for k,v in rank['sample'].items()}
    specifications=[('hypersurface',[kk],{f:-5*base[c]},18),('linear',[d,15*c-f],{d:0,f:15*base[c]},17),
        ('quadratic',[p,q],{a:s.Rational(122,55),b:1,c:2,d:s.Rational(494,55),e:3,f:4},17)]
    strata=[]
    for name,gens,update,dimension in specifications:
        sample=base|update;cons=s.Matrix([u[11],u[12]]+gens);assert cons.subs(sample)==s.zeros(cons.rows,1)
        normal=(cons.jacobian(uu)*action).subs(sample);columns=list(normal.to_DM().rref()[1]);assert len(columns)==cons.rows
        det=normal[:,columns].det();assert det
        strata.append(dict(name=name,generators=list(map(str,gens)),sample={str(k):str(v) for k,v in sample.items()},
            beta_normal_action=[[str(v) for v in row] for row in normal.tolist()],beta_minor_columns=columns,beta_minor=str(det),
            joint_image_dimension=dimension,generic_fixed_beta_dimension=dimension-16))
    # This is the finite involution verified in the hypersurface certificate.
    transform={a:-d,b:c+f/3,c:b-e/2,d:-a,e:2*f/3,f:3*e/2}
    assert s.expand(kk.subs(transform,simultaneous=True)-ll)==0
    assert s.expand(ll.subs(transform,simultaneous=True)-kk)==0
    assert s.expand((15*c-f).subs(transform,simultaneous=True)-3*(5*b-3*e))==0
    assert d.subs(transform,simultaneous=True)==-a
    assert s.expand(p.subs(transform,simultaneous=True)-q)==0 and s.expand(q.subs(transform,simultaneous=True)-p)==0
    paths=[Path(__file__),block_path,family_path,rank_path,ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
        ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json']
    enc=lambda M:[[str(v) for v in row] for row in M.tolist()]
    return dict(id='S61_BINARY_RANK_DROP_COVER_20261101',all_passed=True,field='QQ',core_variables=list(map(str,core)),
        block_matrices=[enc(A),enc(B),enc(Clinear)],prime_ideals=[list(map(str,g)) for g in prime_generators],
        membership_powers=dict(A_KL=2,C_KL=3,B_radical_generators=2),m2_certificate=script,m2_output=proc.stdout,
        m2_binary=str(m2_binary),m2_binary_sha256=hashlib.sha256(m2_binary.read_bytes()).hexdigest(),
        first70_t_compatibility=list(map(str,compat)),t_compatibility_radical=[str(kk),str(ll)],
        full_rank_drop_equivalence='rank(full chart Jacobian)<=34 iff rank(B)<=8 iff beta lies in the union of the five listed prime zero sets.',
        finite_symmetry_core_images={str(k):str(v) for k,v in transform.items()},orbit_strata=strata,
        proof_note='research/notes/S61_BINARY_RANK_DROP_COVER_20261101.md',
        limitations=['Classifies Jacobian rank-drop only on the nonsquare binary orbit, not its higher obstructions or saturated support.',
            'The18D and two17D orbit strata can intersect; their full generic-beta fibers are not removed here.',
            'No full affine-fiber finiteness or component exhaustion.'],
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},seconds=time.time()-start)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--m2-binary',type=Path,default=Path('/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary'))
    args=parser.parse_args();report=verify(args.m2_binary)
    with args.output.open('x') as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_passed=True,strata=[{k:c[k] for k in ['name','joint_image_dimension','generic_fixed_beta_dimension','beta_minor']} for c in report['orbit_strata']],seconds=report['seconds'])))

if __name__=='__main__':main()
