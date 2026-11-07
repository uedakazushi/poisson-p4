#!/usr/bin/env python3
"""Independent complete six-magnetic-source confluence to dependent3H+Q.

No producer or CAS import. Laurent arithmetic uses e and v=e+1 with
denominators cleared before every generic-chart comparison.
"""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import itertools as it
import json
from pathlib import Path
import re
import resource
import time

from classification_minimal_quartic_review import LaurentRing,determinant

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(it.combinations(range(5),2))

def power(R,p,n):
    if n<0:p=R.inverse_monomial(p);n=-n
    out=R.constant(1)
    for _ in range(n):out=R.mul(out,p)
    return out

def sub(R,p,values):
    result={}
    for ex,c in p.items():
        term=R.constant(c)
        for i,n in enumerate(ex):
            if n:term=R.mul(term,power(R,values.get(i,R.variable(R.names[i])),n))
        result=R.add(result,term)
    return result

def derivative_position(R,p,i):
    result={}
    for ex,c in p.items():
        if ex[i]:
            e=list(ex);e[i]-=1;result[tuple(e)]=c*ex[i]
    return result

def antisymmetric(R,P):
    out=dict(P)
    for(i,j),p in P.items():out[j,i]=R.scale(p,-1)
    return out

def strict_jacobi(R,P,indices):
    return {(i,j,k):R.add(*(R.mul(P.get((a,l),{}),R.derivative(P.get((b,c),{}),l))
         for a,b,c in [(i,j,k),(j,k,i),(k,i,j)]for l in indices))
         for i,j,k in it.combinations(indices,3)}

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();started=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,)*2)
    names=['reports/CLASSIFICATION_NIL22_DOUBLE_C_DEPENDENT_LIMIT_20261107.json',
      'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json',
      'reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_DEPENDENT_SOURCE_20261107.json',
      'reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_DEPENDENT_REVIEW_20261107.json',
      'research/classification/THREE_HYPERPLANES_QUADRIC_DEPENDENT.md',
      'research/classification/NIL22_DOUBLE_C_DEPENDENT_LIMIT.md']
    paths=[ROOT/p for p in names]
    claim,original,frozen,accepted=[json.loads(p.read_text())for p in paths[:4]]
    R=LaurentRing([f'm{i}'for i in range(6)]+['epsilon','v','a0','a1']+[f'r{i}'for i in range(6)],
                  ['h','q1','q2','x','y'],invert=['epsilon','v'])
    pos={v:i for i,v in enumerate(R.names)}
    def parse(text):
        # These are the only nonmonomial denominators in the focused packet.
        text=re.sub(r'\bepsilon\*\*2 \+ epsilon\b','epsilon*v',text)
        text=text.replace('(epsilon + 1)','v')
        text=re.sub(r'\bX\b','q1',text);text=re.sub(r'\bY\b','q2',text)
        text=re.sub(r'\bZ\b','x',text);text=re.sub(r'\bW\b','y',text)
        return R.parse(text)
    p=R.parse
    def cleared(poly):
        if not poly:return {}
        ex=[min(e[i]for e in poly)for i in range(len(R.names))]
        clear=R.constant(1)
        for i in [pos['epsilon'],pos['v']]:
            if ex[i]<0:clear=R.mul(clear,power(R,R.variable(R.names[i]),-ex[i]))
        return sub(R,R.mul(clear,poly),{pos['v']:p('epsilon+1')})
    def equal(a,b):assert not cleared(R.add(a,R.scale(b,-1)))
    h,q1,q2,x,y=[R.variable(v)for v in ['h','q1','q2','x','y']]
    mag=p('m0+m1*x+m2*y+m3*x**2+m4*x*y+m5*y**2')
    Q=p('x*(y+1)-1');L=p('x+1');K=p('1+epsilon*(x+1)')
    A=p('-(1+epsilon)*y-(2+3*epsilon)-epsilon*x')
    U=[{},R.scale(Q,-1)]
    V=[R.mul(L,K),p('-(1+epsilon)*(y+1)**2-(1+2*epsilon)*(y+1)-epsilon')]
    affine=antisymmetric(R,{(1,2):R.add(mag,R.mul(A,q1)),(1,4):U[1],(2,3):V[0],(2,4):V[1]})
    for ij in it.combinations(range(1,5),2):equal(affine.get(ij,{}),parse(claim['whole_affine_family'].get(','.join(map(str,ij)),'0')))
    assert all(not v for v in strict_jacobi(R,affine,range(1,5)).values())
    bracket=[R.add(*(R.mul(U[j],R.derivative(V[i],j+3))for j in range(2)),
                   R.scale(R.add(*(R.mul(V[j],R.derivative(U[i],j+3))for j in range(2))),-1))for i in range(2)]
    assert bracket==[R.mul(A,u)for u in U]
    assert R.add(R.mul(U[0],V[1]),R.scale(R.mul(U[1],V[0]),-1))==R.mul(L,K,Q)
    # Homogenize all six brackets, including zero brackets, without any m specialization.
    homogeneous={}
    for ij in PAIRS:
        out={}
        for ex,c in affine.get(ij,{}).items():
            degree=sum(ex[-4:]);assert degree<=2 and ex[R.n]==0
            e=list(ex);e[R.n]=2-degree;out[tuple(e)]=c
        homogeneous[ij]=out
        equal(out,parse(claim['whole_homogeneous_family'].get(','.join(map(str,ij)),'0')))
    homogeneous=antisymmetric(R,homogeneous)
    J=strict_jacobi(R,homogeneous,range(5))
    projective={quad:R.add(*(R.scale(R.mul(R.variable(R.names[R.n+quad[k]]),J[quad[:k]+quad[k+1:]]),(-1)**k)
                      for k in range(4)))for quad in it.combinations(range(5),4)}
    assert all(not z for z in projective.values())
    assert all(sum(ex[R.n:])==4 for z in projective.values()for ex in z)
    endpoint={ij:sub(R,z,{pos['epsilon']:{}})for ij,z in homogeneous.items()}
    source=original['cases'][2];assert source['name']=='nil22_image_double_C'
    def oldparse(text):
        return parse(re.sub(r'\bx([1-4])\b',lambda m:['q1','q2','x','y'][int(m[1])-1],text))
    endpoint_expected=antisymmetric(R,{tuple(map(int,k.split(','))):oldparse(t)for k,t in source['full_six_magnetic_source'].items()})
    assert all(endpoint.get(ij,{})==endpoint_expected.get(ij,{})for ij in PAIRS)
    assert all(not any(ex[i]<0 for i in range(R.n))for z in homogeneous.values()for ex in z)
    D={}
    for missing in range(5):
        i,j,k,l=[a for a in range(5)if a!=missing]
        pf=R.add(R.mul(homogeneous[i,j],homogeneous[k,l]),R.scale(R.mul(homogeneous[i,k],homogeneous[j,l]),-1),R.mul(homogeneous[i,l],homogeneous[j,k]))
        D=R.add(D,R.scale(R.mul(R.variable(R.names[R.n+missing]),pf),(-1)**missing))
    assert D==p('-h*(x+h)*((1+epsilon)*h+epsilon*x)*(x*y+h*x-h*h)')
    equal(D,parse(claim['actual_anticanonical_quintic']))
    # Generic actual affine map. Stored coordinates q1,q2,x,y now denote X,Y,Z,W.
    new=[p('-epsilon*(x+1)'),p('y+2+epsilon*(x+1)/v'),q2,p('-epsilon*v*q1')]
    for z,t in zip(new,claim['accepted_full_eight_source_coordinate_map']):equal(z,parse(t))
    B=[[R.derivative(z,j)for j in range(1,5)]for z in new]
    equal(determinant(R,B),p('-epsilon**2*v'))
    for row,savedrow in zip(B,claim['coordinate_jacobian']):
        for z,t in zip(row,savedrow):equal(z,parse(t))
    back=[p('-y/(epsilon*v)'),x,p('-q1/epsilon-1'),p('q2-2+q1/v')]
    backvalues={pos[v]:z for v,z in zip(['q1','q2','x','y'],back)}
    newvalues={pos[v]:z for v,z in zip(['q1','q2','x','y'],new)}
    for z,name in zip(new,['q1','q2','x','y']):equal(sub(R,z,backvalues),R.variable(name))
    for z,name in zip(back,['q1','q2','x','y']):equal(sub(R,z,newvalues),R.variable(name))
    for z,t in zip(back,claim['full_inverse_coordinate_map']):equal(z,parse(t))
    moved={}
    for i,j in it.combinations(range(4),2):
        moved[i+1,j+1]=sub(R,R.add(*(R.mul(B[i][a],affine.get((a+1,b+1),{}),B[j][b])for a,b in it.product(range(4),repeat=2))),backvalues)
    # Read/audit the immutable accepted source, without executing its SymPy code.
    knownsource=next(x['source']for x in frozen['sources']if x['path'].endswith('.py'))
    function=next(n for n in ast.parse(knownsource).body if isinstance(n,ast.FunctionDef)and n.name=='projective_source')
    expected_function='''def projective_source():
    h,x,y,z,w=xx=s.symbols('h x y z w')
    a0,a1,*rs=params=s.symbols('a0 a1 r0:6')
    f=x*(x-h)
    q=(a0*h+a1*x)*y+f
    magnetic=sum(c*m for c,m in zip(rs,(h*h,h*x,h*y,x*x,x*y,y*y)))
    P=s.zeros(5)
    for ij,value in {(1,3):f,(2,3):y*(a1*y+2*x-h),(2,4):q,(3,4):a1*y*w-magnetic}.items():
        P[ij]=s.expand(value)
        P[ij[::-1]]=-P[ij]
    return P,xx,tuple(params),f,q
'''
    assert ast.dump(function,include_attributes=False)==ast.dump(ast.parse(expected_function).body[0],include_attributes=False)
    rnew=R.scale(R.mul(p('epsilon*v'),sub(R,mag,backvalues)),-1)
    mons=[(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]
    rvalues=[]
    for exponents in mons:
        rvalues.append({ex[:R.n]+(0,)*5:c for ex,c in rnew.items()if ex[R.n:]==(0,*exponents,0,0)})
    assert all(ex[R.n:][0]==0 and ex[R.n:][3:]==(0,0)and ex[R.n:][1:3]in mons for ex in rnew)
    binding={'a0':p('epsilon*v'),'a1':p('v'),**{f'r{i}':z for i,z in enumerate(rvalues)}}
    for name,z in binding.items():equal(z,parse(claim['accepted_full_source_parameter_substitution'][name]))
    known=antisymmetric(R,{(1,3):p('q1*(q1-1)'),(2,3):p('q2*(a1*q2+2*q1-1)'),
         (2,4):p('(a0+a1*q1)*q2+q1*(q1-1)'),(3,4):p('a1*q2*y-r0-r1*q1-r2*q2-r3*q1**2-r4*q1*q2-r5*q2**2')})
    for ij in it.combinations(range(1,5),2):equal(moved[ij],sub(R,known.get(ij,{}),{pos[k]:z for k,z in binding.items()}))
    magnetic_matrix=[[derivative_position(R,z,pos[f'm{i}'])for i in range(6)]for z in rvalues]
    for row,savedrow in zip(magnetic_matrix,claim['all_six_parameter_magnetic_map']):
        for z,t in zip(row,savedrow):equal(z,parse(t))
    equal(determinant(R,magnetic_matrix),p('epsilon**2*v**6'))
    source_open=R.mul(binding['a0'],R.add(binding['a0'],binding['a1']))
    equal(source_open,p('epsilon*v**3'));equal(source_open,parse(claim['source_a0_times_a0_plus_a1']))
    assert accepted['all_passed']
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,no_producer_import=True,no_CAS_import=True,
        entire_sixmag_source_polynomial_curve_and_exact_endpoint=True,
        full_four_strict_Jacobi_and_full350_projective_Jacobi_reconstructed=True,
        actual_signed_Pfaffian_quintic_reconstructed=True,
        actual_generic_affine_map_inverse_determinant_all_six_brackets_checked=True,
        immutable_accepted_eight_source_function_AST_audited=True,
        all_eight_source_parameter_bindings_exact=True,
        all_sixmag_coefficient_map_and_determinant_checked=True,
        coordinate_determinant='-epsilon**2*(epsilon+1)',magnetic_determinant='epsilon**2*(epsilon+1)**6',
        target_source_open='a0*(a0+a1)=epsilon*(epsilon+1)**3',
        generic_parameter_chart='epsilon*(epsilon+1)!=0',source_magnetic_parameter_localizations=[],
        zero_complements='Actual family is polynomial at epsilon=0,-1 for every m. Only the generic affine recognition inverts these two curve factors. Closedness of the already accepted PGL-saturated component includes epsilon0 for each fixed arbitrary m.',
        entire_source_component_containment='R4_DEPENDENT_THREE_HYPERPLANES_QUADRIC_24',new_component_count=0,
        inherited_geometry='The accepted ENTIRE eight-source saturation is this closed 24-dimensional component, independently proved in the public dependent3H+Q review and proofnote. This run verifies literal source membership, without inferring it from divisor resemblance.',
        input_sha256={str(z.relative_to(ROOT)):hashlib.sha256(z.read_bytes()).hexdigest()for z in paths},
        reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_bytes=Path(__file__).read_text(),
        seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ['all_passed','entire_source_component_containment','seconds','maximum_RSS_KiB']}))

if __name__=='__main__':main()
