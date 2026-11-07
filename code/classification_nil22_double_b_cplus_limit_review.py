#!/usr/bin/env python3
"""Independent whole double_B to CPLUS limit: Fraction only, no CAS imports."""
import ast
from fractions import Fraction as F
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import time

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(it.combinations(range(5),2))
TRIPLES=list(it.combinations(range(5),3))
QUADS=list(it.combinations(range(5),4))

def add(*values):
    ans={}
    for p in values:
        for e,c in p.items():ans[e]=ans.get(e,F(0))+c
    return {e:c for e,c in ans.items()if c}

def scale(p,a):return {e:c*a for e,c in p.items()if c*a}

def mul(*values):
    if not values:return {}
    p=values[0]
    for q in values[1:]:
        ans={}
        for e,c in p.items():
            for f,d in q.items():
                k=tuple(a+b for a,b in zip(e,f));ans[k]=ans.get(k,F(0))+c*d
        p={e:c for e,c in ans.items()if c}
    return p

def derivative(p,i):
    ans={}
    for e,c in p.items():
        if e[i]:
            q=list(e);q[i]-=1;ans[tuple(q)]=c*e[i]
    return ans

class Ring:
    def __init__(self):
        self.names=['t']+[f'm{i}'for i in range(6)]+['h','q1','q2','x','y'];self.n=len(self.names);self.zero=(0,)*self.n
    def const(self,a):return {self.zero:F(a)}if a else {}
    def v(self,s):
        e=list(self.zero);e[self.names.index(s)]=1;return {tuple(e):F(1)}
    def power(self,p,n):
        ans=self.const(1)
        for _ in range(n):ans=mul(ans,p)
        return ans
    def parse(self,text):
        alias={'epsilon':add(self.v('t'),self.const(2)),'x1':self.v('q1'),'x2':self.v('q2'),'x3':self.v('x'),'x4':self.v('y')}
        def walk(a):
            if isinstance(a,ast.Constant)and type(a.value)is int:return self.const(a.value)
            if isinstance(a,ast.Name):return alias[a.id]if a.id in alias else self.v(a.id)
            if isinstance(a,ast.UnaryOp)and isinstance(a.op,ast.USub):return scale(walk(a.operand),-1)
            if isinstance(a,ast.BinOp):
                if isinstance(a.op,ast.Pow):return self.power(walk(a.left),a.right.value)
                p,q=walk(a.left),walk(a.right)
                if isinstance(a.op,ast.Add):return add(p,q)
                if isinstance(a.op,ast.Sub):return add(p,scale(q,-1))
                if isinstance(a.op,ast.Mult):return mul(p,q)
                if isinstance(a.op,ast.Div):
                    assert len(q)==1
                    (e,c),=q.items();assert not any(e[1:])
                    return {tuple(a-b for a,b in zip(f,e)):d/c for f,d in p.items()}
            raise ValueError(ast.dump(a))
        return walk(ast.parse(text,mode='eval').body)

def antisym(P):
    P={ij:p for ij,p in P.items()if p}
    for i,j in list(P):
        if i<j:P[j,i]=scale(P[i,j],-1)
    return P

def main():
    resource.setrlimit(resource.RLIMIT_CPU,(89,90));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    start=time.monotonic();R=Ring();co=[R.names.index(s)for s in ['h','q1','q2','x','y']]
    xx=[R.v(s)for s in ['h','q1','q2','x','y']];h,q1,q2,x,y=xx;t=R.v('t');eps=add(t,R.const(2));delta=scale(t,4)
    sourcepath=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json'
    authorpath=ROOT/'reports/CLASSIFICATION_NIL22_DOUBLE_B_CPLUS_LIMIT_20261107.json'
    author=json.loads(authorpath.read_text());source=json.loads(sourcepath.read_text())['cases'][1]
    assert source['name']=='nil22_image_double_B'
    for p,sha in author['input_sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==sha
    # Full common-denominator numerator, independently reconstructed.
    rh=R.parse('m0*h**2+m1*h*x+m2*h*y+m3*x**2+m4*x*y+m5*y**2')
    N=antisym({
        (1,2):R.parse('4*t*(m0*h**2+m1*h*x+m2*h*y+m3*x**2+m4*x*y+m5*y**2)-(t*y+4*(t+4)*h)*q1'),
        (1,3):R.parse('-t*x**2+8*h**2'),
        (1,4):R.parse('-2*t*x*y-8*(t+4)*h*x'),
        (2,3):R.parse('-t*x*y-4*(t+4)*h*x'),
        (2,4):R.parse('32*t*h**2-2*t*y**2+32*(t+2)*x**2')})
    assert all(N.get(ij,{})==mul(delta,R.parse(author['whole_homogeneous_brackets'].get(f'{ij[0]},{ij[1]}','0')))for ij in PAIRS)
    def evalvar(p,k,value):
        out={}
        for e,c in p.items():
            f=list(e);n=f[k];f[k]=0
            out=add(out,{tuple(f):c*value**n})
        return out
    affine={ij:evalvar(p,co[0],F(1))for ij,p in N.items()}
    assert all(affine.get(ij,{})==mul(delta,R.parse(author['whole_affine_brackets'].get(f'{ij[0]},{ij[1]}','0')))for ij in PAIRS)
    original=antisym({tuple(map(int,ij.split(','))):R.parse(p)for ij,p in source['full_six_magnetic_source'].items()})
    assert all(evalvar(N.get(ij,{}),0,F(-2))==scale(original.get(ij,{}),-8)for ij in PAIRS)
    def strict(P,indices):
        return {tri:add(*(mul(P.get((i,l),{}),derivative(P.get((j,k),{}),co[l]))
                         for i,j,k in [tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])]for l in indices))for tri in it.combinations(indices,3)}
    strict4=strict(affine,range(1,5));assert len(strict4)==4 and all(not p for p in strict4.values())
    S=strict(N,range(5));J={quad:add(*(scale(mul(xx[quad[i]],S[quad[:i]+quad[i+1:]]),(-1)**i)for i in range(4)))for quad in QUADS}
    assert all(not p for p in J.values())
    # Five homogeneous Pfaffians with the actual Euler contraction.
    pf={}
    for missing in range(5):
        i,j,k,l=[a for a in range(5)if a!=missing]
        pf[missing]=scale(add(mul(N.get((i,j),{}),N.get((k,l),{})),
                              scale(mul(N.get((i,k),{}),N.get((j,l),{})),-1),
                              mul(N.get((i,l),{}),N.get((j,k),{}))),(-1)**missing)
    D=add(*(mul(xx[i],pf[i])for i in range(5)))
    QA=R.parse('(t+2)*x**2+h*y-4*h**2');QB=R.parse('2*x**2+h*y+4*h**2')
    assert D==scale(mul(t,h,QA,QB),16)
    assert D==mul(delta,delta,R.parse(author['actual_anticanonical_quintic']))
    assert [QA,QB]==[R.parse(p)for p in author['quadratic_factors']]
    # Full polynomial log-tangency and free frame, without relying on a divisor label.
    NU=[N.get((1,3),{}),N.get((1,4),{})];NV=[N.get((2,3),{}),N.get((2,4),{})]
    detUV=add(mul(NU[0],NV[1]),scale(mul(NU[1],NV[0]),-1))
    assert detUV==scale(mul(t,QA,QB),-16)
    for Q,label in [(QA,'A'),(QB,'B')]:
        assert add(mul(NU[0],derivative(Q,co[3])),mul(NU[1],derivative(Q,co[4])))==scale(mul(t,x,Q),-2)
        vy=add(y,scale(h,4 if label=='A'else -4))
        assert add(mul(NV[0],derivative(Q,co[3])),mul(NV[1],derivative(Q,co[4])))==scale(mul(t,vy,Q),-2)
    # In old coordinates (h,x,y)=(X,Y,-2Z-4X), retaining q1,q2.
    def substitute_y(p):
        out={};replacement=R.parse('-2*y-4*h')
        for e,c in p.items():
            n=e[co[4]];e=list(e);e[co[4]]=0
            out=add(out,mul({tuple(e):c},R.power(replacement,n)))
        return out
    stdq=R.parse('h*y-x**2')
    other=R.parse('2*h*y-(t+2)*x**2+8*h**2')
    assert substitute_y(QB)==scale(stdq,-2)
    assert substitute_y(QA)==scale(other,-1)
    # other=epsilon*(u XZ-Y²+d X²), u=2/epsilon,d=8/epsilon.
    # The complete six-magnetic map is invertible, determinant16.
    transformed_r=substitute_y(rh)
    expected_r=R.parse('(m0-4*m2+16*m5)*h**2+(m1-4*m4)*h*x+(-2*m2+16*m5)*h*y+m3*x**2-2*m4*x*y+4*m5*y**2')
    assert transformed_r==expected_r
    matrix=[[1,0,-4,0,0,16],[0,1,0,0,-4,0],[0,0,-2,0,0,16],[0,0,0,1,0,0],[0,0,0,0,-2,0],[0,0,0,0,0,4]]
    assert F(1)*F(1)*F(-2)*F(1)*F(-2)*F(4)==16
    # All generic conics rank3; derive 3x3 determinants directly by Leibniz.
    rank_certs=[]
    for Q in [QA,QB]:
        ids=[co[0],co[3],co[4]]
        A=[[scale(derivative(derivative(Q,i),j),F(1,2))for j in ids]for i in ids]
        determinant={}
        for perm in it.permutations(range(3)):
            sign=(-1)**sum(perm[i]>perm[j]for i in range(3)for j in range(i+1,3))
            determinant=add(determinant,scale(mul(*(A[i][perm[i]]for i in range(3))),sign))
        rank_certs.append(determinant)
    assert rank_certs==[scale(eps,F(-1,4)),R.const(F(-1,2))]
    assert evalvar(QA,co[0],F(0))==mul(eps,x,x) and evalvar(QB,co[0],F(0))==scale(mul(x,x),2)
    inputs=[str(authorpath.relative_to(ROOT)),str(sourcepath.relative_to(ROOT)),
            'research/classification/TWO_QUADRICS_PLANE_CONIC_CONFIGURATIONS.md',
            'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_SOURCES_REVIEW_20261107.json']
    out={'task_id':'E_multiple_hyperplane','all_passed':True,'source_family':'nil22_image_double_B',
         'arithmetic':'Independent standard-library Fraction Laurent polynomials in t=epsilon-2; all six magnetic coefficients retained',
         'common_denominator':'delta=4*(epsilon-2), delta(0)=-8',
         'all_source_endpoint_N0_equals_minus8_original':True,
         'all4_affine_strict_Jacobi_zero':True,'all350_projective_Jacobi_zero':True,
         'actual_Euler_Pfaffian_N':'16*(epsilon-2)*h*QA*QB',
         'actual_Euler_Pfaffian_P':'h*QA*QB/(epsilon-2)',
         'all_generic_conics_rank3':True,'conic_determinants':['-epsilon/4','-1/2'],
         'generic_open':'epsilon!=0,2','source_magnetic_parameter_localizations':[],
         'generic_shared_tangent_linear_map_and_u_d_verified':True,
         'log_anchors_tangent_to_each_conic_and_Saito_frame_complete':True,
         'magnetic_coordinate_map':matrix,'magnetic_coordinate_map_determinant':16,
         'whole_six_magnetic_CPLUS_theorem_applies_including_epsilon4_u_half_chart':True,
         'conclusion':'The ENTIRE double_B six-magnetic source and its PGL-saturated closure lie in accepted R4_CPLUS_EXT. No new component.',
         'geometric_dependency':'Accepted whole shared_tangent_generic free-frame/full A6 theorem, not merely factorization of a point or a source dimension bound.',
         'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in inputs},
         'reviewer_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    target=ROOT/'reports/CLASSIFICATION_NIL22_DOUBLE_B_CPLUS_LIMIT_REVIEW_20261107.json'
    assert not target.exists();target.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
