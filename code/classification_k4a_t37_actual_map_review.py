#!/usr/bin/env python3
"""Independent QQ[theta]/(theta^2+6) full actual T37 -> K4a map."""
import ast
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

ROOT = Path(__file__).resolve().parents[1]
JOB = 'research/classification/jobs/k4a_t37_actual_map_separate_review.json'
OUTPUT = 'reports/CLASSIFICATION_K4A_T37_ACTUAL_MAP_SEPARATE_REVIEW_20261108.json'


class K:
    def __init__(self, a=0, b=0):
        if isinstance(a,K): self.a,self.b=a.a,a.b
        else: self.a,self.b=F(a),F(b)
    def __add__(self,o):
        o=K(o);return K(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return K(-self.a,-self.b)
    def __sub__(self,o):return self+-K(o)
    def __rsub__(self,o):return K(o)+-self
    def __mul__(self,o):
        o=K(o);return K(self.a*o.a-6*self.b*o.b,self.a*o.b+self.b*o.a)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=K(o);n=o.a**2+6*o.b**2
        assert n
        return self*K(o.a/n,-o.b/n)
    def __rtruediv__(self,o):return K(o)/self
    def __pow__(self,n):
        assert isinstance(n,int)
        if n<0:return (1/self)**(-n)
        p=K(1)
        for _ in range(n):p=p*self
        return p
    def __bool__(self):return bool(self.a or self.b)
    def __eq__(self,o):
        o=K(o);return self.a==o.a and self.b==o.b
    def __str__(self):return str(self.a) if not self.b else f'({self.a})+({self.b})*theta'


def add(*terms):
    out={}
    for q in terms:
        for e,c in q.items():
            out[e]=out.get(e,K())+c
            if not out[e]:del out[e]
    return out


def scale(q,c):return {e:K(c)*v for e,v in q.items() if K(c)*v}


def mul(q,r):
    out={}
    for e,c in q.items():
        for f,d in r.items():
            g=tuple(x+y for x,y in zip(e,f));out[g]=out.get(g,K())+c*d
            if not out[g]:del out[g]
    return out


def power(q,n):
    out={(0,)*5:K(1)}
    for _ in range(n):out=mul(out,q)
    return out


def parse(text,names,parameters=()):
    def walk(a):
        if isinstance(a,ast.Constant):return {(0,)*5:K(a.value)} if a.value else {}
        if isinstance(a,ast.Name):
            if a.id in parameters:return {(0,)*5:K(1)}
            e=[0]*5;e[names.index(a.id)]=1;return {tuple(e):K(1)}
        if isinstance(a,ast.UnaryOp):return scale(walk(a.operand),-1) if isinstance(a.op,ast.USub)else walk(a.operand)
        p=walk(a.left)
        if isinstance(a.op,ast.Pow):return power(p,int(ast.literal_eval(a.right)))
        q=walk(a.right)
        if isinstance(a.op,ast.Add):return add(p,q)
        if isinstance(a.op,ast.Sub):return add(p,scale(q,-1))
        if isinstance(a.op,ast.Mult):return mul(p,q)
        if isinstance(a.op,ast.Div):
            assert list(q)==[(0,)*5]
            return scale(p,1/q[(0,)*5])
        raise ValueError(ast.dump(a))
    return walk(ast.parse(str(text),mode='eval').body)


def linear(matrix):
    return [{tuple(int(k==j)for k in range(5)):c for j,c in enumerate(row)if c}for row in matrix]


def subst(q,coords):
    out={}
    for e,c in q.items():
        term={(0,)*5:c}
        for i,n in enumerate(e):term=mul(term,power(coords[i],n))
        out=add(out,term)
    return out


def product(a,b):
    return [[sum((a[i][k]*b[k][j]for k in range(len(b))),K())for j in range(len(b[0]))]for i in range(len(a))]


def inverse_and_det(a):
    b=[list(row)+[K(int(i==j))for j in range(5)]for i,row in enumerate(a)];d=K(1)
    for i in range(5):
        j=next(j for j in range(i,5)if b[j][i])
        if j!=i:b[j],b[i]=b[i],b[j];d=-d
        pivot=b[i][i];d=d*pivot;b[i]=[c/pivot for c in b[i]]
        for j in range(5):
            if j!=i:
                c=b[j][i];b[j]=[s-c*t for s,t in zip(b[j],b[i])]
    return [row[5:]for row in b],d


def pair(q,i,j):
    return q.get((i,j),{})if i<j else scale(q.get((j,i),{}),-1)if j<i else {}


def push(q,m):
    return {(i,j):add(*(scale(pair(q,k,l),m[i][k]*m[j][l])for k in range(5)for l in range(5)if m[i][k]and m[j][l]))for i,j in itertools.combinations(range(5),2)}


def euler_pf(q,coords):
    out={}
    for i in range(5):
        a,b,c,d=[j for j in range(5)if j!=i]
        pf=add(mul(pair(q,a,b),pair(q,c,d)),scale(mul(pair(q,a,c),pair(q,b,d)),-1),mul(pair(q,a,d),pair(q,b,c)))
        out=add(out,scale(mul(coords[i],pf),(-1)**i))
    return out


def serialized(q):return {','.join(map(str,e)):str(c)for e,c in sorted(q.items())}


def main():
    resource.setrlimit(resource.RLIMIT_CPU,(30,31));resource.setrlimit(resource.RLIMIT_AS,(512*1024**2,)*2);begun=time.monotonic()
    job=json.loads((ROOT/JOB).read_bytes())
    for r in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/r['path']).read_bytes()).hexdigest()==r['sha256'],r['path']
    theta=K(0,1);assert theta**2==-6
    t37=json.loads((ROOT/'reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T37.json').read_bytes())
    # Public T37 coordinate names are x0,...,x4, not the map's display names.
    p={tuple(map(int,k.split(','))):parse(v,[f'x{i}'for i in range(5)],t37['period_parameters']+t37['smoothing_parameters'])for k,v in t37['homogeneous_bivector'].items()}
    pk_input=json.loads((ROOT/'reports/CLASSIFICATION_K4A_ORDER_TWO_ZERO_LOCUS_20261108.json').read_bytes())
    pk={tuple(map(int,k.split(','))):parse(v,['c','a','b','u','v'])for k,v in pk_input['generic_actual_brackets'].items()}
    m=[[K(c)for c in row]for row in [[2,-6,6,3,-3],[3,-18,18,6,-6],[0,0,0,0,0],[0,0,0,0,0],[11,-54,54,18,-18]]]
    m[2]=[theta*c for c in [0,-3,-3,2,2]];m[3]=[theta*c for c in [0,-3,-3,1,1]]
    inv,d=inverse_and_det(m);eye=[[K(int(i==j))for j in range(5)]for i in range(5)]
    assert d==2592 and product(m,inv)==eye and product(inv,m)==eye
    lam=-theta/36;assert lam and lam**2==F(-1,216)
    w=[[K(c)*theta for c in row]for row in [[0,0,F(-7,9),1,0],[0,0,-3,2,0],[3,F(-13,4),0,0,F(-7,4)],[3,F(5,4),0,0,F(-3,4)],[0,0,-5,7,0]]]
    yy=linear(m);wx=linear(product(w,m));source_push=push(p,m)
    pairs=list(itertools.combinations(range(5),2));residual={}
    for i,j in pairs:
        residual[i,j]=add(source_push[i,j],scale(subst(pk.get((i,j),{}),yy),-lam),scale(mul(yy[i],wx[j]),-1),mul(yy[j],wx[i]))
    assert all(not q for q in residual.values())
    # Repeat the same exact equality in independent target coordinates, using M^-1.
    xx=linear(inv);target_residual={}
    target_vars=linear(eye);wy=linear(w)
    for i,j in pairs:
        target_residual[i,j]=add(subst(source_push[i,j],xx),scale(pk.get((i,j),{}),-lam),scale(mul(target_vars[i],wy[j]),-1),mul(target_vars[j],wy[i]))
    assert all(not q for q in target_residual.values())
    dt=euler_pf(p,linear(eye));dk=euler_pf(pk,linear(eye));assert dt and dk
    pf_residual=add(scale(dt,d),scale(subst(dk,yy),-(lam**2)))
    assert not pf_residual
    degree2=list(itertools.combinations_with_replacement(range(5),2));degree5=list(itertools.combinations_with_replacement(range(5),5))
    assert len(pairs)*len(degree2)==150 and len(degree5)==126
    source=Path(__file__).read_bytes()
    result={'task_id':'E_reduced_nonnormal','all_independent_checks_passed':True,
      'coefficient_field':'QQ[theta]/(theta^2+6)','theta_squared':-6,'source_T37_all_period_and_smoothing_parameters':1,
      'target':'accepted K4a integercurve e1=t174, same dense24orbit',
      'source_coordinate_order':['h','x','y','z','w'],'target_coordinate_order':['c','a','b','u','v'],
      'M':[[str(c)for c in row]for row in m],'M_inverse':[[str(c)for c in row]for row in inv],
      'det_M':str(d),'both_inverse_products_exact_identity':True,'lambda':str(lam),'lambda_squared':str(lam**2),
      'Euler_W_target_matrix':[[str(c)for c in row]for row in w],
      'actual_identity':'M P_T37(X) M^t = lambda P_K4a(MX) + (MX) wedge W(MX)',
      'all150_forward_coefficient_residuals_zero':True,'all150_inverse_coordinate_residuals_zero':True,
      'Euler_Pfaffian_source':serialized(dt),'Euler_Pfaffian_target':serialized(dk),
      'all126_Pfaffian_coefficients_detM_DT_lambda2_DKMX_equal':True,
      'projective_Jacobi':'Reused accepted full350 source and completegraph zero receipts; exact invertible projective identity preserves it. No new Jacobi run.',
      'T37_actual_membership_of_dense_K4a24orbit':True,'whole_K4a_source_in_T37_via_accepted_curve_and_oneorbit_theorems':True,
      'actual24parent_equals_T37_by_orbitclosure_dimension':True,'new_component_count':0,
      'inputs':job['exact_inputs_and_hashes'],'reviewer_source_bytes':source.decode(),'reviewer_sha256':hashlib.sha256(source).hexdigest(),
      'no_producer_helper_or_CAS_import':True,'seconds':time.monotonic()-begun,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    with (ROOT/OUTPUT).open('x')as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps({k:result[k]for k in ['all_independent_checks_passed','det_M','lambda','seconds','maximum_RSS_KiB']}))


if __name__=='__main__':main()
