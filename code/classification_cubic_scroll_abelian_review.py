#!/usr/bin/env python3
"""Independent whole-source smooth-scroll dimension certificate replay.

Standard-library Fraction/sparse arithmetic only. No producer/helper/CAS
imports and no Groebner calculation. Cone assertions are not reviewed.
"""
import argparse,ast,hashlib,itertools as it,json,math,resource,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(it.combinations(range(5),2));TRIPLES=list(it.combinations(range(5),3));QUADS=list(it.combinations(range(5),4))


def add(*values):
    out={}
    for v in values:
        for e,c in v.items():out[e]=out.get(e,Q(0))+c
    return {e:c for e,c in out.items()if c}


def scale(v,c):return {e:a*c for e,a in v.items()if a*c}


def mul(a,b):
    out={}
    for e,c in a.items():
        for f,d in b.items():
            g=tuple(x+y for x,y in zip(e,f));out[g]=out.get(g,Q(0))+c*d
    return {e:c for e,c in out.items()if c}


def diff(v,i):
    out={}
    for e,c in v.items():
        if e[i]:
            f=list(e);f[i]-=1;out[tuple(f)]=c*e[i]
    return out


def shift(v,i):
    out={}
    for e,c in v.items():f=list(e);f[i]+=1;out[tuple(f)]=c
    return out


def monomials(n,d):
    return sorted([tuple(a.count(i)for i in range(n))for a in it.combinations_with_replacement(range(n),d)],reverse=True)


def parse(text,names):
    zero=(0,)*len(names)
    def walk(node):
        if isinstance(node,ast.Constant)and type(node.value)is int:return {zero:Q(node.value)}if node.value else{}
        if isinstance(node,ast.Name):
            e=list(zero);e[names.index(node.id)]=1;return {tuple(e):Q(1)}
        if isinstance(node,ast.UnaryOp)and isinstance(node.op,ast.USub):return scale(walk(node.operand),-1)
        if isinstance(node,ast.BinOp):
            a=walk(node.left)
            if isinstance(node.op,ast.Pow):
                assert isinstance(node.right,ast.Constant)and type(node.right.value)is int and node.right.value>=0
                out={zero:Q(1)}
                for _ in range(node.right.value):out=mul(out,a)
                return out
            b=walk(node.right)
            if isinstance(node.op,ast.Add):return add(a,b)
            if isinstance(node.op,ast.Sub):return add(a,scale(b,-1))
            if isinstance(node.op,ast.Mult):return mul(a,b)
            if isinstance(node.op,ast.Div):assert set(b)=={zero};return scale(a,1/b[zero])
        raise ValueError(ast.dump(node))
    return walk(ast.parse(text,mode='eval').body)


def rref(rows):
    basis={}
    for row in rows:
        v={j:Q(c)for j,c in(row.items()if isinstance(row,dict)else enumerate(row))if c}
        for p in sorted(basis):
            if p in v:v=add(v,scale(basis[p],-v[p]))
        if v:
            p=min(v);basis[p]=scale(v,1/v[p])
    for p in sorted(basis,reverse=True):
        for q in sorted(basis):
            if q<p and p in basis[q]:basis[q]=add(basis[q],scale(basis[p],-basis[q][p]))
    return basis


def determinant(matrix):
    a=[[Q(c)for c in row]for row in matrix];answer=Q(1)
    for i in range(len(a)):
        k=next((k for k in range(i,len(a))if a[k][i]),None)
        if k is None:return Q(0)
        if k!=i:a[k],a[i]=a[i],a[k];answer=-answer
        c=a[i][i];answer*=c
        for j in range(i+1,len(a)):
            d=a[j][i]/c
            for l in range(i+1,len(a)):a[j][l]-=d*a[i][l]
            a[j][i]=Q(0)
    return answer


def bivectors(quadrics,B):
    bases=[]
    for k in range(len(B[0])):
        P=[[{}for j in range(5)]for i in range(5)]
        for t,(i,j)in enumerate(PAIRS):
            P[i][j]=add(*(scale(quadrics[g],B[10*g+t][k])for g in range(3)));P[j][i]=scale(P[i][j],-1)
        bases.append(P)
    return bases


def jacobi_and_divisor(bases):
    n=len(bases);columns=list(it.combinations_with_replacement(range(n),2));mons4=monomials(5,4);mons5=monomials(5,5)
    loc4={e:i for i,e in enumerate(mons4)};loc5={e:i for i,e in enumerate(mons5)};J=[{}for _ in range(350)];D=[{}for _ in range(126)]
    deriv=[[[[diff(P[i][j],l)for l in range(5)]for j in range(5)]for i in range(5)]for P in bases]
    for col,(a,b)in enumerate(columns):
        orientations=[(a,b)]if a==b else[(a,b),(b,a)];strict={tri:{}for tri in TRIPLES}
        for tri in TRIPLES:
            for i,j,k in [tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])]:
                for l in range(5):
                    strict[tri]=add(strict[tri],*(mul(bases[c][i][l],deriv[d][j][k][l])for c,d in orientations))
        for row,quad in enumerate(QUADS):
            P=add(*(scale(shift(strict[quad[:i]+quad[i+1:]],missing),(-1)**i)for i,missing in enumerate(quad)))
            for e,c in P.items():J[70*row+loc4[e]][col]=c
        total={}
        for missing in range(5):
            i,j,k,l=[v for v in range(5)if v!=missing]
            for first,second,sign in [((i,j),(k,l),1),((i,k),(j,l),-1),((i,l),(j,k),1)]:
                value=add(*(mul(bases[c][first[0]][first[1]],bases[d][second[0]][second[1]])for c,d in orientations))
                total=add(total,scale(shift(value,missing),(-1)**missing*sign))
        for e,c in total.items():D[loc5[e]][col]=c
    return J,D


def decode(rows):return [{int(j):Q(c)for j,c in row.items()}for row in rows]


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();begun=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(170,175));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    oldpath=ROOT/'reports/CLASSIFICATION_CUBIC_SCROLL_ZERO_JACOBI_20261107.json';sourcepath=ROOT/'reports/CLASSIFICATION_CUBIC_SCROLL_ABELIAN_CONORMAL_20261107.json';dimpath=ROOT/'reports/CLASSIFICATION_CUBIC_SCROLL_ABELIAN_DIMENSION_20261107.json'
    old=json.loads(oldpath.read_text());source=json.loads(sourcepath.read_text());dimensions=json.loads(dimpath.read_text());oldcase=next(v for v in old['cases']if v['name']=='smooth_scroll_12');case=next(v for v in source['cases']if v['name']=='smooth_scroll_12');dim=next(v for v in dimensions['cases']if v['name']=='smooth_scroll_12')
    assert source['input_sha256']==hashlib.sha256(oldpath.read_bytes()).hexdigest();assert dimensions['input_sha256']==hashlib.sha256(sourcepath.read_bytes()).hexdigest()
    quadrics=[parse(t,[f'x{i}'for i in range(5)])for t in ['x0*x3-x1*x2','x0*x4-x1*x3','x2*x4-x3**2']]
    mapping=[(1,0,1,0),(0,1,1,0),(2,0,0,1),(1,1,0,1),(0,2,0,1)]
    def pull(poly):
        out={}
        for e,c in poly.items():
            f=tuple(sum(e[i]*mapping[i][j]for i in range(5))for j in range(4));out[f]=out.get(f,Q(0))+c
        return {e:c for e,c in out.items()if c}
    assert all(not pull(q)for q in quadrics)
    derivatives=[[diff(q,i)for i in range(5)]for q in quadrics];rows=[];labels=[]
    for a,b in it.combinations(range(3),2):
        for k in range(5):
            columns=[]
            for g in range(3):
                for i,j in PAIRS:columns.append(pull(mul(derivatives[g][k],add(mul(derivatives[a][i],derivatives[b][j]),scale(mul(derivatives[a][j],derivatives[b][i]),-1)))))
            for e in sorted(set().union(*(set(v)for v in columns)),reverse=True):
                rows.append([v.get(e,Q(0))for v in columns]);labels.append(dict(generator_pair=[a,b],ambient_partial=k,parametrization_exponents=list(e)))
    assert len(rows)==159 and rows==[[Q(c)for c in row]for row in case['full_linear_matrix']] and labels==case['all_generator_pair_all_firstpartial_equations']
    RR=rref(rows);assert len(RR)==18
    free=[j for j in range(30)if j not in RR];B=[[Q(0)for _ in free]for j in range(30)]
    for k,j in enumerate(free):
        B[j][k]=1
        for p in RR:B[p][k]=-RR[p].get(j,Q(0))
    assert free==case['free_original30_indices'];assert B==[[Q(v)for v in row]for row in case['complete_original30_by_reduced_source_matrix']]
    assert all(sum(c*B[j][k]for j,c in enumerate(row))==0 for row in rows for k in range(12))
    cert=case['rank_minor'];det=determinant([[rows[i][j]for j in cert['columns']]for i in cert['rows']]);assert det==Q(cert['determinant'])and det
    identity=[[Q(i==j)for j in range(30)]for i in range(30)];J30,D30=jacobi_and_divisor(bivectors(quadrics,identity));assert J30==decode(oldcase['all350_Jacobi_rows']);assert D30==decode(oldcase['all126_actual_D_rows'])
    J,D=jacobi_and_divisor(bivectors(quadrics,B));assert J==decode(case['all350_reduced_Jacobi_rows']);assert D==decode(case['all126_reduced_actual_D_rows']);assert len(rref(J))==24
    names=[f'r{i}'for i in range(12)];monomials12=list(it.combinations_with_replacement(range(12),2));gens=[]
    for index,text in zip(dim['original350_row_indices'],dim['generators']):
        g=parse(text,names);target={tuple((i,j).count(k)for k in range(12)):c for col,c in J[index].items()for i,j in [monomials12[col]]};assert set(g)==set(target)
        ratio=g[next(iter(g))]/target[next(iter(g))];assert ratio and g==scale(target,ratio);gens.append(g)
    assert len(gens)==24
    lifted=[parse(text,names)for text in dim['lifted_polynomials']];sums=[{}for _ in lifted]
    for item in dim['ideal_membership_lifts']:sums[item['polynomial']]=add(sums[item['polynomial']],mul(gens[item['generator']],parse(item['coefficient'],names)))
    assert len(lifted)==62 and len(dim['ideal_membership_lifts'])==342 and sums==lifted
    leading=[max(g,key=lambda e:(sum(e),tuple(-a for a in reversed(e))))for g in lifted];assert leading==[tuple(e)for e in dim['leading_monomial_exponents']]
    supports=[sum(1<<i for i,a in enumerate(e)if a)for e in leading];counter=[]
    for variables in it.combinations(range(12),5):
        mask=sum(1<<i for i in variables);witness=next((i for i,edge in enumerate(supports)if not edge&mask),None);assert witness is not None;counter.append(witness)
    cover=next(v for v in it.combinations(range(12),6)if all(edge&sum(1<<i for i in v)for edge in supports))
    assert dim['leading_support_minimum_vertex_cover']==6 and dim['fixed_Poisson_source_dimension_upper']==6 and dim['whole_projective_PGL_image_dimension_upper']==23
    result=dict(task_id='G',all_passed=True,scope='Smooth cubic scroll rankzero plus irreducible reduced quintic only; no cone, plane, rank2 restriction, E coverage or new component claim.',
        whole30_source_reconstructed=True,linear_shape=[159,30],exact_linear_rank=18,complete_kernel_dimension=12,nonzero_rank_minor=str(det),all350_original_and_reduced_Jacobi_reconstructed=True,all126_original_and_reduced_actual_D_reconstructed=True,reduced_Jacobi_row_rank=24,
        all24_generators_original_Jacobi_constant_multiples=True,all62_polynomials_and342_original_ideal_lifts_verified=True,all62_grevlex_leading_monomials_verified=True,all792_size5_sets_fail_cover=True,size5_counterexample_leading_indices=counter,size6_cover=list(cover),fixed_source_dimension_upper=6,preserving_group_dimension=6,whole_projective_PGL_image_dimension_upper=23,
        dimension_logic='The62 literal leading monomials belong to the initial ideal of the full Jacobi ideal. Their supports have minimum vertex cover6, so dim<=12-6. The preserving group and scalar give free seven-dimensional fibers on PGL5 times nonzero source, independently of individual stabilizers.',
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in [Path(__file__),oldpath,sourcepath,dimpath]},reviewer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,fixed_source_dimension_upper=6,whole_projective_PGL_image_dimension_upper=23,seconds=result['seconds'])))


if __name__=='__main__':main()
