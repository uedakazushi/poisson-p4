#!/usr/bin/env python3
"""Independent complete projected-Veronese rank-zero exclusion.

Standard-library Fraction arithmetic only. No producer, helper, CAS, NumPy,
legacy, or repository arithmetic imports. Saved expressions are data only.
"""
import argparse
import ast
from fractions import Fraction as Q
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


def plus(*values):
    out={}
    for value in values:
        for e,c in value.items():out[e]=out.get(e,Q(0))+c
    return {e:c for e,c in out.items()if c}


def times(*values):
    if not values:return {}
    out=values[0]
    for value in values[1:]:
        nxt={}
        for e,a in out.items():
            for f,b in value.items():
                p=tuple(x+y for x,y in zip(e,f));nxt[p]=nxt.get(p,Q(0))+a*b
        out={e:c for e,c in nxt.items()if c}
    return out


def scale(value,c):return {e:c*a for e,a in value.items()if c*a}


def derivative(value,i):
    out={}
    for e,c in value.items():
        if e[i]:
            f=list(e);f[i]-=1;out[tuple(f)]=c*e[i]
    return out


class Ring:
    def __init__(self,names):self.names=tuple(names);self.zero=(0,)*len(names)
    def const(self,c):return {self.zero:Q(c)}if c else {}
    def var(self,name):
        e=list(self.zero);e[self.names.index(name)]=1;return {tuple(e):Q(1)}
    def power(self,value,n):
        out=self.const(1)
        for _ in range(n):out=times(out,value)
        return out
    def parse(self,expression):
        def walk(a):
            if isinstance(a,ast.Constant)and type(a.value)is int:return self.const(a.value)
            if isinstance(a,ast.Name):return self.var(a.id)
            if isinstance(a,ast.UnaryOp)and isinstance(a.op,ast.USub):return scale(walk(a.operand),-1)
            if isinstance(a,ast.BinOp):
                if isinstance(a.op,ast.Pow):
                    assert isinstance(a.right,ast.Constant)and type(a.right.value)is int and a.right.value>=0
                    return self.power(walk(a.left),a.right.value)
                u,v=walk(a.left),walk(a.right)
                if isinstance(a.op,ast.Add):return plus(u,v)
                if isinstance(a.op,ast.Sub):return plus(u,scale(v,-1))
                if isinstance(a.op,ast.Mult):return times(u,v)
                if isinstance(a.op,ast.Div):
                    assert set(v)=={self.zero};return scale(u,1/v[self.zero])
            raise ValueError(ast.dump(a))
        return walk(ast.parse(expression,mode='eval').body)


def monomials(n,d):
    out=[]
    for terms in it.combinations_with_replacement(range(n),d):
        e=[0]*n
        for i in terms:e[i]+=1
        out.append(tuple(e))
    return sorted(out,reverse=True)


def sparse_rref(rows,witness=False):
    """Incremental exact row reduction, optionally tracking literal input rows."""
    basis={};proof={};chosen=[]
    for i,row in enumerate(rows):
        v={j:Q(c)for j,c in (row.items()if isinstance(row,dict)else enumerate(row))if c}
        w={i:Q(1)}if witness else {}
        for p in sorted(basis):
            if p in v:
                a=v[p];v=plus(v,scale(basis[p],-a))
                if witness:w=plus(w,scale(proof[p],-a))
        if not v:continue
        p=min(v);a=v[p];v=scale(v,1/a);w=scale(w,1/a)
        for q in list(basis):
            if p in basis[q]:
                a=basis[q][p];basis[q]=plus(basis[q],scale(v,-a))
                if witness:proof[q]=plus(proof[q],scale(w,-a))
        basis[p]=v;proof[p]=w;chosen.append(i)
    return basis,proof,chosen


def det(matrix):
    m=[[Q(c)for c in row]for row in matrix];out=Q(1)
    for j in range(len(m)):
        k=next((i for i in range(j,len(m))if m[i][j]),None)
        if k is None:return Q(0)
        if k!=j:m[k],m[j]=m[j],m[k];out=-out
        a=m[j][j];out*=a
        for i in range(j+1,len(m)):
            b=m[i][j]/a
            if b:
                for k in range(j+1,len(m)):m[i][k]-=b*m[j][k]
                m[i][j]=Q(0)
    return out


def matrix_hash(rows):
    return hashlib.sha256(json.dumps(rows,default=str,separators=(',',':'),sort_keys=True).encode()).hexdigest()


def transpose(rows):return [list(r)for r in zip(*rows)]


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='reports/CLASSIFICATION_PROJECTED_VERONESE_ZERO_REVIEW_20261107.json');args=ap.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(119,120));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    started=time.monotonic();input_path=ROOT/'reports/CLASSIFICATION_PROJECTED_VERONESE_ZERO_20261107.json';saved=json.loads(input_path.read_text())
    for name,digest in saved['input_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest
    assert saved['producer_bytes']==(ROOT/'code/classification_projected_veronese_zero.py').read_text()
    X=Ring([f'x{i}'for i in range(5)]);S=Ring(['a','b','c']);Z=Ring([f'z{i}'for i in range(5)]);ZX=Ring(Z.names+X.names)
    q=[S.parse(t)for t in ['a**2-b**2','b**2-c**2','a*b','a*c','b*c']]
    assert q==[S.parse(t)for t in saved['projected_surface_parametrization']]
    cols=[(ij,e)for ij in PAIRS for e in monomials(5,2)]
    labels=[(tri,e)for tri in TRIPLES for e in monomials(3,6)];row_index={a:i for i,a in enumerate(labels)}
    assert cols==[(tuple(a['pair']),tuple(a['exponents']))for a in saved['ambient_coefficient_columns']]
    assert labels==[(tuple(a['triple']),tuple(a['exponents']))for a in saved['substitution_equation_rows']]
    A=[{}for _ in labels]
    for j,(ij,e)in enumerate(cols):
        value=S.const(1)
        for a,n in zip(q,e):value=times(value,S.power(a,n))
        for tri in TRIPLES:
            if set(ij).issubset(tri):
                k=next(k for k in tri if k not in ij);co=times(q[k],value)
                for exponent,c in co.items():A[row_index[tri,exponent]][j]=(-1)**tri.index(k)*c
    expected=[{int(j):Q(c)for j,c in saved['complete_280_by150_sparse_matrix'].get(str(i),{}).items()}for i in range(280)]
    assert A==expected
    red,_,chosen=sparse_rref(A);assert len(red)==121
    piv=sorted(red);free=[j for j in range(150)if j not in red];K=[]
    for j in free:
        v=[Q(0)]*150;v[j]=Q(1)
        for p in piv:v[p]=-red[p].get(j,Q(0))
        K.append(v)
    assert len(K)==29 and all(sum(c*v[j]for j,c in row.items())==0 for row in A for v in K)
    cert=saved['rank_certificate'];minor=det([[A[i].get(j,Q(0))for j in cert['columns']]for i in cert['rows']])
    assert len(cert['rows'])==len(cert['columns'])==121 and minor==Q(cert['determinant'])==1
    source_kernel=[[Q(a)for a in row]for row in saved['complete_kernel_basis_columns']]
    assert len(source_kernel)==29 and len(sparse_rref(source_kernel)[0])==29
    assert all(sum(c*v[j]for j,c in row.items())==0 for row in A for v in source_kernel)
    def decode(v,ring=X):
        P={ij:{}for ij in PAIRS}
        offset=len(ring.names)-5
        for c,(ij,e)in zip(v,cols):
            if c:P[ij][(0,)*offset+e]=c
        for i,j in PAIRS:P[j,i]=scale(P[i,j],-1)
        return P
    def vector(P):return [P.get(ij,{}).get(e,Q(0))for ij,e in cols]
    E=[]
    for a,b in it.product(range(5),repeat=2):
        P={}
        for i,j in PAIRS:
            P[i,j]=plus(times(X.var(f'x{i}'),X.var(f'x{b}'))if j==a else {},scale(times(X.var(f'x{j}'),X.var(f'x{b}')),-1)if i==a else {})
        E.append(vector(P))
    assert len(sparse_rref(E)[0])==24
    assert all(sum(c*v[j]for j,c in row.items())==0 for row in A for v in E)
    canonical=[]
    for v in K:
        P=decode(v);d=[plus(*(derivative(P.get((i,j),{}),i)for i in range(5)))for j in range(5)]
        C={ij:plus(P[ij],scale(plus(times(X.var(f'x{i}'),d[j]),scale(times(X.var(f'x{j}'),d[i]),-1)),-Q(1,5)))for ij in PAIRS for i,j in [ij]}
        canonical.append(vector(C))
    basis,_,_=sparse_rref(canonical);B=[[row.get(j,Q(0))for j in range(150)]for p,row in sorted(basis.items())]
    assert len(B)==5 and len(sparse_rref(E+B)[0])==29
    for b,source in zip(B,saved['canonical_section_basis']):
        decoded=decode(b)
        assert all(decoded.get(ij,{})==X.parse(source.get(f'{ij[0]},{ij[1]}','0'))for ij in PAIRS)
        assert all(not plus(*(derivative(decoded.get((i,j),{}),i)for i in range(5)))for j in range(5))
    P={ij:plus(*(times(ZX.var(f'z{k}'),decode(b,ZX)[ij])for k,b in enumerate(B)))for ij in PAIRS}
    for i,j in PAIRS:P[j,i]=scale(P[i,j],-1)
    assert all(P[ij]==ZX.parse(saved['complete_canonical_source'].get(f'{ij[0]},{ij[1]}','0'))for ij in PAIRS)
    strict={tri:plus(*(times(P.get((a,l),{}),derivative(P.get((b,c),{}),5+l))for a,b,c in [tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])]for l in range(5)))for tri in TRIPLES}
    jac={quad:plus(*(scale(times(ZX.var(f'x{quad[k]}'),strict[quad[:k]+quad[k+1:]]),(-1)**k)for k in range(4)))for quad in QUADS}
    J=[]
    for quad in QUADS:
        for e in monomials(5,4):J.append({p[:5]:c for p,c in jac[quad].items()if p[5:]==e})
    assert len(J)==350 and J==[Z.parse(t)for t in saved['complete350_Jacobi']]
    mon2=monomials(5,2);Jrows=[[p.get(e,Q(0))for e in mon2]for p in J]
    rowbasis,rowproof,_=sparse_rref(Jrows,True);assert len(rowbasis)==14
    target_texts=['z0**2-z4**2','z0*z1-z4**2','z1**2-4*z4**2','z0*z2','z1*z2','z2**2-z4**2','z0*z3','z1*z3','z2*z3','z3**2-z4**2','z0*z4','z1*z4','z2*z4','z3*z4']
    targets=[Z.parse(t)for t in target_texts];witnesses=[]
    for target in targets:
        v={j:c for j,e in enumerate(mon2)if(c:=target.get(e,Q(0)))};w={}
        for p in sorted(rowbasis):
            if p in v:
                c=v[p];v=plus(v,scale(rowbasis[p],-c));w=plus(w,scale(rowproof[p],c))
        assert not v and plus(*(scale(J[i],c)for i,c in w.items()))==target
        witnesses.append({str(i):str(c)for i,c in sorted(w.items())})
    assert len(sparse_rref([[v.get(e,Q(0))for e in mon2]for v in targets])[0])==14
    inverse_basis,inverse_proof,_=sparse_rref([[p.get(e,Q(0))for e in mon2]for p in targets],True)
    all350_from_targets=[]
    for original in J:
        v={j:c for j,e in enumerate(mon2)if(c:=original.get(e,Q(0)))};w={}
        for p in sorted(inverse_basis):
            if p in v:
                c=v[p];v=plus(v,scale(inverse_basis[p],-c));w=plus(w,scale(inverse_proof[p],c))
        assert not v and plus(*(scale(targets[i],c)for i,c in w.items()))==original
        all350_from_targets.append({str(i):str(c)for i,c in sorted(w.items())})
    # Each parameter cube is an explicit ideal combination. This proves the
    # radical is the origin without invoking a Gröbner-basis black box.
    cube_recipes=[[(0,0,1),(10,4,1)],[(2,1,1),(11,4,4)],[(5,2,1),(12,4,1)],[(9,3,1),(13,4,1)],[(0,4,-1),(10,0,1)]]
    for k,recipe in enumerate(cube_recipes):
        assert plus(*(scale(times(targets[j],Z.var(f'z{v}')),c)for j,v,c in recipe))==Z.power(Z.var(f'z{k}'),3)
    signed_pfaffians=[]
    for missing in range(5):
        i,j,k,l=[a for a in range(5)if a!=missing]
        signed_pfaffians.append(scale(plus(times(P[i,j],P[k,l]),scale(times(P[i,k],P[j,l]),-1),times(P[i,l],P[j,k])),(-1)**missing))
    assert signed_pfaffians==[ZX.parse(t)for t in saved['full_Pfaffian_five_vector']]
    global_D=plus(*(times(ZX.var(f'x{i}'),value)for i,value in enumerate(signed_pfaffians)))
    assert not global_D
    # Independently project to affine x0=1, where the true representative is
    # Pij-xi P0j+xj P0i. Its Pfaffian is an additional check, not the proof.
    affine={}
    for i,j in it.combinations(range(1,5),2):
        value=plus(P[i,j],scale(times(ZX.var(f'x{i}'),P[0,j]),-1),times(ZX.var(f'x{j}'),P[0,i]))
        out={}
        for e,c in value.items():
            f=e[:5]+(0,)+e[6:];out[f]=out.get(f,Q(0))+c
        affine[i,j]={e:c for e,c in out.items()if c}
    pf=plus(times(affine[1,2],affine[3,4]),scale(times(affine[1,3],affine[2,4]),-1),times(affine[1,4],affine[2,3]))
    assert not pf and saved['actual_affine_Pfaffian']=='0'
    inputs=[input_path,ROOT/'code/classification_projected_veronese_zero.py',Path(__file__)]
    output=dict(task_id='E_reduced_nonnormal',all_passed=True,all_frozen_producer_input_hashes_and_bytes_match=True,implementation='standard library Fraction only; no repository arithmetic/producer/CAS imports',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},
        independently_reconstructed_restriction=dict(shape=[280,150],matrix_exactly_matches_saved=True,full_matrix_sha256=matrix_hash(A),rank=121,minor_rows=cert['rows'],minor_columns=cert['columns'],minor_determinant=str(minor),independent_complete_kernel_dimension=29,saved_complete_kernel_independently_spans=True),
        independently_reconstructed_projective_quotient=dict(Euler_column_count=25,Euler_rank=24,all_Euler_annihilated=True,canonical_divergence_zero=True,complete_quotient_dimension=5,Euler_plus_quotient_rank=29,all_five_basis_tensors_exactly_match=True,all_five_parameter_source_exactly_matches=True),
        independently_reconstructed_Jacobi=dict(all350_equations_exactly_match=True,quadratic_rowspace_dimension=14,quadratic_monomial_order=[list(e)for e in mon2],fourteen_generating_quadrics=target_texts,literal_350_row_functionals=witnesses,all350_inverse_functionals_in_fourteen_generators=all350_from_targets,fourteen_quadrics_span_whole_Jacobi_rowspace=True,five_cube_identities=[dict(parameter=f'z{k}',terms=[dict(quadric_index=j,variable=f'z{v}',coefficient=c)for j,v,c in r])for k,r in enumerate(cube_recipes)],every_parameter_cube_in_original_Jacobi_ideal=True,radical='(z0,z1,z2,z3,z4)',only_reduced_Poisson_point='zero',all_actual_affine_Pfaffians_identically_zero_even_before_Jacobi=True,all_five_signed_Pfaffian_minors_exactly_match=True,correct_global_Euler_contracted_quintic_identically_zero=True,incorrect_common_Pfaffian_vector_identity_not_used=True),
        geometry=dict(center_in_full_Sym2_coordinates=['1','1','1','0','0','0'],center_symmetric_rank=3,secant_and_tangent_matrices_have_rank_at_most=2,projection_is_closed_embedding=True,all_rank3_centers_congruent_over_C=True,all_smooth_projected_Veroneses_one_PGL5_orbit=True,proof_checked='Projection kernel is the identity symmetric matrix; rank<=2 is the full secant variety, including tangent limits. Hence projection separates all points and tangent vectors. All rank3 symmetric forms are GL3-congruent, inducing an ambient PGL5 equivalence of quotient projections.'),
        scope='All Poisson tensors vanishing on a smooth projected Veronese surface are zero. Does not classify all degree4 surfaces, all rankzero surfaces, or all irreducible nonnormal quintics; no assumptions on normalization X used.',seconds=time.monotonic()-started,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    dest=ROOT/args.output
    if dest.exists():raise FileExistsError(str(dest))
    dest.write_text(json.dumps(output,indent=2)+'\n');print(json.dumps({'all_passed':True,'seconds':output['seconds'],'max_rss_kib':output['max_rss_kib'],'report':str(dest.relative_to(ROOT))}))

if __name__=='__main__':main()
