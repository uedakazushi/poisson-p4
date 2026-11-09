#!/usr/bin/env python3
"""Independent Fraction replay of the COMPLETE K4a one-parameter polynomial slice.

No producer, SymPy, CAS, or arithmetic helper is imported.  Existing independently
reviewed D/G/kernel certificates are reused, and the same certified 101 minor is
checked.  The sole new polynomial calculation is the literal full graph, all
350 original equations, and all degrees zero through ten.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import ast, hashlib, json, resource, signal, time

ROOT=Path(__file__).resolve().parents[1]
INPUT='reports/CLASSIFICATION_K4A_COMPLETE_TORUS_SLICE_20261108.json'
AMBIENT='reports/CLASSIFICATION_K4A_CANONICAL_FULL_AMBIENT_20261108.json'
NORMAL='reports/CLASSIFICATION_K4A_EXISTING_MATRIX_NORMAL_READBACK_20261108.json'
PRIOR='reports/CLASSIFICATION_K4A_CANONICAL_FULL_AMBIENT_SEPARATE_REVIEW_20261108.json'
PROOF='research/classification/K4A_EQUIVARIANT_FINITE_SLICE_PROOF.md'
OUT='reports/CLASSIFICATION_K4A_COMPLETE_TORUS_SLICE_SEPARATE_REVIEW_20261108.json'
JOB='research/classification/jobs/k4a_complete_torus_slice_separate_review.json'
PAIRS=list(combinations(range(5),2));TRIPLES=list(combinations(range(5),3));QUADS=list(combinations(range(5),4))
ZERO=(0,)*6


def add(*polys):
    out={}
    for p in polys:
        for e,c in p.items():out[e]=out.get(e,F())+c
    return{e:c for e,c in out.items()if c}
def scale(p,c):return{e:c*z for e,z in p.items()if c*z}
def mul(p,q):
    out={}
    for e,c in p.items():
        for f,d in q.items():
            ef=tuple(a+b for a,b in zip(e,f));out[ef]=out.get(ef,F())+c*d
    return{e:c for e,c in out.items()if c}
def der(p,i):
    out={}
    for e,c in p.items():
        if e[i]:f=list(e);f[i]-=1;out[tuple(f)]=c*e[i]
    return out
def coord(i):
    e=[0]*6;e[i]=1;return{tuple(e):F(1)}
def parse(text,names):
    def visit(n):
        if isinstance(n,ast.Constant):return{ZERO:F(n.value)}if n.value else{}
        if isinstance(n,ast.Name):return coord(names.index(n.id))
        if isinstance(n,ast.UnaryOp):return scale(visit(n.operand),-1 if isinstance(n.op,ast.USub)else 1)
        if isinstance(n,ast.BinOp):
            if isinstance(n.op,ast.Add):return add(visit(n.left),visit(n.right))
            if isinstance(n.op,ast.Sub):return add(visit(n.left),scale(visit(n.right),-1))
            if isinstance(n.op,ast.Mult):return mul(visit(n.left),visit(n.right))
            if isinstance(n.op,ast.Div):
                de=visit(n.right);assert set(de)=={ZERO};return scale(visit(n.left),1/de[ZERO])
            if isinstance(n.op,ast.Pow):
                assert isinstance(n.right,ast.Constant)and isinstance(n.right.value,int)and n.right.value>=0
                value={ZERO:F(1)}
                for _ in range(n.right.value):value=mul(value,visit(n.left))
                return value
        raise ValueError(ast.dump(n))
    return visit(ast.parse(text,mode='eval').body)
def jacobi(P):
    strict={}
    for triple in TRIPLES:
        cyclic=[triple,(triple[1],triple[2],triple[0]),(triple[2],triple[0],triple[1])]
        strict[triple]=add(*(mul(P.get((i,l),{}),der(P.get((j,k),{}),l))for i,j,k in cyclic for l in range(5)))
    return{quad:add(*(scale(mul(coord(quad[k]),strict[quad[:k]+quad[k+1:]]),(-1)**k)for k in range(4)))for quad in QUADS}
def determinant(A):
    A=[list(row)for row in A];n=len(A);value=F(1)
    for j in range(n):
        k=next((k for k in range(j,n)if A[k][j]),None)
        if k is None:return F()
        if k!=j:A[k],A[j]=A[j],A[k];value=-value
        pivot=A[j][j];value*=pivot;active=[k for k in range(j+1,n)if A[j][k]]
        for i in range(j+1,n):
            z=A[i][j]
            if z:
                ratio=z/pivot
                for k in active:A[i][k]-=ratio*A[j][k]
                A[i][j]=F()
    return value


def main():
    out=ROOT/OUT
    if out.exists():raise FileExistsError(out)
    started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,)*2)
    def timeout(signum,frame):raise TimeoutError('one literal complete-slice Fraction review budget90s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(90)
    job=json.loads((ROOT/JOB).read_text())
    for item in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/item['path']).read_bytes()).hexdigest()==item['sha256'],item['path']
    report=json.loads((ROOT/INPUT).read_text());old=json.loads((ROOT/AMBIENT).read_text());normal=json.loads((ROOT/NORMAL).read_text());prior=json.loads((ROOT/PRIOR).read_text())
    assert prior['all_independent_arithmetic_passed']and prior['full49dim_D_kernel_and_all51_tangent_columns_verified']
    columns=old['full_original_coefficient_order'];rows=old['full_original_equation_order'];names=old['coordinate_order'];assert names==['c','a','b','u','v']
    weights=[2,1,1,0,0];wc=[sum(weights[i]for i in col['pair'])-sum(a*b for a,b in zip(col['exponent'],weights))-1 for col in columns]
    wr=[sum(weights[i]for i in row['quad'])-sum(a*b for a,b in zip(row['exponent'],weights))-2 for row in rows]
    assert wc==report['full150_coordinate_weights']and wr==report['full350_equation_weights']and(min(wc),max(wc))==(-5,2)
    D=old['complete_sparse_matrices']['full350_D'];ds={tuple(map(int,key.split(','))):F(z)for key,z in D['entries'].items()}
    assert all(wr[i]==wc[j]for i,j in ds)
    cert=old['certificates']['full350_D'];R=cert['rows'];U=cert['columns'];assert R==report['selected101_equation_rows']and U==report['selected101_coordinate_columns']
    block=[[ds.get((i,j),F())for j in U]for i in R];det=determinant(block);assert det==F(cert['determinant'])==F(report['existing101_block_determinant'])and det
    K=old['full_D_kernel_sparse'];nindex=normal['kernel_basis_index'];assert nindex==23
    N=[F(K['entries'].get(str(i)+','+str(nindex),'0'))for i in range(150)];assert all(wc[i]==-1 for i,z in enumerate(N)if z)
    assert all(not sum((value*N[j]for(i0,j),value in ds.items()if i0==i),F())for i in range(350))
    M=old['complete_sparse_matrices']['Euler_GL5_scalar51'];ann={int(i):F(z)for i,z in normal['annihilator_of_all51_actual_columns'].items()}
    assert all(not sum((a*F(M['entries'].get(str(i)+','+str(j),'0'))for i,a in ann.items()),F())for j in range(51))
    annN=sum((z*N[i]for i,z in ann.items()),F());assert annN==F(8,87)==F(report['N_non_group_exact_functional_value'])
    point={tuple(map(int,key.split(','))):parse(z,names)for key,z in old['complete_actual_brackets'].items()}
    pvec=[point.get(tuple(col['pair']),{}).get(tuple(col['exponent'])+(0,),F())for col in columns]
    graph={(int(i),int(j)):F(z)for key,z in report['graph_coefficients_sparse'].items()for i,j in[key.split(',')]}
    assert report['graph_coefficient_array_shape']==[150,6]
    assert [graph.get((i,0),F())for i in range(150)]==pvec
    assert [graph.get((i,1),F())for i in range(150)]==N
    assert all(not z or wc[i]==0 for i,z in enumerate(pvec))
    for(i,j),z in graph.items():
        assert 0<=i<150 and 0<=j<=5 and z
        if j>=2:assert i in U and wc[i]==-j
    for step in report['all_four_fixed_recursion_steps']:
        degree=step['degree'];assert degree in[2,3,4,5]
        layer={str(i):str(graph[i,degree])for i in range(150)if(i,degree)in graph};assert layer==step['entire150_layer_sparse']
        assert all(graph.get((U[int(k)],degree),F())==F(z)for k,z in step['implicit_coefficients_sparse'].items())
    P={ij:{}for ij in PAIRS}
    for(i,j),z in graph.items():
        col=columns[i];P[tuple(col['pair'])][tuple(col['exponent'])+(j,)]=z
    for i,j in PAIRS:P[j,i]=scale(P[i,j],-1)
    allJ=jacobi(P);full={}
    for i,row in enumerate(rows):
        poly=allJ[tuple(row['quad'])]
        for degree in range(11):
            value=poly.get(tuple(row['exponent'])+(degree,),F())
            if value:full[i,degree]=value
        assert all(e[5]<=10 and sum(e[:5])==4 for e in poly)
    producerfull={tuple(map(int,key.split(','))):F(z)for key,z in report['all350_full_equation_coefficients_sparse'].items()}
    assert full==producerfull
    assert all(not value for(i,j),value in full.items()if i in R)
    assert not full and report['all350_exact_literal_graph_Jacobi_zero']
    result={'task_id':'E_reduced_nonnormal','all_independent_arithmetic_passed':True,'no_producer_SymPy_or_CAS_imported':True,'same_existing101_minor_reverified':str(det),'prior_full_D101_G48_kernel_certificate_reused':PRIOR,'all150_input350_output_weights_verified':True,'N_kernel_and_non_group_functional_verified':str(annN),'graph_all150_coefficients_and_four_weight_layers_verified':True,'all101_implicit_and350_full_equations_through_degree10_zero':True,'new_point_operator_rank_profile_or_jet':False,'formal_full_slice':'A1','raw_smooth_local_dimension':49,'projective_smooth_unique_component_dimension':24,'orbitclosure23_is_not_full_component':True,'component_name_or_new_count':'root recognition/distinctness obligation, not asserted','full_graph_sparse_verified':report['graph_coefficients_sparse'],'inputs':[{'path':item['path'],'sha256':item['sha256']}for item in job['exact_inputs_and_hashes']],'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'reviewer_source_bytes':Path(__file__).read_text(),'seconds':time.monotonic()-started,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'PASS':True,'all350_literal_polynomials_zero':True,'projective_dimension':24,'seconds':result['seconds'],'RSS_KiB':result['maximum_RSS_KiB']}))


if __name__=='__main__':main()
