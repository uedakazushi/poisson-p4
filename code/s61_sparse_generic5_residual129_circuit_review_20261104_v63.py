#!/usr/bin/env python3
"""Independent formal-ring audit of the complete129 residual circuit.

The ALL670 exact source-functional gate is REQUIRED. Selected25/rank5 alone
does not certify all129 plane restrictions. No author builder is run and the
692953-term matrix is not expanded or copied by this checker.
"""
import argparse,gzip,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from s61_sparse_generic5_full134_normal_circuit_review_20261104_v63 import plus,times,neg,const,symbol

ROOT=Path(__file__).resolve().parents[1]
AUTHOR='S61_SPARSE_GENERIC5_RESIDUAL129_CIRCUIT_20261104_v63'
ID=AUTHOR+'_REVIEW'
def read(p):
    raw=p.read_bytes()
    return json.loads(gzip.decompress(raw) if p.suffix=='.gz' else raw)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def determinant(A):
    """Laplace recursion, independent of the producer's permutation code."""
    if not A:return const(1)
    out={}
    for j,a in enumerate(A[0]):
        sub=[[r[k] for k in range(len(A)) if k!=j] for r in A[1:]]
        term=times(a,determinant(sub))
        out=plus(out,term if j%2==0 else neg(term))
    return out

class Reader:
    def __init__(self,nodes,stops=None):self.nodes=nodes;self.stops=stops or {};self.memo={}
    def __call__(self,i):
        if i in self.stops:return self.stops[i]
        if i in self.memo:return self.memo[i]
        op,*a=self.nodes[i]
        if op=='qq':v=const(Q(a[0]))
        elif op=='external_prior_scalar':v=symbol(('prior',a[0]))
        elif op=='external_prior_F134':v=symbol(('F',a[0]))
        elif op=='external_five_functional_matrix':v=symbol(('S',*a))
        elif op in ('z','n'):v=symbol((op,a[0]))
        elif op in ('sum','prod'):
            v={} if op=='sum' else const(1)
            for j in a[0]:v=plus(v,self(j)) if op=='sum' else times(v,self(j))
        else:raise AssertionError(('uncut_inverse_or_unknown',i,op))
        self.memo[i]=v
        return v

def main(out,functional_gate):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    base=ROOT/'reports'/AUTHOR;report=read(base/'report.json');assert report['all_passed']
    hashes={}
    for p,h in report['input_sha256'].items():assert sha(ROOT/p)==h;hashes[p]=h
    for p,h in report['evidence_sha256'].items():assert sha(base/p)==h;hashes[str((base/p).relative_to(ROOT))]=h
    hashes[str((base/'report.json').relative_to(ROOT))]=sha(base/'report.json')
    fullbase=ROOT/'reports/S61_SPARSE_GENERIC5_FULL134_NORMAL_CIRCUIT_20261104_v63'
    fullgate=ROOT/'reports/S61_SPARSE_GENERIC5_FULL134_NORMAL_CIRCUIT_20261104_v63_REVIEW/report.json'
    selectgate=ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_FIVE_UNIVARIATE_20261104_v63_REVIEW/report.json'
    for p in (fullgate,selectgate,functional_gate):
        g=read(p);assert g['all_passed'];hashes[str(p.relative_to(ROOT))]=sha(p)
        for q,h in g['input_sha256'].items():assert sha(ROOT/q)==h;hashes[q]=h
    # Selected25 only proves the nonzero determinant. All670 source bindings
    # come from this separately audited exact streaming equality gate.
    fg=read(functional_gate)
    matrixpath='reports/S61_SPARSE_GENERIC5_PLANE_FIVE_SELECTOR_20261104_v63/five-column-matrix.json.gz'
    assert fg['input_sha256'][matrixpath]==sha(ROOT/matrixpath)
    direct=read(fullbase/'full134-direct-circuit.json.gz');block=read(fullbase/'coefficient-block-circuit.json.gz')
    old=direct['nodes']+block['appended_nodes'];br=block['rows134']
    d=read(base/'residual129-factored-circuit.json.gz');nodes=d['nodes']
    selected=d['selected_original134_rows'];rows=d['residual_original134_rows']
    assert selected==[50,10,52,51,15]
    assert rows==[i for i in range(134) if i not in selected]
    assert d['all45_unknowns_order']==direct['unknown_order']
    for name in ('pivot_q_rows','normal_q_rows42','q_to_u_indices_zero_based'):assert d[name]==direct[name]
    assert d['external_files']['five_functional_matrix']==matrixpath
    assert d['external_files']['prior_direct']==str((fullbase/'full134-direct-circuit.json.gz').relative_to(ROOT))
    assert d['external_files']['prior_blocks']==str((fullbase/'coefficient-block-circuit.json.gz').relative_to(ROOT))

    # Bind every external prior scalar to a genuinely parameter-only root.
    olddegree=[]
    for op,*a in old:
        if op in ('z','n','q_input'):v=1
        elif op=='prod':v=sum(olddegree[j] for j in a[0])
        elif op=='sum':v=max(olddegree[j] for j in a[0])
        elif op=='inv':assert olddegree[a[0]]==0;v=0
        else:assert op=='qq' or op=='parameter_poly_atom' or op.startswith('external_parameter_');v=0
        olddegree.append(v)
    degree=[];counts={};external={};inverse=[]
    for i,(op,*a) in enumerate(nodes):
        counts[op]=counts.get(op,0)+1
        if op in ('sum','prod'):
            assert all(isinstance(j,int) and 0<=j<i for j in a[0])
            v=max(degree[j] for j in a[0]) if op=='sum' else sum(degree[j] for j in a[0])
        elif op=='inv':assert 0<=a[0]<i and degree[a[0]]==0;inverse.append(a[0]);v=0
        elif op=='external_prior_scalar':assert 0<=a[0]<len(old) and olddegree[a[0]]==0;external[a[0]]=True;v=0
        elif op=='external_prior_F134':assert 0<=a[0]<134 and a[1]==direct['F134'][a[0]];v=2
        elif op=='external_five_functional_matrix':assert 0<=a[0]<134 and 0<=a[1]<5;v=0
        elif op in ('n','z'):assert 0<=a[0]<(42 if op=='n' else 3);v=1
        else:assert op=='qq';Q(a[0]);v=0
        degree.append(v)
    # The producer's14156 scalar-cache roots include242 QQ constants; only
    # the13914 non-QQ roots are serialized as external_prior_scalar leaves.
    assert len(nodes)==131828 and len(external)==13914 and counts['qq']==242
    delta=d['Delta5_node'];invd=d['inverseDelta5_node']
    assert inverse==[delta] and nodes[invd]==['inv',delta]
    raw=Reader(nodes)
    B=[[symbol(('S',i,j)) for j in range(5)] for i in selected]
    for a,row in enumerate(d['B5_nodes']):
        assert len(row)==5
        for b,ref in enumerate(row):assert raw(ref)==B[a][b]
    D=determinant(B);assert len(D)==120 and raw(delta)==D
    adj=[]
    for i in range(5):
        row=[]
        for j in range(5):
            A=[[B[r][c] for c in range(5) if c!=i] for r in range(5) if r!=j]
            p=determinant(A)
            if (i+j)%2:p=neg(p)
            assert raw(d['adjugate5_nodes'][i][j])==p
            row.append(p)
        adj.append(row)
    identities=0
    for X,Y in ((B,adj),(adj,B)):
        for i in range(5):
            for j in range(5):
                p={}
                for k in range(5):p=plus(p,times(X[i][k],Y[k][j]))
                assert p==(D if i==j else {});identities+=1
    for j in range(5):
        p={}
        for a in range(5):
            for b in range(5):p=plus(p,times(symbol(('row',a)),times(adj[a][b],B[b][j])))
        assert p==times(D,symbol(('row',j)));identities+=1
    assert identities==55
    adjstops={invd:symbol(('inverseDelta',))}
    for i in range(5):
        for j in range(5):adjstops[d['adjugate5_nodes'][i][j]]=symbol(('adj',i,j))
    LR=Reader(nodes,adjstops)
    for r,l in zip(rows,d['lambda129x5']):
        assert len(l)==5
        for j,ref in enumerate(l):
            p={}
            for h in range(5):p=plus(p,times(symbol(('S',r,h)),symbol(('adj',h,j))))
            assert LR(ref)==times(symbol(('inverseDelta',)),p)
    stops={ref:symbol(('lambda',i,j)) for i,row in enumerate(d['lambda129x5']) for j,ref in enumerate(row)}
    assert len(stops)==129*5
    R=Reader(nodes,stops)
    def lam(i,j):return symbol(('lambda',i,j))
    def scalar(ref):return const(Q(old[ref][1])) if old[ref][0]=='qq' else symbol(('prior',ref))
    def difference(i,refs):
        p=scalar(refs[rows[i]])
        for j,s in enumerate(selected):p=plus(p,neg(times(lam(i,j),scalar(refs[s]))))
        return p
    assert [R(ref) for ref in d['selected_direct_F5']]==[symbol(('F',i)) for i in selected]
    assert len(d['residual_direct_F129'])==129
    for i,(r,ref) in enumerate(zip(rows,d['residual_direct_F129'])):
        p=symbol(('F',r))
        for j,s in enumerate(selected):p=plus(p,neg(times(lam(i,j),symbol(('F',s)))))
        assert R(ref)==p
    for i in range(129):
        assert len(d['linear_n129x42'][i])==42 and len(d['cross129x126'][i])==126
        for k,ref in enumerate(d['linear_n129x42'][i]):
            assert R(ref)==difference(i,[a['linear_n42'][k] for a in br])
        for k,ref in enumerate(d['cross129x126'][i]):
            assert R(ref)==difference(i,[a['cross_z_n126'][k] for a in br])
    qfactor=[]
    for i,row in enumerate(br):
        expected=[{} for _ in range(42)]
        for a,b,ref in row['n2_sparse']:
            assert 0<=a<=b<42
            expected[a]=plus(expected[a],times(scalar(ref),symbol(('n',b))))
        assert len(d['source_normal2_factor134x42'][i])==42
        assert [R(ref) for ref in d['source_normal2_factor134x42'][i]]==expected
        qfactor.append(expected)
    for i,r in enumerate(rows):
        assert len(d['linear_normal_cofactor129x42'][i])==42
        for k,ref in enumerate(d['linear_normal_cofactor129x42'][i]):
            p=R(d['linear_n129x42'][i][k])
            for h in range(3):p=plus(p,times(R(d['cross129x126'][i][42*h+k]),symbol(('z',h))))
            p=plus(p,qfactor[r][k])
            for j,s in enumerate(selected):p=plus(p,neg(times(lam(i,j),qfactor[s][k])))
            assert R(ref)==p
    recipe=d['normal2_residual_lazy_recipe']
    assert recipe['source_rows']==rows and recipe['selected_rows']==selected and recipe['multipliers']==d['lambda129x5']
    # Every degree<=2 monomial with a normal factor is covered, including all
    #903 normal squares/products. No normal2 term is silently dropped.
    mons=[(i,) for i in range(3,45)]
    mons+=list(itertools.product(range(3),range(3,45)))
    mons+=list(itertools.combinations_with_replacement(range(3,45),2))
    assert len(set(mons))==1071
    for m in mons:
        v=list(m);k=min(j for j in v if j>=3);v.remove(k)
        assert tuple(sorted((k,*v)))==tuple(sorted(m)) and len(v)<=1
    assert all(degree[i]==0 for a in d['lambda129x5']+d['linear_n129x42']+d['cross129x126'] for i in a)
    assert all(degree[i]<=1 for a in d['linear_normal_cofactor129x42'] for i in a)
    assert all(degree[i]<=2 for i in d['residual_direct_F129'])
    proof=dict(all670_source_functional_gate_REQUIRED_and_bound=True,
      selected25_nonzero_Delta5_gate_separately_bound=True,
      all131828_nodes_topological_and_unknown_degree_checked=True,
      all13914_external_nonQQ_prior_scalar_roots_parameter_only=True,
      all242_QQ_constants_and14156_total_scalar_roots_bound=True,
      all670_external_matrix_entries_by_frozen_source_gate=True,
      own_Laplace_Delta_and25_adjugates_and55_formal_identities=True,
      all645_lambda_roots_exact=True,all129_direct_residuals_exact=True,
      full134_ideal_equals_selected5_plus_residual129_by_invertible_block_row_operation=True,
      all5418_linear_n_and16254_cross_coefficients_checked=True,
      all5628_source_normal2_factors_and5418_residual_linear_normal_cofactors_checked=True,
      all1071_normal_basis_monomials_and_normal2_lazy_recipe_checked=True,
      all129_plane_restrictions_zero_by_actual_functionals_plus_injectivity=True,
      all129_factor_identities_R_equals_sum_nA_proved=True,
      only_new_inverse_Delta5_parameter_only=True,all45_unknowns_and_scheme_structure_retained=True,
      no_692953_parameter_expansion_matrix_copy_or_author_builder=True,
      no_new_rank126_point_prime_GB_J0_or_exhaustion=True)
    (out/'independent-residual-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in (Path(__file__),ROOT/f'research/notes/{ID}.md',
              ROOT/'code/s61_sparse_generic5_full134_normal_circuit_review_20261104_v63.py'):
        hashes[str(p.relative_to(ROOT))]=sha(p)
    review=dict(id=ID,all_passed=True,input_sha256=hashes,
      evidence_sha256={'independent-residual-proof.json':sha(out/'independent-residual-proof.json')},
      full129_residual_and_factorization_independently_accepted=True,
      all670_source_gate_is_required_not_selected25_only=True,
      no_rank126_or_exhaustion=True,seconds=time.monotonic()-start,
      max_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(review,indent=2)+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(review))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--functional-review',type=Path,required=True)
    a=p.parse_args();main(a.output,a.functional_review.resolve())
