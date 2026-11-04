#!/usr/bin/env python3
"""Read-only, independent exact audit of the frozen full134 scalar circuits.

No producer is imported or run.  Parameter polynomials are independently
decoded and shifted.  Scalar circuit comparison takes place in a stronger
formal polynomial ring whose symbols are primitive parameter atoms and the
already audited E entries; this never expands the large parameter fractions.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
import gzip, hashlib, itertools, json, math, resource, signal, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ID = 'S61_SPARSE_GENERIC5_FULL134_NORMAL_CIRCUIT_20261104_v63_REVIEW'
AUTHOR = ID.removesuffix('_REVIEW')
ZERO = (0,) * 5

def read(p):
    raw = p.read_bytes()
    return json.loads(gzip.decompress(raw) if p.suffix == '.gz' else raw)

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def plus(a, b):
    c = dict(a)
    for m, x in b.items():
        c[m] = c.get(m, F(0)) + x
        if not c[m]: del c[m]
    return c

def times(a, b):
    c = {}
    for m, x in a.items():
        for n, y in b.items():
            k = tuple(sorted(m + n))
            c[k] = c.get(k, F(0)) + x*y
    return {k:x for k,x in c.items() if x}

def neg(a): return {k:-v for k,v in a.items()}
def const(x): return {():F(x)} if x else {}
def symbol(x): return {(x,):F(1)}

# Independent five-parameter polynomial arithmetic (different monomial law
# from the formal-symbol ring above).
def pp_add(a, b): return plus(a, b)
def pp_mul(a, b):
    c = {}
    for m, x in a.items():
        for n, y in b.items():
            k = tuple(u+v for u,v in zip(m,n))
            c[k] = c.get(k,F(0)) + x*y
    return {k:x for k,x in c.items() if x}

def native(p):
    return {tuple(t[1:]):F(p['content'])*int(t[0]) for t in p['terms']}

def translate(row, C):
    """Expand each original monomial by choosing center or new variable.

    This generic substitution, rather than a derivative formula, is used in
    both directions.  Coefficients remain exact five-variable QQ polynomials.
    """
    out = {}
    for mon, coef in row.items():
        partial = {():coef}
        for i in mon:
            nxt = {}
            for m, p in partial.items():
                nxt[m] = pp_add(nxt.get(m,{}), pp_mul(p,C[i]))
                k = tuple(sorted(m+(i,)))
                nxt[k] = pp_add(nxt.get(k,{}),p)
            partial = {m:p for m,p in nxt.items() if p}
        for m,p in partial.items(): out[m] = pp_add(out.get(m,{}),p)
    return {m:p for m,p in out.items() if p}

class FormalReader:
    def __init__(self,nodes,stops=None):
        self.nodes=nodes
        self.stops=stops or {}
        self.memo={}
    def __call__(self,i):
        if i in self.stops: return self.stops[i]
        if i in self.memo: return self.memo[i]
        op,*a=self.nodes[i]
        if op=='qq': p=const(F(a[0]))
        elif op=='parameter_poly_atom': p=symbol(('atom',a[0]))
        elif op.startswith('external_parameter_'):
            p=symbol((op.removeprefix('external_parameter_'),*a))
        elif op in ('z','n','q_input'): p=symbol((op,a[0]))
        elif op in ('sum','prod'):
            p={} if op=='sum' else const(1)
            for j in a[0]: p=plus(p,self(j)) if op=='sum' else times(p,self(j))
        elif op=='inv':
            q=self(a[0])
            assert q
            if set(q)=={()}: p=const(1/q[()])
            else:
                # A formal inverse symbol is keyed by its ENTIRE independently
                # normalized denominator.  We never invert an unknown.
                key=tuple((m,str(c)) for m,c in sorted(q.items()))
                p=symbol(('inverse',key))
        else: raise AssertionError(op)
        self.memo[i]=p
        return p

def main(out):
    started=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('60s audit cap')))
    signal.alarm(60)
    out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/AUTHOR
    parent=read(folder/'report.json')
    assert parent['all_passed']
    hashes={}
    for name,h in parent['input_sha256'].items():
        assert sha(ROOT/name)==h, name
        hashes[name]=h
    for name,h in parent['evidence_sha256'].items():
        assert sha(folder/name)==h,name
        hashes[str((folder/name).relative_to(ROOT))]=h
    hashes[str((folder/'report.json').relative_to(ROOT))]=sha(folder/'report.json')
    prior=ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_RANK3_CIRCUIT_20261104_v63_REVIEW/report.json'
    gate=read(prior)
    assert gate['all_passed'] and gate['rank_L18_plane3_normal42']
    # The prior review independently reconstructed the original physical134
    # source and audited every original1608 adapter operation.
    for path in (
        'reports/S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63/pulled134-input.json.gz',
        'reports/S61_SPARSE_GENERIC5_PLANE_RANK3_CIRCUIT_20261104_v63/normal42-adapter-SLP.json',
    ):
        assert sha(ROOT/path)==gate['input_sha256'][path]
    hashes[str(prior.relative_to(ROOT))]=sha(prior)
    d=read(folder/'full134-direct-circuit.json.gz')
    b=read(folder/'coefficient-block-circuit.json.gz')
    atoms=read(folder/'scalar-atoms-full.json.gz')
    oldatoms=read(folder/'scalar-atoms-direct.json.gz')
    assert oldatoms['primitive_integer_atoms']==atoms['primitive_integer_atoms'][:len(oldatoms['primitive_integer_atoms'])]
    source=ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_RANK3_CIRCUIT_20261104_v63'
    original=read(source/'plane-W-original.json.gz')
    slp=read(source/'normal42-adapter-SLP.json')
    ranks=read(source/'rank-certificate.json.gz')
    pulled=read(ROOT/'reports/S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63/pulled134-input.json.gz')
    assert atoms['variables']==oldatoms['variables']==pulled['parameter_variables']
    piv=d['pivot_q_rows']; normal=d['normal_q_rows42']
    assert piv==slp['pivot_q_rows'] and normal==[i for i in range(45) if i not in piv]
    assert len(piv)==3 and len(normal)==42
    assert d['unknown_order']==['z0','z1','z2']+[f'n{i}' for i in range(42)]
    assert d['q_to_u_indices_zero_based']==pulled['q_to_u_indices_zero_based']==slp['q_to_u_indices_zero_based']
    assert b['prefix_node_count']==len(d['nodes'])
    nodes=d['nodes']+b['appended_nodes']

    # Every saved scalar atom is literally a primitive integer polynomial.
    atomkeys={}
    for i,enc in enumerate(atoms['primitive_integer_atoms']):
        p={tuple(t[1:]):int(t[0]) for t in enc}
        assert len(p)==len(enc) and p and all(len(k)==5 and min(k)>=0 and v for k,v in p.items())
        assert [tuple(t[1:]) for t in enc]==sorted(p)
        assert math.gcd(*(abs(v) for v in p.values()))==1 and p[min(p)]>0
        key=tuple(sorted(p.items()))
        assert key not in atomkeys
        atomkeys[key]=i
    used_atoms=set()
    def atom(p):
        if not p:return {}
        if set(p)=={ZERO}: return const(p[ZERO])
        den=math.lcm(*(x.denominator for x in p.values()))
        integer={k:x.numerator*(den//x.denominator) for k,x in p.items()}
        gcd=math.gcd(*(abs(x) for x in integer.values()))
        if integer[min(integer)]<0:gcd=-gcd
        key=tuple((k,x//gcd) for k,x in sorted(integer.items()))
        aid=atomkeys[key];used_atoms.add(aid)
        return {(('atom',aid),):F(gcd,den)}

    # The adapter's rational expressions are independently compared as formal
    # polynomials in immutable external leaves and formal inverse factors.
    imp=FormalReader(d['nodes'])
    sv=[]
    for op,*a in slp['nodes']:
        if op=='qq': p=const(F(a[0]))
        elif op in ('original_W_poly','rank_pivot_rf','selected_row_factor_poly'):
            if op=='original_W_poly': nonzero=bool(original['W45x3'][a[0]][a[1]])
            elif op=='rank_pivot_rf': nonzero=bool(ranks['nonzero_pivot_factors'][a[0]]['numerator'])
            else: nonzero=bool(ranks['selected_row_factors'][a[0]]['terms'])
            p=symbol((op,*a)) if nonzero else {}
        elif op=='add':p=plus(sv[a[0]],sv[a[1]])
        elif op=='sub':p=plus(sv[a[0]],neg(sv[a[1]]))
        elif op=='mul':p=times(sv[a[0]],sv[a[1]])
        elif op=='inv':
            q=sv[a[0]];assert q
            p=const(1/q[()]) if set(q)=={()} else symbol(('inverse',tuple((m,str(c)) for m,c in sorted(q.items()))))
        elif op in ('z','n','center'):p=symbol((op,*a))
        else: raise AssertionError(op)
        sv.append(p)
    for i in range(45):
        for j in range(3):
            want=const(i==piv[j]) if i in piv else sv[slp['E45x3_nodes'][i][j]]
            assert imp(d['E45x3'][i][j])==want,(i,j)
    C=[native(p) for p in original['center45']]
    for c,ref in zip(C,d['center45']):assert imp(ref)==atom(c)

    # E is now treated as a scalar formal symbol, never as a syntactic proof
    # that its rational value is nonzero. Equal root IDs receive equal symbols.
    stops={}
    for row in d['E45x3']:
        for ref in row:
            if nodes[ref][0]!='qq':stops[ref]=symbol(('Eroot',ref))
    R=FormalReader(nodes,stops)
    E=[[R(ref) for ref in row] for row in d['E45x3']]
    for i,ref in enumerate(d['forward_q45']):
        want=atom(C[i])
        for j in range(3):want=plus(want,times(E[i][j],symbol(('z',j))))
        if i in normal:want=plus(want,symbol(('n',normal.index(i))))
        assert R(ref)==want
    for j,ref in enumerate(d['inverse_z3']):
        assert R(ref)==plus(symbol(('q_input',piv[j])),neg(atom(C[piv[j]])))
    for k,ref in enumerate(d['inverse_n42']):
        i=normal[k];want=plus(symbol(('q_input',i)),neg(atom(C[i])))
        for j in range(3):want=plus(want,neg(times(E[i][j],R(d['inverse_z3'][j]))))
        assert R(ref)==want

    # Both complete linear inverse products, with the 42 original coordinate
    # directions, are checked independently as sparse matrices.
    A=[]
    for i in range(45):
        row={j:E[i][j] for j in range(3) if E[i][j]}
        if i in normal:row[3+normal.index(i)]=const(1)
        A.append(row)
    B=[{p:const(1)} for p in piv]
    for i in normal:
        row={i:const(1)}
        for j in range(3):
            if E[i][j]:row[piv[j]]=neg(E[i][j])
        B.append(row)
    for X,Y in ((A,B),(B,A)):
        for i in range(45):
            row={}
            for k,x in X[i].items():
                for j,y in Y[k].items():row[j]=plus(row.get(j,{}),times(x,y))
            assert {j:p for j,p in row.items() if p}=={i:const(1)}

    # Bind every original term to the direct134 output without expanding any
    # forward coordinate: it is a polynomial in 45 formal q-forward symbols.
    Freader=FormalReader(nodes,{ref:symbol(('forward_q',i)) for i,ref in enumerate(d['forward_q45'])})
    rows=[];source_terms=0
    for idx,encoded in enumerate(pulled['equations']):
        row={};seen=set()
        for m,e,c in encoded:
            m=tuple(m);e=tuple(e);c=F(c)
            assert len(m)<=2 and tuple(sorted(m))==m and all(0<=i<45 for i in m)
            assert len(e)==5 and min(e)>=0 and c and (m,e) not in seen
            seen.add((m,e));row.setdefault(m,{})[e]=c;source_terms+=1
        assert row
        want={}
        for m,p in row.items():
            term=atom(p)
            for i in m:term=times(term,symbol(('forward_q',i)))
            want=plus(want,term)
        assert Freader(d['F134'][idx])==want,idx
        rows.append(row)
    assert len(rows)==len(d['F134'])==len(b['rows134'])==134

    # Independently expand centered quadratic tensors and then their 3+42
    # linear substitution. This computes every block in a single monomial map,
    # rather than copying the producer's four separate block formulas.
    count=0;shift_terms=0;block_entries=0
    for ridx,(row,saved) in enumerate(zip(rows,b['rows134'])):
        shifted=translate(row,C)
        reverse=translate(shifted,[{m:-v for m,v in c.items()} for c in C])
        assert reverse==row,ridx
        for m in [()]+[(i,) for i in range(45)]:
            assert reverse.get(m,{})==row.get(m,{})
            count+=1
        shift_terms+=sum(len(p) for m,p in shifted.items() if len(m)<2)
        expected={}
        for mon,p in shifted.items():
            work={():atom(p)}
            for i in mon:
                nxt={}
                for m,c in work.items():
                    for j,e in A[i].items():
                        k=tuple(sorted(m+(j,)))
                        nxt[k]=plus(nxt.get(k,{}),times(c,e))
                work=nxt
            for m,c in work.items():expected[m]=plus(expected.get(m,{}),c)
        actual={}
        assert len(saved['plane_z_degree_le2'])==10 and len(saved['linear_n42'])==42 and len(saved['cross_z_n126'])==126
        plane=[()]+[(i,) for i in range(3)]+list(itertools.combinations_with_replacement(range(3),2))
        for m,ref in zip(plane,saved['plane_z_degree_le2']):actual[m]=R(ref)
        for i,ref in enumerate(saved['linear_n42']):actual[(3+i,)]=R(ref)
        for i,ref in enumerate(saved['cross_z_n126']):actual[(i//42,3+i%42)]=R(ref)
        seen=set()
        for i,j,ref in saved['n2_sparse']:
            assert 0<=i<=j<42 and (i,j) not in seen
            seen.add((i,j));actual[(3+i,3+j)]=R(ref)
        assert {m:c for m,c in actual.items() if c}=={m:c for m,c in expected.items() if c},ridx
        # Include all missing n^2 zeros in the checked dense component count.
        block_entries+=10+42+126+903
    assert count==6164 and shift_terms==parent['shifted_scalar_terms']==16570

    # All 1081 basis monomials were covered by the tensor substitution rule:
    # constant,45 linear,1035 quadratic. Check its index-independent expansion
    # separately, using the exact saved E roots as formal scalar coefficients.
    basis_count=1+45
    for i in range(45):
        for j in range(i,45):
            direct={}
            for a,x in A[i].items():
                for c,y in A[j].items():
                    m=tuple(sorted((a,c)));direct[m]=plus(direct.get(m,{}),times(x,y))
            assembled={}
            for a in range(3):
                for c in range(a,3):
                    p=times(E[i][a],E[j][c])
                    if a!=c:p=plus(p,times(E[i][c],E[j][a]))
                    assembled[(a,c)]=p
            for v,other in ((i,j),(j,i)):
                if v in normal:
                    for a in range(3):
                        m=(a,3+normal.index(v));assembled[m]=plus(assembled.get(m,{}),E[other][a])
            if i in normal and j in normal:assembled[tuple(sorted((3+normal.index(i),3+normal.index(j))))]=const(1)
            assert {m:p for m,p in direct.items() if p}=={m:p for m,p in assembled.items() if p}
            basis_count+=1
    assert basis_count==1081

    # Degree and truth propagation are checked independently over the entire
    # saved DAG. Unknown sums remain UNDETERMINED; no field cancellation is
    # decided from syntax, and no unknown-dependent inverse is permitted.
    deg=[];truth=[];inv_count=0
    for i,(op,*a) in enumerate(nodes):
        refs=a[0] if op in ('sum','prod') else ([a[0]] if op=='inv' else [])
        assert all(isinstance(j,int) and 0<=j<i for j in refs)
        if op=='qq':di=0;ti=bool(F(a[0]))
        elif op=='parameter_poly_atom':assert 0<=a[0]<len(atomkeys);di=0;ti=True
        elif op.startswith('external_parameter_'):
            name=op.removeprefix('external_parameter_')
            if name=='original_W_poly':ti=bool(original['W45x3'][a[0]][a[1]])
            elif name=='rank_pivot_rf':ti=bool(ranks['nonzero_pivot_factors'][a[0]]['numerator'])
            else:assert name=='selected_row_factor_poly';ti=bool(ranks['selected_row_factors'][a[0]]['terms'])
            di=0
        elif op in ('z','n','q_input'):di=1;ti=True
        elif op=='prod':di=sum(deg[j] for j in refs);ti=False if any(truth[j] is False for j in refs) else (True if all(truth[j] is True for j in refs) else None)
        elif op=='sum':di=max(deg[j] for j in refs);ti=None
        else:assert op=='inv' and deg[a[0]]==0 and truth[a[0]] is not False;di=0;ti=True;inv_count+=1
        deg.append(di);truth.append(ti)
    assert inv_count==6 and all(deg[i]<=2 for i in d['F134'])
    def stats(refs):
        refs=list(refs); c=Counter(truth[i] for i in refs)
        return dict(entries=len(refs),proven_zero=c[False],proven_nonzero=c[True],undetermined=c[None],nonzero_lower=c[True],nonzero_upper=len(refs)-c[False])
    blocks=b['rows134']
    pattern=dict(n2=dict(entries=134*903,exact_nonzero=sum(len(t['n2_sparse']) for t in blocks),exact_zero=134*903-sum(len(t['n2_sparse']) for t in blocks)),
        linear_n=stats(i for t in blocks for i in t['linear_n42']),
        cross_z_n=stats(i for t in blocks for i in t['cross_z_n126']),
        plane=stats(i for t in blocks for i in t['plane_z_degree_le2']))
    for key,value in pattern.items():assert value==parent['patterns'][key]
    assert len(atomkeys)==4296 and len(nodes)==228650
    # No extra source polynomial is hidden in the atom table.
    common=read(ROOT/'reports/S61_SPARSE_GENERIC5_SMALL_TAIL_FIELD_20261104_v63/shared-denominator-kernel-plane.json.gz')['common_denominator_P']
    P=atom({tuple(t[1:]):F(t[0]) for t in common})
    sigma=[]
    for j,ref in enumerate(d['sigma_coefficients']):
        want=times(P,times(sv[slp['adjugate3x3_nodes'][2][j]],sv[slp['inverse_Gamma3_node']]))
        assert imp(ref)==want
        sigma.append(ref)
    sigreader=FormalReader(nodes,{d['sigma_node']:symbol(('sigma',))})
    s=symbol(('sigma',));s2=times(s,s)
    want=[plus(plus(const(1),s),times(const(-2),s2)),times(const(F(4,3)),plus(s2,neg(s))),times(const(F(1,3)),plus(times(const(2),s2),s))]
    assert [sigreader(i) for i in d['known5_idempotents_nodes']]==want
    for a,x in enumerate((F(0),F(-1,2),F(1))):
        assert [sum(c*x**len(m) for m,c in p.items()) for p in want]==[F(a==j) for j in range(3)]
    assert used_atoms==set(range(4296)),len(used_atoms)
    proof=dict(all134_original_pulled_polynomials_bound=True,original_scalar_terms=source_terms,
        all4296_primitive_integer_parameter_atoms_independently_bound=True,
        all134_direct_outputs_independently_reconstructed=True,
        prior1608_adapter_import_and135_E_entries_bound=True,
        all45_forward_and_inverse_maps_and4050_inverse_entries=True,
        exact_QQ_reverse_centering_identities=count,independently_shifted_scalar_terms=shift_terms,
        four_block_dense_entries=block_entries,all1081_basis_identities=True,
        all228650_nodes_unknown_degrees_and_all6_parameter_only_inverses=True,
        coefficient_patterns=pattern,known5_sigma_idempotents_only=True,
        no_unresolved_field_sum_promoted=True,no_new_plane_rank_cross_rank_or_exhaustion=True,
        no_author_builder_or_source_solver_run=True,no_new_point_prime_GB_J0=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in (Path(__file__),ROOT/f'research/notes/{ID}.md'):
        hashes[str(p.relative_to(ROOT))]=sha(p)
    report=dict(id=ID,all_passed=True,input_sha256=hashes,
        evidence_sha256={'independent-proof.json':sha(out/'independent-proof.json')},
        scope='Full134 equivalent affine45=3+42 problem over QQ(H,L,s0,q0,t0); no plane containment, full fiber or cross-rank claim.',
        seconds=time.monotonic()-started,max_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'],max_rss_KiB=report['max_rss_KiB'],proof=proof)))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    main(ap.parse_args().output)
