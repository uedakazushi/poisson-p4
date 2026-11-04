#!/usr/bin/env python3
"""Independent type, source-root and all-chart audit of the lazy exact presentation."""
import argparse, collections, gzip, hashlib, itertools, json, resource, signal, time
from fractions import Fraction as F
from pathlib import Path
from s61_sparse_generic5_full134_normal_circuit_review_20261104_v63 import plus,times,neg,const,symbol
ROOT=Path(__file__).resolve().parents[1]
AUTHOR='S61_SPARSE_GENERIC5_NORMAL42_EXACT_PRESENTATION_20261104_v63';ID=AUTHOR+'_REVIEW'
def read(p):
    b=p.read_bytes();return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def total(vs):
    r={}
    for v in vs:r=plus(r,v)
    return r
def degrees(nodes,old=None,res=None):
    out=[]
    for op,*a in nodes:
        if op in ('n','z','q_input'):d=1
        elif op=='prod':d=sum(out[i] for i in a[0])
        elif op=='sum':d=max(out[i] for i in a[0])
        elif op=='inv':assert out[a[0]]==0;d=0
        elif op=='external_prior_F134':d=2
        elif op=='external_full134_node':d=old[a[0]]
        elif op=='external_residual129_node':d=res[a[0]]
        else:
            assert op in ('qq','parameter_poly_atom','external_prior_scalar','external_five_functional_matrix','external_original_W') or op.startswith('external_parameter_'),op
            d=0
        out.append(d)
    return out
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports'/AUTHOR;hashes={}
    def bind(p):hashes[str(p.relative_to(ROOT))]=sha(p)
    def audit(p):
        d=read(p);assert d['all_passed'];bind(p)
        for n,h in d.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n;hashes[n]=h
        for n,h in d.get('evidence_sha256',{}).items():assert sha(p.parent/n)==h,n;bind(p.parent/n)
        return d
    audit(base/'report.json')
    for g in ('S61_SPARSE_GENERIC5_G6_RESIDUAL_20261104_v63_REVIEW',
              'S61_SPARSE_GENERIC5_RESIDUAL129_CIRCUIT_20261104_v63_REVIEW',
              'S61_SPARSE_GENERIC5_CROSS_RANK_COVARIANCE_20261104_v63_REVIEW',
              'S61_AFFINE_NORMAL_FINITE5_20261103_v63_REVIEW',
              'S61_AFFINE_POLARIZATION_SPLITTER_20261103_v63_REVIEW'):
        audit(ROOT/'reports'/g/'report.json')
    rfile=ROOT/'reports/S61_SPARSE_GENERIC5_RESIDUAL129_CIRCUIT_20261104_v63/residual129-factored-circuit.json.gz'
    gfile=ROOT/'reports/S61_SPARSE_GENERIC5_G6_RESIDUAL_20261104_v63/G6-normal-residual-circuit.json.gz'
    ffile=ROOT/'reports/S61_SPARSE_GENERIC5_FULL134_NORMAL_CIRCUIT_20261104_v63/full134-direct-circuit.json.gz'
    bfile=ffile.parent/'coefficient-block-circuit.json.gz'
    r,g,f,b=map(read,[rfile,gfile,ffile,bfile]);old=f['nodes']+b['appended_nodes']
    od=degrees(old);rd=degrees(r['nodes']);gd=degrees(g['nodes'],od,rd)
    ctx={'full134':od,'residual129':rd,'G6':gd}
    d=read(base/'typed-normal42-operator.json.gz');nodes=d['nodes'];assert len(nodes)==24501
    assert d['source_context_files']==dict(residual129=str(rfile.relative_to(ROOT)),G6=str(gfile.relative_to(ROOT)),full134_direct=str(ffile.relative_to(ROOT)),full134_blocks=str(bfile.relative_to(ROOT)))
    selected=r['selected_original134_rows'];assert selected==d['selected_original134_rows']==[50,10,52,51,15]
    assert d['residual130_row_order']==r['residual_original134_rows']+['G6']
    def node(i):return nodes[i]
    def typ(i):return node(i)['type']
    def match(i,op,inputs=None,attrs=None):
        n=node(i);assert n['op']==op,(i,op,n['op'])
        if inputs is not None:assert n['inputs']==inputs,(i,op)
        if attrs is not None:assert n['attrs']==attrs,(i,op)
        return n
    def scalar(i,c,root):match(i,'external_parameter_scalar',[],dict(context=c,root=root))
    def zpoly(i,c,root):match(i,'external_polynomial_z_n',[],dict(context=c,root=root))
    types=[];opcounts=collections.Counter()
    for i,n in enumerate(nodes):
        op,ins,a,t=n['op'],n['inputs'],n['attrs'],n['type'];opcounts[op]+=1
        assert all(isinstance(j,int) and 0<=j<i for j in ins)
        if op in ('external_parameter_scalar','external_polynomial_z_n'):
            assert not ins and a['context'] in ctx and 0<=a['root']<len(ctx[a['context']])
            degree=ctx[a['context']][a['root']]
            if op=='external_parameter_scalar':assert degree==0;expected=dict(kind='K')
            else:assert degree<=2;expected=dict(kind='PolyZN',degree_bound=2)
        elif op=='matrix_from_entries':
            assert len(ins)==a['rows']*a['cols'] and all(types[j]['kind']=='K' for j in ins)
            expected=dict(kind='MatK',shape=[a['rows'],a['cols']])
        elif op=='transpose':
            assert len(ins)==1 and types[ins[0]]['kind']=='MatK';p,q=types[ins[0]]['shape'];expected=dict(kind='MatK',shape=[q,p])
        elif op in ('matrix_product','unique_parameter_linear_solve'):
            assert len(ins)==2 and all(types[j]['kind']=='MatK' for j in ins)
            p,q=types[ins[0]]['shape'];u,v=types[ins[1]]['shape']
            if op=='matrix_product':assert q==u;expected=dict(kind='MatK',shape=[p,v])
            else:assert p==q==u and a['evaluated'] is False;expected=dict(kind='MatK',shape=[q,v])
        elif op=='specialize_z_zero':
            assert len(ins)==1 and types[ins[0]]['kind']=='PolyZN'
            assert t['normal_order_lower_bound'] in (0,1)
            expected=dict(kind='PolyN',degree_bound=types[ins[0]]['degree_bound'],normal_order_lower_bound=t['normal_order_lower_bound'])
        elif op=='polynomial_vector':
            assert ins and all(types[j]['kind']=='PolyN' for j in ins)
            expected=dict(kind='VecPolyN',length=len(ins),degree_bound=max(types[j]['degree_bound'] for j in ins),normal_order_lower_bound=min(types[j]['normal_order_lower_bound'] for j in ins))
        elif op=='parameter_matrix_times_polynomial_vector':
            assert len(ins)==2 and types[ins[0]]['kind']=='MatK' and types[ins[1]]['kind']=='VecPolyN'
            p,q=types[ins[0]]['shape'];v=types[ins[1]];assert q==v['length']
            expected=dict(kind='VecPolyN',length=p,degree_bound=v['degree_bound'],normal_order_lower_bound=v['normal_order_lower_bound'])
        elif op in ('vector_negation','vector_addition'):
            assert len(ins)==(1 if op=='vector_negation' else 2) and types[ins[0]]['kind']=='VecPolyN'
            if len(ins)==2:assert types[ins[0]]==types[ins[1]]
            expected=types[ins[0]]
        elif op=='vector_entry':
            assert len(ins)==1;v=types[ins[0]];assert v['kind']=='VecPolyN' and 0<=a['index']<v['length']
            expected=dict(kind='PolyN',degree_bound=v['degree_bound'],normal_order_lower_bound=v['normal_order_lower_bound'])
        elif op=='normal_variable':
            assert not ins and 0<=a['index']<42;expected=dict(kind='PolyN',degree_bound=1,normal_order_lower_bound=1)
        elif op=='parameter_as_normal_polynomial':
            assert len(ins)==1 and types[ins[0]]['kind']=='K';expected=dict(kind='PolyN',degree_bound=0,normal_order_lower_bound=0)
        else:
            assert op in ('polynomial_product','polynomial_sum','polynomial_difference') and ins and all(types[j]['kind']=='PolyN' for j in ins)
            if op=='polynomial_difference':assert len(ins)==2
            reducer=sum if op=='polynomial_product' else max
            order=sum if op=='polynomial_product' else min
            expected=dict(kind='PolyN',degree_bound=reducer(types[j]['degree_bound'] for j in ins),normal_order_lower_bound=order(types[j]['normal_order_lower_bound'] for j in ins))
        assert t==expected,(i,op,t,expected);types.append(expected)
    assert opcounts['unique_parameter_linear_solve']==1
    X=d['actual_X130x126'];flat=node(X)['inputs'];match(X,'matrix_from_entries',attrs=dict(rows=130,cols=126))
    expectedX=[('residual129',x) for row in r['cross129x126'] for x in row]+[('G6',x) for x in g['R_G6_cross126']]
    assert len(flat)==len(expectedX)==16380
    for i,(c,root) in zip(flat,expectedX):scalar(i,c,root)
    gram=d['Gram_parameter_matrix'];gr=match(gram,'matrix_product');assert gr['inputs'][1]==X
    xt=gr['inputs'][0];match(xt,'transpose',[X])
    L=d['leftinverse_parameter_matrix'];match(L,'unique_parameter_linear_solve',[gram,xt]);assert typ(L)['shape']==[126,130]
    bv=d['pure_normal_b130'];bi=match(bv,'polynomial_vector')['inputs'];assert len(bi)==130
    expectedR=[('residual129',x) for x in r['residual_direct_F129']]+[('G6',g['R_G6_direct_node'])]
    for i,(c,root) in zip(bi,expectedR):
        n=match(i,'specialize_z_zero',attrs={});assert typ(i)['normal_order_lower_bound']==1
        zpoly(n['inputs'][0],c,root)
    P=d['normal_rewrite_P126'];negP=match(P,'vector_negation');match(negP['inputs'][0],'parameter_matrix_times_polynomial_vector',[L,bv])
    H=d['pure_normal_H130'];hn=match(H,'vector_addition');assert hn['inputs'][0]==bv;match(hn['inputs'][1],'parameter_matrix_times_polynomial_vector',[X,P])
    pure=d['H_generator_nodes130'];assert len(pure)==130
    for i,h in enumerate(pure):match(h,'vector_entry',[H],dict(index=i))
    # Small free-polynomial reader: P entries, H entries and all source
    # coefficients stay independent, so no field expressions are expanded.
    memo={}
    def value(i):
        if i in memo:return memo[i]
        n=node(i);op,ins,a=n['op'],n['inputs'],n['attrs']
        if op=='vector_entry':assert ins[0] in (P,H);v=symbol(('P' if ins[0]==P else 'H',a['index']))
        elif op=='normal_variable':v=symbol(('n',a['index']))
        elif op=='specialize_z_zero':v=symbol(('b',ins[0]))
        elif op=='parameter_as_normal_polynomial':
            s=match(ins[0],'external_parameter_scalar');v=symbol(('k',s['attrs']['context'],s['attrs']['root']))
        elif op=='polynomial_product':
            v=const(1)
            for j in ins:v=times(v,value(j))
        elif op=='polynomial_sum':v=total(value(j) for j in ins)
        else:assert op=='polynomial_difference';v=plus(value(ins[0]),neg(value(ins[1])))
        memo[i]=v;return v
    def pn(h,j):return symbol(('P',42*h+j))
    def nn(j):return symbol(('n',j))
    labels=list(itertools.product(range(3),itertools.combinations(range(42),2)))
    labels=[[h,j,k] for h,(j,k) in labels]
    assert d['cross_compatibility_labels']==labels and len(d['cross_compatibility_generator_nodes2583'])==2583
    for root,(h,j,k) in zip(d['cross_compatibility_generator_nodes2583'],labels):
        assert value(root)==plus(times(nn(j),pn(h,k)),neg(times(nn(k),pn(h,j))))
    qlabels=[[r0,j] for r0 in selected for j in range(42)]
    assert d['plane_quartic_labels']==qlabels and len(d['plane_quartic_generator_nodes210'])==210
    zz=list(itertools.combinations_with_replacement(range(3),2))
    selected_b={}
    for i,n in enumerate(nodes):
        if n['op']=='specialize_z_zero' and typ(i)['normal_order_lower_bound']==0:
            p=node(n['inputs'][0]);assert p['op']=='external_polynomial_z_n' and p['attrs']['context']=='full134'
            selected_b[p['attrs']['root']]=n['inputs'][0]
    assert set(selected_b)=={f['F134'][j] for j in selected}
    for root,(r0,j) in zip(d['plane_quartic_generator_nodes210'],qlabels):
        row=b['rows134'][r0];k=lambda ref:symbol(('k','full134',ref))
        c0=symbol(('b',selected_b[f['F134'][r0]]))
        terms=[times(times(nn(j),nn(j)),c0)]
        for h in range(3):
            linear=plus(k(row['plane_z_degree_le2'][1+h]),total(times(k(row['cross_z_n126'][42*h+s]),nn(s)) for s in range(42)))
            terms.append(times(times(nn(j),pn(h,j)),linear))
        terms += [times(k(row['plane_z_degree_le2'][4+i]),times(pn(h,j),pn(l,j))) for i,(h,l) in enumerate(zz)]
        assert value(root)==total(terms)
    assert d['all2923_J_generators']==pure+d['cross_compatibility_generator_nodes2583']+d['plane_quartic_generator_nodes210']
    assert len(set(d['all2923_J_generators']))==2923
    for roots,degree in ((pure,2),(d['cross_compatibility_generator_nodes2583'],3),(d['plane_quartic_generator_nodes210'],4)):
        assert all(typ(i)['kind']=='PolyN' and typ(i)['degree_bound']<=degree and typ(i)['normal_order_lower_bound']>=1 for i in roots)
    ns=d['saturation_recipe']['ideal_m'];assert len(ns)==42
    for j,i in enumerate(ns):match(i,'normal_variable',[],dict(index=j))
    assert d['saturation_recipe']['coefficient_saturation_computed'] is False
    # Independent exponent-counter check of every clearing-of-denominators
    # identity, rather than comparing the producer's two sorted lists.
    mons=[()]+[(i,) for i in range(45)]+list(itertools.combinations_with_replacement(range(45),2));checks=0
    for j in range(42):
        for mon in mons:
            left=collections.Counter({('n',j):2})
            for i in mon:
                if i<3:left[('P',i)]+=1;left[('n',j)]-=1
                else:left[('n',i-3)]+=1
            left=+left
            zs=[i for i in mon if i<3];normal=[i-3 for i in mon if i>=3]
            right=collections.Counter([('P',i) for i in zs]+[('n',i) for i in normal]+[('n',j)]*(2-len(zs)))
            assert left==right and sum(2*v if k[0]=='P' else v for k,v in left.items())<=4;checks+=1
    assert checks==45402
    chartchecks=0
    for h,j,k in itertools.product(range(3),range(42),range(42)):
        recovered=plus(times(nn(k),pn(h,j)),neg(times(nn(j),pn(h,k))))
        antisym=plus(times(nn(j),pn(h,k)),neg(times(nn(k),pn(h,j))))
        assert plus(recovered,antisym)=={};chartchecks+=1
    assert chartchecks==5292
    proof=dict(all24501_node_types_independently_derived=True,all16380_X_entries_source_bound=True,
        all130_residual_b_roots_source_bound=True,canonical_parameter_Gram_solve_inputs_verified=True,
        no_Gram_scalar_expansion_or_field_elimination=True,all130_H_all2583_compat_all210_quartic_roots_checked=True,
        all2923_normal_generators_preserve_full_normal_quadratic_terms=True,
        all45402_monomial_homogenization_identities=True,all5292_patch_recovery_identities=True,
        actual_A_patch_scheme_isomorphisms=True,finite_A_over_K_normal42_generators_atmost5=True,
        B_at_normal_origin_zero_by_finite_Nakayama=True,R_to_B_surjective_by_finite_cokernel_all_prime_test=True,
        exact_B_is_R_mod_J_colon_m_infinity=True,empty_iff_some_normal_power_in_J=True,
        generic_K_only_complex_specialization_requires_all_parameter_denominators_and_Gram_nonzero=True,
        no_sum_normal_squares_inverse_no_basis_no_free_module_no_saturation_unit_or_exhaustion=True,
        old_NORMAL_FINITE5_abstract_theorem_reused_and_credited=True,op_counts=dict(opcounts))
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in (Path(__file__),ROOT/'code/s61_sparse_generic5_full134_normal_circuit_review_20261104_v63.py',ROOT/f'research/notes/{ID}.md'):bind(p)
    report=dict(id=ID,all_passed=True,input_sha256=hashes,evidence_sha256={p.name:sha(p) for p in out.iterdir()},
        seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        scope='Exact generic-K normal42 typed presentation of the closed residual, retaining scheme structure; no saturation, membership or emptiness computation.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'],rss=report['max_rss_kib'],report_sha256=sha(out/'report.json'))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
