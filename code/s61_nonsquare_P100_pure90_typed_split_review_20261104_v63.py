#!/usr/bin/env python3
"""Independent typed/source and universal algebra audit of the global pure90 split."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_P100_PURE90_TYPED_SPLIT_20261104_v63_REVIEW'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);hashes={};cache={}
    def read(p):
        name=str(p.relative_to(ROOT));hashes[name]=sha(p)
        if name not in cache:cache[name]=json.loads(p.read_text())
        return cache[name]
    def gate(stem,standard=True):
        p=ROOT/'reports'/stem/'report.json';r=read(p)
        if standard:assert r['all_passed']
        for name,h in r.get('input_sha256',{}).items():assert sha(ROOT/name)==h,name;hashes[name]=h
        for name,h in r.get('evidence_sha256',{}).items():
            q=p.parent/name;assert sha(q)==h,name;hashes[str(q.relative_to(ROOT))]=h
        return r
    for stem in ['S61_NONSQUARE_P100_PURE15_COVER_20261104_v63_REVIEW',
        'S61_NONSQUARE_P100_H1_20261104_v63_REVIEW',
        'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63_REVIEW',
        'S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63_REVIEW',
        'S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63_REVIEW',
        'S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW']:
        gate(stem)
    b=ROOT/'reports';stem='S61_NONSQUARE_P100_PURE90_TYPED_SPLIT_20261104_v63'
    r=gate(stem,False);assert r['all_shape_and_source_bindings_passed'] and r['new_explicit_coefficient_terms']==0
    op=read(b/stem/'operator.json');assert sha(b/stem/'operator.json')==r['operator_sha256']=='9524cf9b5d8899d8ea174ab2f3dc7d4d05632b80380f544ea8d7a840382a184b'
    assert op['explicit_new_polynomial_or_RF_coefficient_terms']==0 and op['external_Bezout_atoms_remain_external']
    nodes=op['nodes'];assert len(nodes)==r['typed_nodes']==101
    references=[]
    def resolve(ref):
        p=ROOT/ref['path'];value=read(p);assert hashes[ref['path']]==ref['sha256']
        for key in ref['json_pointer'].split('/')[1:]:
            key=key.replace('~1','/').replace('~0','~')
            value=value[int(key)] if isinstance(value,list) else value[key]
        return value
    def visit(v):
        if isinstance(v,dict):
            if {'path','sha256','json_pointer'}<=v.keys():resolve(v);references.append(v)
            else:
                for child in v.values():visit(child)
        elif isinstance(v,list):
            for child in v:visit(child)
    visit(op);assert len(references)==75
    shapes=[]
    leaf_ops={'external_row_matrix','external_polynomial_scalar','external_lazy_full_weight4_source_matrix',
        'arbitrary_regular_row','external_certified_unit_scalar','lazy_polynomial_then_NF_pure_part',
        'external_lazy_scalar_root_row','certified_regular_unit_formula','external_regular_RF_matrix',
        'certified_sparse_source_embedding'}
    for i,n in enumerate(nodes):
        assert n['id']==i;shape=n['shape'];assert len(shape)==2 and min(shape)>0
        args=n.get('args',[]);assert all(0<=j<i for j in args)
        o=n['op']
        if o=='identity_matrix':assert shape[0]==shape[1]
        elif o=='zero_matrix':pass
        elif o in leaf_ops:pass
        elif o=='matrix_product':assert len(args)==2 and shapes[args[0]][1]==shapes[args[1]][0] and shape==[shapes[args[0]][0],shapes[args[1]][1]]
        elif o=='central_scalar_multiple':assert shapes[args[0]]==[1,1] and shape==shapes[args[1]]
        elif o=='matrix_sum':assert args and all(shapes[j]==shape for j in args)
        elif o=='negate':assert len(args)==1 and shape==shapes[args[0]]
        elif o=='select_columns':
            c=n['columns'];assert len(args)==1 and all(0<=k<shapes[args[0]][1] for k in c) and len(c)==len(set(c));assert shape==[shapes[args[0]][0],len(c)]
        elif o=='select_rows':
            rows=n['rows'];assert len(args)==1 and all(0<=k<shapes[args[0]][0] for k in rows) and len(rows)==len(set(rows));assert shape==[len(rows),shapes[args[0]][1]]
        elif o=='weighted_source_selector':
            assert shape==[8,15] and len(n['selected_rows'])==len(set(n['selected_rows']))==8
            assert all(0<=j<15 for j in n['selected_rows']) and len(n['row_unit_powers5'])==8
            assert all(len(p)==5 and all(isinstance(e,int) and e>=0 for e in p) for p in n['row_unit_powers5'])
            assert n['unit_factor_names']==['z','z+20','z+60','203z+2700','z+12']
        elif o=='polynomial_adjugate_8x8':assert shape==[8,8]==shapes[args[0]]
        elif o=='certified_unit_inverse':assert shape==[1,1]==shapes[args[0]]
        elif o=='transpose':assert shape==shapes[args[0]][::-1]
        elif o=='horizontal_block':assert all(shapes[j][0]==shape[0] for j in args) and sum(shapes[j][1] for j in args)==shape[1]
        elif o=='vertical_block':assert all(shapes[j][1]==shape[1] for j in args) and sum(shapes[j][0] for j in args)==shape[0]
        elif o=='block_matrix':
            assert sum(n['row_partition'])==shape[0] and sum(n['column_partition'])==shape[1]
            for a,row in enumerate(n['blocks']):
                for bb,j in enumerate(row):assert j<i and shapes[j]==[n['row_partition'][a],n['column_partition'][bb]]
        elif o=='accepted_weighted4_coordinate_transport':assert len(args)==1 and shape==shapes[args[0]]==[319,239]
        elif o=='apply_external_linear_operator':assert len(args)==1 and shapes[args[0]]==[319,239] and shape==[319,122]
        elif o=='delta7_source_preconditioner':assert shape==[112,112] and n['delta']<i and shapes[n['delta']]==[1,1]
        else:raise AssertionError(o)
        shapes.append(shape)
    ex=op['exports'];context=op['context']
    M,A,B,K,KA,KB,R=[ex[k] for k in ['M15x98','A15x8','B15x90','K8x15','KA8x8','KB8x90','residual15x90']]
    assert nodes[M]['op']=='external_row_matrix' and shapes[M]==[15,98]
    actual_rows=[resolve(ref) for ref in nodes[M]['row_refs']]
    assert len(actual_rows)==15 and all(len(row)==98 for row in actual_rows)
    selected=context['selected_local98_columns8'];tail=context['remaining_local98_columns90']
    assert selected==[0,1,2,3,4,5,6,9] and tail==[i for i in range(98) if i not in selected]
    assert nodes[A]['args']==nodes[B]['args']==[M] and nodes[A]['columns']==selected and nodes[B]['columns']==tail
    cover=read(b/'S61_NONSQUARE_P100_PURE15_COVER_20261104_v63/native.json')
    binding=read(b/'S61_NONSQUARE_P100_PURE15_COVER_20261104_v63/source-binding.json')
    old=read(b/'S61_NONSQUARE_P100_PURE8_GENERIC_RANK_20261104_v63/native.json')
    source_old=old['selected_rows'];assert source_old==[0,2,3,5,7,6,1,4]
    oldpowers=binding['old_selected_source_binding']['selected_rowLCD_powers']
    assert len(op['five_minors'])==5
    terms=[]
    for n,(item,contrib) in enumerate(zip(op['five_minors'],cover['actual_minor_polynomial_Bezout_contributors'])):
        idx=item['basis_index'];assert idx==contrib['basis_index']==[0,7,9,10,15][n]
        rows=source_old[:];powers=[p[:] for p in oldpowers]
        if idx:
            cand=cover['candidates'][idx-1];light=cand['K2_light_position'];position=cand['replaced_old_selected_row_position']
            rows[position]=8+light;powers[position]=binding['new_sources_in_order'][light]['known_unit_rowLCD_powers5']
            assert item['candidate']==dict(K2_light_position=light,replace_position=position)
            minor=cand['split']['full']
        else:assert item['candidate'] is None;minor=cover['old_minor_strip']['full']
        p,s,adj,d,ai=[item[k] for k in ['weighted_P','S','Adj','D','Ai']]
        assert nodes[p]['selected_rows']==item['source_rows_in_order']==rows
        assert nodes[p]['row_unit_powers5']==item['LCDpowers_in_order']==powers
        assert nodes[s]['op']=='matrix_product' and nodes[s]['args']==[p,A]
        assert nodes[adj]['op']=='polynomial_adjugate_8x8' and nodes[adj]['args']==[s]
        assert resolve(nodes[d]['reference'])==minor and nodes[d]['certified_equal_to_determinant_of']==s
        assert resolve(nodes[ai]['reference'])==contrib['coefficient']
        terms.append((ai,adj,p))
    # Validate the exact exported expression, not only its shape.
    def formula(i):
        n=nodes[i];o=n['op']
        if o in ['matrix_product','matrix_sum','central_scalar_multiple','negate','transpose','horizontal_block','vertical_block']:
            return (o,tuple(formula(j) for j in n['args']))
        if o in ['select_columns','select_rows']:
            return (o,formula(n['args'][0]),tuple(n['columns'] if o=='select_columns' else n['rows']))
        return ('leaf',i)
    def leaf(i):return formula(i)
    def prod(a,b):return ('matrix_product',(a,b))
    def sc(a,b):return ('central_scalar_multiple',(a,b))
    def sub(a,b):return ('matrix_sum',(a,('negate',(b,))))
    kin=nodes[K];assert kin['op']=='central_scalar_multiple';ei,sumnode=kin['args']
    assert nodes[ei]['op']=='certified_unit_inverse'
    e=nodes[ei]['args'][0];assert resolve(nodes[e]['reference'])==cover['common_actual_minor_unit_polynomial']
    assert nodes[e]['unit_monomial_powers5']==[21,93,6,0,1]
    expected_K=sc(leaf(ei),('matrix_sum',tuple(sc(leaf(ai),prod(leaf(adj),leaf(p))) for ai,adj,p in terms)))
    assert formula(K)==expected_K
    assert formula(KA)==prod(formula(K),formula(A)) and formula(KB)==prod(formula(K),formula(B))
    assert formula(R)==prod(sub(('leaf',next(n['id'] for n in nodes if n['op']=='identity_matrix' and n['shape']==[15,15])),prod(formula(A),formula(K))),formula(B))
    U,Ui=[nodes[ex[k]] for k in ['target_GL98','target_GL98_inverse']]
    assert U['acts_after_column_order']==Ui['acts_after_column_order']==selected+tail
    assert U['blocks'][0][0]==Ui['blocks'][0][0] and nodes[U['blocks'][0][0]]['op']=='identity_matrix'
    assert nodes[U['blocks'][0][1]]['op']=='negate' and nodes[U['blocks'][0][1]]['args']==[KB]
    assert Ui['blocks'][0][1]==KB and U['blocks'][1]==Ui['blocks'][1]
    # Every full source and section argument is an actual accepted input.
    cat=read(b/'S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63/certificate.json')
    H=read(b/'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63/certificate.json')
    NF=read(b/'S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63/certificate.json')
    adapter=read(b/'S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63/adapter.json')
    oracle=read(ROOT/adapter['oracle_path']);assert sha(ROOT/adapter['oracle_path'])==adapter['oracle_sha256']
    catalogue=[cat['source_rows'][i] for i in cat['retained_row_indices']]
    physical=next(n for n in nodes if n['op']=='external_lazy_full_weight4_source_matrix')
    assert physical['source_catalogue_rows']==catalogue and cat['retained_shape']==[319,239]
    assert catalogue[174:286]==H['all112_cubic_source319_definitions']
    assert cat['genuine_quartic_original_rows_retained']==oracle['genuine_quartic_original_rows']
    trans=next(n for n in nodes if n['op']=='accepted_weighted4_coordinate_transport')
    nf=next(n for n in nodes if n['op']=='apply_external_linear_operator')
    assert trans['args']==[physical['id']] and nf['args']==[trans['id']]
    assert resolve(trans['old_physical_monomials_ref'])==cat['weight4_columns']
    assert resolve(trans['new_joint_monomials_ref'])==H['parent_joint_target239_monomials']
    assert len(resolve(trans['physical_light_inverse_ref']))==7
    for mons,heavy in [(cat['weight4_columns'],5),(H['parent_joint_target239_monomials'],7)]:
        assert len(mons)==len(set(map(tuple,mons)))==239
        assert all(sum(m)+(m[heavy])==4 for m in mons)
        assert [sum(m[heavy]==k for m in mons) for k in range(3)]==[210,28,1]
    assert resolve(nf['reference'])==NF['full_target_projection239_to122']
    f0,f1,f2=[next(n['id'] for n in nodes if n['op']=='select_columns' and n['args']==[nf['id']] and n['columns']==H['target122_Hdegree_classes'][str(k)]) for k in [0,1,2]]
    assert H['target122_Hdegree_classes']['0']==context['pure98_target122_labels']
    assert [context['pure98_target122_labels'][i] for i in selected]==context['selected_original122_labels8']
    chi0=next(n for n in nodes if n['op']=='lazy_polynomial_then_NF_pure_part')
    chi1=next(n for n in nodes if n['op']=='external_lazy_scalar_root_row')
    assert resolve(chi0['polynomial_ref'])==H['chi_full_joint_weight4_polynomial_root']
    assert resolve(chi0['polynomial_nodes_ref'])==H['chi_polynomial_nodes']
    assert resolve(chi0['NF_ref'])==NF['full_target_projection239_to122']
    assert resolve(chi1['reference'])==H['chi_H1_post117_scalar_roots']
    u2=next(n['id'] for n in nodes if n['op']=='external_certified_unit_scalar')
    u2i=next(n['id'] for n in nodes if n['op']=='certified_unit_inverse' and n['args']==[u2])
    assert {n['args'][0] for n in nodes if n['op']=='certified_unit_inverse'}=={e,u2}
    delta=next(n for n in nodes if n['op']=='delta7_source_preconditioner')
    assert nodes[delta['delta']]['formula']=='-64000*C^2 restrictedtoP100'
    L=next(n for n in nodes if n['op']=='external_regular_RF_matrix')
    h1=read(b/'S61_NONSQUARE_P100_H1_RIGHT_INVERSE_20261104_v63/native.json')
    assert resolve(L['reference'])==h1['right_inverse112x23']
    assert h1['all23_source_identities_both_source_frames_verified']
    pi_original=('select_columns',('select_rows',formula(nf['id']),tuple(range(174,286))),tuple(H['target122_Hdegree_classes']['0']))
    pi_actual=prod(leaf(delta['id']),pi_original)
    section=prod(('transpose',(leaf(L['id']),)),pi_actual)
    h2pure=sub(formula(f0),prod(sc(leaf(u2i),formula(f2)),leaf(chi0['id'])))
    h2h1=sub(formula(f1),prod(sc(leaf(u2i),formula(f2)),leaf(chi1['id'])))
    expected_Gamma=sub(h2pure,prod(h2h1,section))
    gamma=ex['full_Gamma319x98'];assert formula(gamma)==expected_Gamma
    expected_g90=sub(('select_columns',formula(gamma),tuple(tail)),prod(('select_columns',formula(gamma),tuple(selected)),formula(KB)))
    assert formula(ex['full_Gamma319x90'])==expected_g90
    assert formula(ex['full_redundant334x90'])==('vertical_block',(expected_g90,formula(R)))
    J=nodes[ex['source_embedding15x319']];emb=J['row_recipes'];assert len(emb)==15 and shapes[J['id']]==[15,319]
    assert catalogue[emb[0]['original319_source_index']]==['quadric',5,[0,0,0,0,0,1,0,0]] and emb[0]['original319_source_index']==148
    assert [v['original319_source_index'] for v in emb[1:8]]==list(range(174,181))
    for j,rec in enumerate(emb[8:]):
        assert rec['original319_source_indices']==[174+7*s+j for s in range(16)]
        k2=read(b/'S61_NONSQUARE_P100_PURE_CUBIC_KERNEL_20261104_v63/pure-cubic84-2.json')
        assert resolve(rec['original16_cofactor_ref'])==k2['old16_source_cofactors']
        assert resolve(rec['actual_delta7_cofactor_ref'])==k2['actual_delta7_source_cofactors']
    # Universal noncommutative word reduction using KA=I. No scalar coefficient
    # products are expanded, especially none of the external18MB cover atoms.
    def word(w):
        w=list(w)
        while True:
            p=next((i for i in range(len(w)-1) if w[i:i+2]==['K','A']),None)
            if p is None:return tuple(w)
            del w[p:p+2]
    def poly(terms):
        result={}
        for c,w in terms:
            k=word(w);result[k]=result.get(k,0)+c
            if not result[k]:del result[k]
        return result
    assert not poly([(1,('K','B')),(-1,('K','A','K','B'))])
    assert poly([(1,()),(-2,('A','K')),(1,('A','K','A','K'))])==poly([(1,()),(-1,('A','K'))])
    proof=dict(all_passed=True,all101_node_types_shapes_and_backward_dependencies_checked=True,
        all75_external_JSON_pointers_resolved_and_hashes_current=True,
        weighted_minor_selectors5_bound_to_actual_old_and_new_LCDs=True,
        global_KA_identity='E^-1 sum ai Adj(Si) Pi A = E^-1 sum ai det(Si) I8 = I8',
        inverse_only_E_and_U2=True,old_unknown_determinant_never_inverted=True,
        target_GL_unipotent_det1_and_general_row_full_tail_preserved=True,
        all319_source_catalogue_rows_and33_true_quartics_retained=True,
        physical_Tn_to_joint_before_NF117=True,H2_chi_full_pure_and_H1_parts_retained=True,
        actual_delta7_on_PiOrig_before_Lactual_transpose=True,
        full_Gamma_formula_checked=True,actual_source_embedding15_by319_checked=True,
        all15_selected_residuals_retained_no_free7_basis_assumption=True,
        redundant334_by90_same_cokernel_proof='The accepted actual NF117/H2/H1 sections remove141 target directions and give Gamma. M=J Gamma by the independently accepted q5H, literalG18 and fullK2 source identities. Appending M changes no image. Right target GL gives [A,R], and the regular source combinations K give [I8,0]. Removing precisely these8 targets leaves all319 Gamma90 and all15 R, with the same cokernel. The original weighted354-to319 reduction uses the five unit genuine relations times7 light variables, not a pointwise fullcubic assumption.',
        P100_restriction='a=-10(z+20)/(9z^2),e=(z+100)/(9z),d=-3z(z+20)/100; F=-(203z+2700)/(9z),B=z(z+60)/12,C=-25(z+60)/(9z),J=(z+12)/6 are chart units',
        source_coefficient_functions_are_accepted_external_recipes_not_evaluated_rows=True,
        no_new_explicit_coefficients_full90_evaluation_minor_rank_GB_point_or_prime=True,
        no_full239_or_exhaustion_claim=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in [Path(__file__),ROOT/f'research/notes/{ID}.md']:hashes[str(p.relative_to(ROOT))]=sha(p)
    report=dict(id=ID,all_passed=True,global_KA_leftinverse_over_whole_P100=True,
        actual_redundant334_by90_same_cokernel=True,new_coefficient_terms=0,
        no_residual90_rank_or_full239_claim=True,input_sha256=hashes,
        evidence_sha256={p.name:sha(p) for p in out.iterdir()},seconds=time.monotonic()-start,
        max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'],report_sha256=sha(out/'report.json'))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
