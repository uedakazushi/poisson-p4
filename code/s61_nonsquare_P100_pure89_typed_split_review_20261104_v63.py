#!/usr/bin/env python3
"""Independent finite-type/source and pure-algebra audit of global pure89 split.
No coefficient constructor, Bezout, CRT, norm, NF or rank is evaluated.
"""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PRIMARY='S61_NONSQUARE_P100_PURE89_TYPED_SPLIT_20261104_v63';ID=PRIMARY+'_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/PRIMARY;rp=folder/'report.json';A=read(rp)
    assert sha(rp)=='14b8b8669300bc96c172edf15a0b7ec9c5ae6645addaadc494912ad9b7793653' and A['all_structural_and_binding_checks_passed']
    files=[rp,Path(__file__)];cache={}
    def load(p):
        if p not in cache:cache[p]=read(p);files.append(p)
        return cache[p]
    for q,h in A['input_sha256'].items():assert sha(ROOT/q)==h;load(ROOT/q)
    for q,h in A['source_sha256'].items():assert sha(ROOT/q)==h;files.append(ROOT/q)
    for q,h in A['evidence_sha256'].items():assert sha(folder/q)==h;load(folder/q)
    O=load(folder/'operator.json');nodes=O['nodes'];ex=O['exports'];ctx=O['context']
    assert len(nodes)==192 and not O['typed_coefficient_operators_executed'] and O['explicit_new_parameter_polynomial_or_RF_terms']==0
    parent=ROOT/'reports/S61_NONSQUARE_P100_PURE90_TYPED_SPLIT_20261104_v63_REVIEW/report.json';pr=load(parent);assert pr['all_passed']
    for q,h in pr['input_sha256'].items():assert sha(ROOT/q)==h
    refs=[]
    def pointer(r):
        p=ROOT/r['path'];assert sha(p)==r['sha256'];v=load(p)
        for k in r['json_pointer'].split('/')[1:]:
            k=k.replace('~1','/').replace('~0','~');v=v[int(k)] if isinstance(v,list) else v[k]
        refs.append(r);return v
    shapes=[]
    for i,n in enumerate(nodes):
        assert n['id']==i and len(n['shape'])==2 and min(n['shape'])>0
        s=n['shape'];args=n.get('args',[]);assert all(0<=j<i for j in args);t=n['op']
        if t=='external_exact_operator':pointer(n['reference'])
        elif t=='transpose':assert s==shapes[args[0]][::-1]
        elif t=='matrix_product':assert shapes[args[0]][1]==shapes[args[1]][0] and s==[shapes[args[0]][0],shapes[args[1]][1]]
        elif t=='sum':assert len(args)==2 and shapes[args[0]]==shapes[args[1]]==s
        elif t=='negate':assert shapes[args[0]]==s
        elif t=='central_scalar_multiple':assert shapes[args[0]]==[1,1] and shapes[args[1]]==s
        elif t in ['universal_polynomial_adjugate','universal_polynomial_determinant']:
            assert shapes[args[0]][0]==shapes[args[0]][1] and not n['evaluated']
            assert s==(shapes[args[0]] if t.endswith('adjugate') else [1,1])
        elif t=='identity':assert s[0]==s[1]
        elif t=='coordinate_column_embedding':assert len(n['columns'])==s[1] and len(set(n['columns']))==s[1] and all(0<=j<s[0] for j in n['columns'])
        elif t=='known_unit_weighted_row_selector':assert s==[9,15] and len(n['rows'])==len(set(n['rows']))==9 and len(n['LCD_powers5'])==9
        elif t=='saved_known_unit_content_times_five_factor_product':assert s==[1,1] and n['not_unknown_factor'] and len(n['powers'])==5
        elif t=='certified_regular_unit_inverse':assert shapes[args[0]]==s==[1,1] and not n['evaluated']
        elif t=='canonical_QQ_polynomial_extendedgcd_family':assert s==[5,1] and len(args)==5 and all(shapes[j]==[1,1] for j in args) and not n['evaluated'] and n['not_a_saved_QQ_Bezout']
        elif t=='component_polynomial':assert s==[1,1] and shapes[args[0]]==[5,1] and 0<=n['index']<5
        elif t=='canonical_integer_lift_Fp_CRT_idempotent':assert s==[1,1] and n['prime']==32003 and not n['evaluated'] and n['not_a_saved_idempotent']
        elif t=='canonical_multiplication_mod_monic_polynomial':assert s==[15,15] and all(shapes[j]==[1,1] for j in args) and not n['evaluated']
        elif t=='constant_standard_column':assert s==[15,1] and n['index']==0
        elif t=='coefficient_column_to_QQ_polynomial':assert s==[1,1] and shapes[args[0]]==[15,1] and not n['evaluated']
        elif t=='certified_exact_division_in_curve_unit_ring':assert s==[1,1] and all(shapes[j]==s for j in args) and not n['evaluated']
        elif t=='unit_lower_block_frame':assert s==[98,98] and shapes[args[0]]==[89,9] and n['determinant']==1
        elif t=='vertical_block':assert sum(shapes[j][0] for j in args)==s[0] and all(shapes[j][1]==s[1] for j in args)
        elif t in ['constant_one','saved_Fp_monic_polynomial']:assert s==[1,1]
        else:raise AssertionError(t)
        shapes.append(s)
    def node(i,t,args=None):
        n=nodes[i];assert n['op']==t
        if args is not None:assert n['args']==args
        return n
    def flatten(i):
        return flatten(nodes[i]['args'][0])+flatten(nodes[i]['args'][1]) if nodes[i]['op']=='sum' else [i]
    gamma=ex['Gamma319x98'];aa=ex['A98x319'];pi=ex['Pi9x98'];ec=ex['E_C'];et=ex['E_T'];ff=ex['F9x319']
    node(aa,'transpose',[gamma]);node(pi,'transpose',[ec]);node(ff,'matrix_product',[pi,aa])
    C9=ctx['C9_oldlocal98'];T89=ctx['T89_oldlocal98']
    assert C9==[0,1,2,3,4,5,6,9,25] and T89==[i for i in range(98) if i not in C9]
    assert node(ec,'coordinate_column_embedding')['columns']==C9 and node(et,'coordinate_column_embedding')['columns']==T89
    parentop=load(ROOT/'reports/S61_NONSQUARE_P100_PURE90_TYPED_SPLIT_20261104_v63/operator.json')
    pctx=parentop['context'];assert pointer(nodes[gamma]['reference'])==parentop['nodes'][parentop['exports']['full_Gamma319x98']]
    assert ctx['Gamma_formula']==pctx['general_full_source_formula'] and 'NOTzero' in ctx['Gamma_formula']
    assert pctx['H1_source_frame'].startswith('UsePiActual=delta7*PiOrig')
    assert pctx['source_embeddings15'][0]['original319_source_index']==148
    assert [r['original319_source_index'] for r in pctx['source_embeddings15'][1:8]]==list(range(174,181))
    assert [pctx['pure98_target122_labels'][i] for i in C9]==ctx['C9_original122'] and ctx['C9_original122'][-1]==30
    oldsrc=load(ROOT/'reports/S61_NONSQUARE_P100_NEXT9_COVER_20261104_v63/source-binding.json')
    old=load(ROOT/'reports/S61_NONSQUARE_P100_NEXT9_COVER_20261104_v63/native.json')
    gcd=load(ROOT/'reports/S61_NONSQUARE_P100_NEXT9_QQ_GCD_20261104_v63/native.json')
    fac=load(ROOT/'reports/S61_NONSQUARE_P100_H15_MODP_FACTOR_INVENTORY_20261104_v63/native.json')
    assert gcd['gcd_degree']==15 and gcd['all5_literal_gcd_times_quotient_equals_input_verified'] and gcd['literal_modular_Bezout_verified'] and gcd['cofactor_reduction_degrees_preserved'] and gcd['oldprime_cofactor_gcd_is1']
    assert fac['Hbar_squarefree'] and fac['full_product_equals_Hbar'] and all(r['multiplicity']==1 for r in fac['factors'])
    RH=ex['R_H'];parts=flatten(RH);assert len(parts)==5
    inverse_inputs=set();betas=[]
    for j,(binding,case,d,part) in enumerate(zip(ctx['old_five_bindings'],oldsrc['five_cases'],old['determinants5'],parts)):
        assert case['nine_local98_target_columns']==C9
        selector=node(binding['source_selector'],'known_unit_weighted_row_selector')
        assert selector['rows']==case['nine_source_rows_in_order'] and selector['LCD_powers5']==case['effective_LCD_powers5']
        unit=node(binding['unit'],'saved_known_unit_content_times_five_factor_product');assert F(d['QQ_unit_content'])!=0 and unit['powers']==d['known_unit_multiplicities5']
        assert unit['QQ_content_reference']==f'/determinants5/{j}/QQ_unit_content'
        X=binding['actual_projected_matrix'];node(X,'matrix_product');assert nodes[X]['args'][0]==ff
        si=nodes[X]['args'][1];node(si,'matrix_product');jt,wit=nodes[si]['args'];assert nodes[jt]['op']=='transpose' and nodes[wit]['args']==[selector['id']]
        jref=nodes[nodes[jt]['args'][0]]['reference'];assert pointer(jref)==parentop['nodes'][parentop['exports']['source_embedding15x319']]
        node(binding['determinant'],'universal_polynomial_determinant',[X])
        partn=node(part,'central_scalar_multiple');coef,lift=partn['args'];node(coef,'matrix_product');beta,iU=nodes[coef]['args']
        bn=node(beta,'component_polynomial');assert bn['index']==j;betas.append(bn['args'][0]);node(iU,'certified_regular_unit_inverse',[unit['id']]);inverse_inputs.add(unit['id'])
        node(lift,'matrix_product');assert nodes[lift]['args'][0]==si;node(nodes[lift]['args'][1],'universal_polynomial_adjugate',[X])
        assert pointer(nodes[binding['quotient']]['reference'])==gcd['all5_quotients'][j]
    assert len(set(betas))==1
    beta=node(betas[0],'canonical_QQ_polynomial_extendedgcd_family');assert beta['args']==[r['quotient'] for r in ctx['old_five_bindings']]
    G=ex['G'];finite_terms=flatten(G);r9=node(ex['R9'],'sum');assert len(finite_terms)==7
    # Rfinite is the second term's coefficient-independent matrix factor.
    second=node(r9['args'][1],'central_scalar_multiple');b=second['args'][0];Rfin=second['args'][1];rparts=flatten(Rfin);assert len(rparts)==7
    for j,(binding,term,rpart) in enumerate(zip(ctx['finite_seven_bindings'],finite_terms,rparts)):
        r=load(ROOT/f'reports/S61_NONSQUARE_P100_H15_FULL239_FINITE_TEST_20261104_v63/factor{j}.json')
        assert r['rank']==98 and len(set(r['selected_source319_rows']))==98
        sj=binding['selected_source_embedding'];node(sj,'coordinate_column_embedding');assert nodes[sj]['columns']==r['selected_source319_rows']
        bj=binding['global_selected_matrix'];node(bj,'matrix_product',[aa,sj]);dj=binding['global_Dj'];node(dj,'universal_polynomial_determinant',[bj])
        ej=binding['integer_CRT_lift'];en=node(ej,'canonical_integer_lift_Fp_CRT_idempotent');assert en['factor_index']==j and en['factor']==fac['factors'][j]['coefficients']
        hb=node(en['args'][0],'saved_Fp_monic_polynomial');assert hb['coefficients']==fac['Hbar']
        node(term,'matrix_product',[ej,dj]);node(rpart,'central_scalar_multiple',[ej,binding['projected_adjugate_source_lift']])
        lift=node(binding['projected_adjugate_source_lift'],'matrix_product');assert lift['args'][1]==ec
        sa=node(lift['args'][0],'matrix_product');assert sa['args'][0]==sj;node(sa['args'][1],'universal_polynomial_adjugate',[bj])
    norm=ex['norm15'];normn=node(norm,'universal_polynomial_determinant');MG=normn['args'][0]
    # The modulus is the exact monic H, not any unknown determinant.
    mg=node(MG,'canonical_multiplication_mod_monic_polynomial');assert mg['args'][0]==G and pointer(nodes[mg['args'][1]]['reference'])==gcd['exact_gcd_monic']
    bpoly=node(ex['b'],'coefficient_column_to_QQ_polynomial');bvec=node(bpoly['args'][0],'central_scalar_multiple');inorm=bvec['args'][0]
    node(inorm,'certified_regular_unit_inverse',[norm]);inverse_inputs.add(norm)
    prod=node(bvec['args'][1],'matrix_product');node(prod['args'][0],'universal_polynomial_adjugate',[MG]);node(prod['args'][1],'constant_standard_column')
    a=node(ex['a'],'certified_exact_division_in_curve_unit_ring');assert a['args'][1]==mg['args'][1]
    numerator=node(a['args'][0],'sum');node(numerator['args'][0],'constant_one');ng=node(numerator['args'][1],'negate');node(ng['args'][0],'matrix_product',[ex['b'],G])
    first=node(r9['args'][0],'central_scalar_multiple',[ex['a'],RH]);assert second['args'][0]==ex['b']
    assert {n['args'][0] for n in nodes if n['op']=='certified_regular_unit_inverse'}==inverse_inputs
    node(ex['F_R9'],'matrix_product',[ff,ex['R9']])
    AR=next(n['id'] for n in nodes if n['op']=='matrix_product' and n['args']==[aa,ex['R9']])
    targetP=node(ex['target_projection'],'sum');node(targetP['args'][0],'identity');pn=node(targetP['args'][1],'negate');node(pn['args'][0],'matrix_product',[AR,pi])
    residual=node(ex['residual89x319'],'matrix_product');assert residual['args'][1]==aa
    projected=node(residual['args'][0],'matrix_product');assert projected['args'][1]==ex['target_projection'];node(projected['args'][0],'transpose',[et])
    node(ex['residual319x89'],'transpose',[ex['residual89x319']])
    assert ctx['all319_sources_and_pure_tails_retained'] and O['canonical_QQ_Bezout_CRT_norm_and_exact_division_are_unexpanded']
    proof={'all192_node_types_and_sources_checked':True,'external_references_resolved':len(refs),
      'fixed_C9':[0,1,2,3,4,5,6,9,25],'ninth_original122_target':30,
      'oldfive_actual_weighted_selectors_QQcontent_and_cofactor_gcd_bound':True,
      'finite7_selected_actual_global98matrices_source_bound':True,
      'constructor_proof':['Canonical QQ polynomial extendedgcd of certified coprime Qi gives beta; universal adjugates give FRH=H I9.',
        'Squarefree Hbar seven-factor CRT integerlifts give Gbar unit using actual seven selected98 minors; all source/LCD coefficients reduce under the accepted finitefield gate.',
        'Thus the15-dimensional multiplication norm of G is a nonzeroQQconstant. Its adjugate gives bG=1 modH.',
        'After clearing only known curve-unit denominators D, gcd(H,D)=1 makes the numerator of1-bG exactly divisible byH. This defines regular a without invertingH.',
        'R9=aRH+bRfinite satisfies FR9=(aH+bG)I9=I9.',
        'P=I-AR9Pi projects onto kerPi and preserves imageA; PiT P and ET induce inverse cokernel maps. All319source generators remain.'],
      'only_global_inverse_inputs':'five known curve units/nonzeroQQcontents and norm15 constant',
      'gA_delta_and_full_pure_tails_inherited_exactly_from_accepted_PURE90_source_operator':True,
      'ordered_target_frame_det1_no_free310_basis':True,'samecoker319_by89':True,
      'no_executed_QQ_Bezout_CRT_norm_coefficient_NF_rank_or_full239_claim':True}
    ep=out/'independent-typed-and-patching-proof.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
    files.append(ROOT/'research/notes'/f'{ID}.md')
    report={'id':ID,'all_passed':True,'global_fixed9_projection_and_finite_typed_R9_accepted':True,
      'samecoker319_by89_accepted':True,'coefficient_operators_not_executed':True,'no_full239_or_free310_claim':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},'evidence_sha256':{ep.name:sha(ep)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['id','all_passed','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
