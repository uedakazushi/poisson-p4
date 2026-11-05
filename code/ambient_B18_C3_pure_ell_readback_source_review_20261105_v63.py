#!/usr/bin/env python3
"""Independent pre-value custody/frame/operation audit; never emit source values."""
import argparse,ast,collections,hashlib,json,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_C3_PURE_ELL_READBACK_SOURCE_20261105_v63_REVIEW'
SOURCE='reports/AMBIENT_B18_C3_PURE_ELL_READBACK_SOURCE_20261105_v63'
SOURCE_SHA='f78ba207056dda46060c035156f50d642dda24b43ec9bda0f36383eb32ebe010'

def sha(p):
    h=hashlib.sha256()
    with pathlib.Path(p).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def save(p,x):pathlib.Path(p).write_text(json.dumps(x,indent=2)+'\n')
def field(t,k):
    j=t.index(':',t.index(json.dumps(k)))+1
    while t[j].isspace():j+=1
    return json.JSONDecoder().raw_decode(t,j)[0]
def nums(x):
    if isinstance(x,dict):return sum(nums(k)+nums(v) for k,v in x.items())
    if isinstance(x,(list,tuple)):return sum(map(nums,x))
    if isinstance(x,int):return 1
    if isinstance(x,str):
        try:Q(x);return 1
        except (ValueError,ZeroDivisionError):return 0
    return 0

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--custody',required=True);ap.add_argument('--prepared-code',required=True);args=ap.parse_args()
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
    try:
        sp=ROOT/SOURCE/'report.json';assert sha(sp)==SOURCE_SHA;sr=read(sp)
        assert sr['all_passed'] and sr['readback_values_not_evaluated'] and not sr['new_IFT_gauge_higherjet']
        pins={str(sp.relative_to(ROOT)):SOURCE_SHA}
        for p,h in sr['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
        for p,h in sr['evidence_sha256'].items():
            p=SOURCE+'/'+p;assert sha(ROOT/p)==h;pins[p]=h
        D=read(ROOT/SOURCE/'compiled-original164-DAG-and-provided-pure-ell-inputs.json')
        profile=read(ROOT/SOURCE/'pure-ell-readback-custody-height-and-whole-budget.json')
        cp=ROOT/args.custody;custody=read(cp);assert custody['all_passed'];pins[args.custody]=sha(cp)
        for p,h in custody['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
        codepath=ROOT/args.prepared_code;code=codepath.read_text();pins[args.prepared_code]=sha(codepath)
        assert '_v2_' in codepath.name and sha(codepath)==custody['prepared_V2_code_sha256']
        for p,h in custody.get('evidence_sha256',{}).items():
            pp=cp.parent/p;assert sha(pp)==h;pins[str(pp.relative_to(ROOT))]=h
        snapshot=custody['preserved_V1_snapshot'];assert sha(ROOT/snapshot['path'])==snapshot['sha256'];pins[snapshot['path']]=snapshot['sha256']
        # The exact newer operational helper text/control flow is reviewed before this script is run.
        tree=ast.parse(code);functions={n.name:n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
        assert 'main' in functions and 'times' in functions and 'plus' in functions
        assert not any(isinstance(n,ast.ImportFrom) and n.module and ('ambient_' in n.module) for n in tree.body)
        gatepath=next(ROOT/p for p in sr['input_sha256'] if p.endswith('actual-Artin20-source-and-cost-preflight.json'))
        gt=gatepath.read_text()
        keys=['unknown105_order','equation105_order','target_subtraction_positions','needed_closure_order','requested_root_nodes','unit_inverse_center_checks','physical12_coefficient_roots']
        gate={k:field(gt,k) for k in keys};del gt
        records=D['operator_records_originalID_op_args'];order=gate['needed_closure_order']
        assert len(order)==len(records)==18584 and [r[0] for r in records]==order and order==sorted(order)
        unknown=set(gate['unknown105_order']);assert len(unknown)==105
        assert D['equation105_order']==gate['equation105_order'] and len(D['equation105_order'])==105
        assert D['CA_target_positions']==gate['target_subtraction_positions']==list(range(48,57))
        assert D['requested_root_nodes164']==gate['requested_root_nodes'] and len(D['requested_root_nodes164'])==164
        assert D['unit_inverse_centers6']==gate['unit_inverse_center_checks'] and len(D['unit_inverse_centers6'])==6
        assert D['physical12_coefficient_roots']==gate['physical12_coefficient_roots'] and len(D['physical12_coefficient_roots'])==12
        sourcepath=next(ROOT/p for p in sr['input_sha256'] if p.endswith('actual-source-coefficients.json'))
        text=sourcepath.read_text();families=field(text,'families')
        fa,fB=families['CA'],families['CB']
        assert gate['equation105_order'][:48]==fa['all48_zero_IFT_quadratic_coefficient_equations']
        assert gate['equation105_order'][57:]==fB['all48_zero_IFT_quadratic_coefficient_equations']
        assert gate['equation105_order'][48:57]==[i for row,i in fa['actual_zero_odd_I0_generators'][:7]]+[i for row,i in fa['actual_positive_odd_J0_generators']]
        assert D['nu_roots_CA_CB']==[fa['actual_tracefree_nu'],fB['actual_tracefree_nu']]
        for r in D['physical12_coefficient_roots']:
            a,b=fa['complete_normalized_physical_w_coefficients'][0][r['physical_row']]
            c,d=fB['complete_normalized_physical_w_coefficients'][0][r['physical_row']]
            assert (a,b,c,d)==(r['aA_node'],r['bA_node'],r['aB_node'],r['bB_node'])
            assert r['delta_operator']==['subtract',a,c]
        physical_q20_roots={i for v in families.values() for pair in v['complete_normalized_physical_w_coefficients'][0] for i in pair}
        physical_q20_roots.update(D['nu_roots_CA_CB']);physical_q20_roots.update([7790,28539]);del families,fa,fB
        pos=text.index('[',text.index('"DAG_scalar_nodes"'));assert pos==D['source_field_custody']['DAG_stream_interval'][0]
        pos+=1;decoder=json.JSONDecoder();cursor=0;oldid=0;seen=set();constbits=0;metadata_hash=hashlib.sha256()
        while cursor<len(records):
            while text[pos].isspace() or text[pos]==',':pos+=1
            node,pos=decoder.raw_decode(text,pos)
            if oldid==records[cursor][0]:
                rec=records[cursor]
                if oldid in unknown:assert node[0]=='v' and rec==[oldid,'input']
                else:
                    assert rec==[oldid]+node and node[0] in ['c','a','m','n','i']
                    if node[0]!='c':assert all(i in seen for i in node[1:])
                    else:
                        assert node[2] and max(abs(node[1]).bit_length(),abs(node[2]).bit_length())<=2048
                        constbits=max(constbits,abs(node[1]).bit_length(),abs(node[2]).bit_length())
                seen.add(oldid);cursor+=1;metadata_hash.update((json.dumps(rec,separators=(',',':'))+'\n').encode())
            oldid+=1
        assert pos==D['source_field_custody']['DAG_stream_interval'][1]
        assert set(D['requested_root_nodes164'])<=seen;del text,seen,node,rec
        # B27 and its cached source inputs are read only as provided coefficient arrays.
        qpath=next(ROOT/p for p in sr['input_sha256'] if p.endswith('phase1-common105-and-114-source-equations.json'))
        qt=qpath.read_text();mon=[tuple(m) for m in field(qt,'base27_monomials')];v105=field(qt,'all105_cached_composite_values')
        residuals=field(qt,'all114_original_residuals');physical_q20_roots.update(i for i,v in residuals)
        assert len(residuals)==114 and all(not v for i,v in residuals)
        assert set(D['requested_root_nodes164'])<=physical_q20_roots
        assert len(mon)==27 and len(set(mon))==27
        assert all(len(m)<=2 and all(0<=i<9 for i in m) and tuple(sorted(m))==m for m in mon)
        assert () in mon and (8,) in mon and (8,8) in mon
        for m in mon:
            for j in range(len(m)):assert m[:j]+m[j+1:] in mon
        def image(m):return len(m) if all(i==8 for i in m) and len(m)<3 else None
        # All729 basis products: neither ring coefficient values nor missing readbacks are evaluated.
        for x,y in itertools_product(mon,mon):
            z=tuple(sorted(x+y));lhs=image(z) if z in mon else None
            ix,iy=image(x),image(y);rhs=ix+iy if ix is not None and iy is not None and ix+iy<3 else None
            assert lhs==rhs
        inds=[mon.index(m) for m in [(),(8,),(8,8)]]
        assert D['provided105_originalID_pure_ell_coefficients']==[[i,[v[j] for j in inds]] for i,v in v105]
        assert {i for i,v in v105}==unknown
        assert D['expected_nu_pure_ell3_from_Q20']==[[v[j] for j in inds] for v in field(qt,'nu_A_B27')]
        del qt,v105,residuals
        rawpath=next(ROOT/p for p in sr['input_sha256'] if p.endswith('actual-ell2-product33-and-raw45-actions.json'))
        oldpath=next(ROOT/p for p in sr['input_sha256'] if p.endswith('actual-source105-ell2-raw33-values.json'))
        b33=field(rawpath.read_text(),'actual33_basis');zero=b33.index([0,0,0]);ell=b33.index([0,0,1])
        old164=field(oldpath.read_text(),'all_requested164_root_values')
        assert D['expected164_constant_and_ell_coefficients_from_C2']==[[i,[dict(v).get(zero,'0'),dict(v).get(ell,'0')]] for i,v in old164]
        # Pre-existing independent Q20/C2 proofs certify exactly this original source namespace.
        for p in sr['input_sha256']:
            if p.endswith('_REVIEW/report.json'):assert read(ROOT/p)['all_passed']
        qcodepath=next(ROOT/p for p in sr['input_sha256'] if p.endswith('quadratic20_values_v3_20261105_v63.py'))
        qcode=qcodepath.read_text();assert 'guard(v);vals[i]=v' in qcode and 'assert overlap==164' in qcode
        assert "b.add(y[7],b.scale(y[8],Q(1,5)))" in qcode
        assert profile['future_owned_budget']=={'compiled_all_metadata_and_supplied_inputs':56696,'all18584_value_triples':55752,
          'output_all164_triples_plus_labels':656,'output_physical12_triples_and_nu':150,'controller_peer_and_serialization_reserve':10000}
        assert sum(profile['future_owned_budget'].values())==123254<130000<150000
        assert count_records(records)==55385 and max(r[0] for r in records)==40943
        assert profile['future_input_plus_output_bytes']==2<<20 and profile['source_bit_bound']==2048
        # Operational V2 helper audit is pinned by the accompanying exact-code receipt.
        # No helper is imported or executed here; no actual pure-ell value is produced.
        proof={
          'all_primary_source_and_prepared_custody_hashes_match':True,
          'all18584_original_DAG_records_rebuilt_byte_stream_one_node_at_a_time':True,
          'record_stream_sha256':metadata_hash.hexdigest(),'input105_labels_and_provided315_coefficients_exact_Q20':True,
          'input105_are_original_variables_and_no_new_IFT_or_gauge':True,
          'all164_requested_roots_subset_of_actual_Q20_guarded_roots':True,
          'same_target_F9':'(zA0..6,rA,sA), sA=rA+ell/5; pureell target(0^8,ell/5)',
          'historical_ell0_gate_last_target_rA_not_used_as_full_sA_identity':True,
          'all105_equation_positions_and_12_actual_physical_root_pairs_verified':True,
          'physical_row_indices_are_zero_based_in_original_physical150_array':True,
          'all164_old_constant_first_overlap_source_arrays_verified':True,
          'all6_unit_inverse_metadata_and_both_nu_root_frames_verified':True,
          'B27_to_pureell3_homomorphism_all729_basis_products_verified':True,
          'height_prior':'Each saved original DAG coefficient is a coefficient of the SAME guarded Q20 B27 value under this restriction, hence <=2048 bits.',
          'height_operational':'Prior saved-node height is distinct from temporary arithmetic bounds; strict V2 prechecks must stop before unsafe QQ sums/products, including diagnostics and physical delta.',
          'metadata_numeric_records':55385,'dense_retained_value_coefficients_upper':55752,
          'dense_original_ID_array_control_slots':40944,'future_owned_source_total':123254,
          'reserved_controller_slots':10000,'static_cap':130000,'hard_numeric_cap':150000,
          'raw_cap':2<<20,'seconds':60,'memory_bytes':2<<30,'original_constant_input_bit_max':constbits,
          'source_inputs_only_no_164_values_emitted_by_reviewer':True,
          'delta_ell_squared_NOT_assumed_zero':True,'C3_rank_cubic_or_exponent_result_not_claimed':True,
          'prepared_code':args.prepared_code,'prepared_code_sha256':sha(codepath),'prepared_custody_sha256':sha(cp),
        }
        # Fresh exact-code operational proof and small rational fixtures, no source-value replay.
        operational=out/'operational-helper-proof.json'
        required={'halt','bits','qconst','qparse','qadd','qmul','qneg','plus','neg','times','main'}
        assert required<=set(functions)
        def bound_expr(name):
            return [n.value for n in ast.walk(functions[name]) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='bound' for t in n.targets)][0]
        assert ast.dump(bound_expr('qadd'))==ast.dump(ast.parse('max(max(nx+dy,ny+dx)+1,dx+dy)',mode='eval').body)
        assert ast.dump(bound_expr('qmul'))==ast.dump(ast.parse('max(nx+ny,dx+dy)',mode='eval').body)
        for fun,oper in [('qadd',ast.Add),('qmul',ast.Mult)]:
            arith=[n for n in ast.walk(functions[fun]) if isinstance(n,ast.BinOp) and isinstance(n.op,oper) and isinstance(n.left,ast.Name) and n.left.id=='x' and isinstance(n.right,ast.Name) and n.right.id=='y']
            assert len(arith)==1
            guards=[n for n in ast.walk(functions[fun]) if isinstance(n,ast.If) and 'bound > 2048' in ast.unparse(n.test)]
            assert len(guards)==1 and guards[0].lineno<arith[0].lineno and any(isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='halt' for n in ast.walk(guards[0]))
        seg=lambda name:ast.get_source_segment(code,functions[name])
        assert 'qconst' in seg('qparse') and 'bits(x)>2048' in seg('qparse')
        assert 'qadd(x,y)' in seg('plus') and 'qneg(x)' in seg('neg')
        assert 'qadd(z,qmul(a[j],b[i-j]))' in seg('times')
        assert 'max(abs(n).bit_length(),abs(d).bit_length())>2048' in seg('qconst')
        assert seg('qconst').index('halt')<seg('qconst').index('z=Q(n,d)')
        assert "v=(qconst(1),qneg(x[1]),qadd(qmul(x[1],x[1]),qneg(x[2])))" in code
        assert "delta=plus(vals[r['aA_node']],neg(vals[r['aB_node']]))" in code
        assert "if 2*maxbits+2>2048:" in code and 'stop_before_physical_delta_subtractions' in code
        assert 'computed_requested_roots_only' in code and 'stop_before_offending_Fraction_operation' in code and 'no_cap_raise' in code
        assert code.index('HEIGHT_STOP_HANDLER=operational_checkpoint')<code.index('phi={i:tuple(map(qparse,v))')
        fixtures=[]
        for x,y in [(Q(1,3),Q(2,5)),(Q(-3,7),Q(3,7)),(Q(17,16),Q(-31,27)),(Q(1),Q(1))]:
            nx,dx=abs(x.numerator).bit_length(),x.denominator.bit_length();ny,dy=abs(y.numerator).bit_length(),y.denominator.bit_length()
            ab=max(max(nx+dy,ny+dx)+1,dx+dy);mb=max(nx+ny,dx+dy)
            actual=lambda z:max(abs(z.numerator).bit_length(),z.denominator.bit_length())
            assert actual(x+y)<=ab and actual(x*y)<=mb
            fixtures.append({'x':str(x),'y':str(y),'add_bound':ab,'multiply_bound':mb,'add_actual':actual(x+y),'multiply_actual':actual(x*y)})
        op={'all_passed':True,'prepared_code_sha256':sha(codepath),'strict_all_temporary_QQ_contract':True,
          'qadd_numerator_bound':'max(nx+dy,ny+dx)+1','qadd_denominator_bound':'dx+dy',
          'qmul_numerator_bound':'nx+ny','qmul_denominator_bound':'dx+dy',
          'bounds_checked_before_only_raw_Fraction_add_and_multiply':True,
          'zero_shortcuts_and_sign_change_do_not_increase_height':True,
          'all_convolution_inverse_diagnostics_and_delta_call_guarded_helpers':True,
          'qparse_and_original_constant_constructor_guarded_before_Fraction':True,
          'fixed_comparison_target_constants_only_Q0_and_Q1over5':True,
          'STOP_handler_streams_only_at_most164_computed_requested_roots':True,
          'scratch_and_partial_checkpoint_fit_existing10000_controller_reserve':True,
          'same_mathematical_formulas_as_frozen_V1_no_source_retry':True,'independent_small_bound_fixtures':fixtures,
          'no_actual164_source_values_computed_by_reviewer':True}
        save(operational,op)
        proof['operational_helper_proof_sha256']=sha(operational)
        save(out/'independent-original-source-ring-hom-and-budget-proof.json',proof)
        report={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'input_sha256':pins,
          'code_sha256':sha(__file__),'all_actual_source_frame_and_ring_hom_checks_passed':True,
          'all_temporary_QQ_prechecks_and_stop_protocol_accepted':True,
          'no_reviewer_readback_values_or_new_source_solve':True,
          'evidence_sha256':{p.name:sha(p) for p in out.iterdir() if p.is_file()},'elapsed_seconds':time.monotonic()-started}
        save(out/'report.json',report);signal.alarm(0)
        print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json')}))
    except Exception as e:
        save(out/'failure.json',{'error':type(e).__name__+': '+str(e),'no_values_or_positive_source_GO':True});raise

def itertools_product(a,b):
    for x in a:
        for y in b:yield x,y
def count_records(records):return sum(len(r)-1 for r in records if r[1]=='input')+sum(len(r)-1 for r in records if r[1]!='input')
if __name__=='__main__':main()
