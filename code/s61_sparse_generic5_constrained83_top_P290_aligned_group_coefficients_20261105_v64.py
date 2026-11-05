#!/usr/bin/env python3
"""ONE authorized1--4-row common-P coefficient group, thin-source streaming."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,gc,hashlib,json,resource,signal,subprocess,time
from functools import cmp_to_key
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import encode,decode
ROOT=Path(__file__).resolve().parents[1]

def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for z in iter(lambda:f.read(65536),b''):h.update(z)
    return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,indent=2)+'\n').encode()
def dump(p,v):p.write_bytes(enc(v))
def tree(p):return sum(x.stat().st_size for x in p.rglob('*')if x.is_file())
def numeric(v):
    if isinstance(v,bool):return 0
    if isinstance(v,(int,float)):return 1
    if isinstance(v,list):return sum(map(numeric,v))
    if isinstance(v,dict):return sum(map(numeric,v.values()))
    return 0
def fingerprint(p):
    h=hashlib.sha256();h.update((str(len(p['terms']))+'\n').encode());q=Q(p['content'])
    for t in p['terms']:
        c=q*int(t[0]);e=[p['shift'][j]+t[j+1]for j in range(5)];h.update((str(c)+' '+' '.join(map(str,e))+'\n').encode())
    return h.hexdigest()
def putliteral(f,p):
    f.write(str(Q(p['content']))+' '+' '.join(map(str,p['shift']))+' '+str(len(p['terms']))+' '+fingerprint(p)+'\n')
    for t in p['terms']:f.write(' '.join(map(str,t))+'\n')
def main(a):
    ID=a.id;GATE=a.gate
    begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out=ROOT/'reports'/ID;work=ROOT/'work'/ID
    if a.execute:
        assert a.prepared_sha and a.addendum_sha;execute(a,out,work,begin);return
    out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    gd=ROOT/'reports'/GATE;gr=gd/'report.json';assert sha(gr)==a.gate_sha;g=read(gr);assert g['all_passed']and g['FULL_safe_coefficient_job']
    bp=gd/'whole-budget.json';fp=gd/'sixteen-source-support-height-cell-plans.json'
    for p in(bp,fp):assert sha(p)==g['evidence_sha256'][p.name]
    b=read(bp);frame=read(fp);plans=frame['all16_cell_plans'];RAW=b['fresh_RAW_upper'];LIMIT=b['whole_owned_with_fixed_metadata_upper'];PAIRS=b['coefficient_term_pair_operations_upper'];selected=frame['selected_rows'];assert RAW<=32*1024*1024 and LIMIT<=150000 and len(selected)<=4
    prior=ROOT/frame['previous16_coefficient_cells_unchanged'];prior_report=read(prior);assert prior_report['all_passed'];prior_api=prior.parent/'evaluated16-cells-and-full177-frame.json';assert sha(prior_api)==prior_report['evidence_sha256'][prior_api.name];previous=read(prior_api);history=previous.get('all_evaluated_cells',previous.get('evaluated_cells16',[]));history=[{**z,'literal_cell_file':z.get('literal_cell_file',str((prior.parent/z['file']).relative_to(ROOT)))}for z in history];assert len(history)==4*len({z['row']for z in history}) and not set(selected)&{z['row']for z in history}
    active={ROOT/p:h for p,h in g['input_sha256'].items()};active[gr]=sha(gr)
    for p,h in g['evidence_sha256'].items():active[gd/p]=h
    for p,h in active.items():assert sha(p)==h
    commonp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_CACHED_KERNEL_COMMON_SHAPE_20261105_v64/fourteen-source-denominators-and-factored-P290.json'
    common=read(commonp);records=sorted(common['all14_fixed_denominator_equalities'],key=lambda z:(z['row_original0'],z['component']));keys={(z['row_original0'],z['component']):i for i,z in enumerate(records)}
    Ppath=commonp.parent/common['actual_common_polynomial'];P=read(Ppath);assert len(P['terms'])==290
    cpp=ROOT/'code/s61_sparse_generic5_constrained83_top_P290_aligned_streamcache_coefficients_20261105_v64.cpp';extra={Path(__file__),cpp,ROOT/f'research/notes/{ID}.md',ROOT/'code/s61_sparse_generic5_row_stream_kernel_20261104_v63.hpp',ROOT/'code/s61_sparse_generic5_primitive_binary_native_20261105_v64.hpp',ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py',ROOT/'code/s61_sparse_generic5_constrained83_aligned_group_storage_addendum_20261105_v64.py'};extra={p for p in extra if p not in active}
    base=sum(p.stat().st_size for p in active);external=sum(p.stat().st_size for p in extra);control=work/'control.txt';prep_peak=0
    with control.open('w')as f:
        matrix=ROOT/frame['literal_coefficient_source_file'];f.write(f"{base} {external} {LIMIT} {RAW} {PAIRS} {matrix} {frame['literal_coefficient_source_sha256']}\n")
        # Exactly one historical row payload is parsed at a time. Never parse
        # the whole old kernel or the current full matrix index.
        for row in sorted({z['row_original0']for z in records}):
            pair=[z for z in records if z['row_original0']==row];d=pair[0]['common_numerator_recipe']['source_row_payload'];p=ROOT/d['path']
            with p.open('rb')as src:src.seek(d['offset']);payload=src.read(d['bytes'])
            assert hashlib.sha256(payload).hexdigest()==d['payload_sha256'];v=json.loads(payload)
            for z in pair:
                rec=z['common_numerator_recipe'];a={'content':rec['QQ_scalar'],'shift':rec['parameter_exponent'],'terms':v[z['component']]['numerator']};putliteral(f,a)
                prep_peak=max(prep_peak,numeric(plans)+numeric(common)+numeric(P)+numeric(v)+numeric(frame['selected_rows'])+len(a['terms'])+8192);assert prep_peak<LIMIT
            del v,a,payload
        putliteral(f,P);f.write(str(len(selected))+'\n');support_index={}
        for row in selected:
            cp=[p for p in plans if p['row']==row];assert len(cp)==4;refs={z['source_column']:z['source_cell_reference']for p in cp for z in p['steps']}
            assert all(sum(z['terms']for z in refs.values())==p['source_row_input_terms']for p in cp)
            f.write(f'{row} {len(refs)}\n')
            for col,z in sorted(refs.items()):f.write(' '.join(map(str,[col,z['payload_offset'],z['payload_bytes'],z['terms'],z['payload_sha256']]))+'\n')
            f.write('4\n')
            for p in cp:
                assert p['column']in[7,13,27,33];f.write(' '.join(map(str,[p['column'],p['RAW_literal_polynomial_upper'],p['output_support_upper'],p['exact_integer_scalar_common_denominator'],p['integer_numerator_coefficient_abs_upper'],*p['actual_accumulator_shift'],len(p['steps'])]))+'\n')
                for z in p['steps']:
                    k=14 if z['kind']=='P_times_oldfree'else keys[tuple(z['cached_numerator_key'])]
                    f.write(' '.join(map(str,[z['source_column'],k,z['n_a'],z['n_b'],z['product_support_upper'],z['union_support_upper'],z['previous_union_support_upper'],z['whole_owned_scratch_upper']]))+'\n')
                sp=ROOT/p['support_file'];assert sha(sp)==p['support_sha256'];support=read(sp);assert len(support)==p['output_support_upper'];support.sort(key=cmp_to_key(lambda x,y:-1 if(sum(x),x)>(sum(y),y)else(1 if(sum(x),x)<(sum(y),y)else 0)));f.write(str(len(support))+'\n');offset=f.tell()
                for e in support:f.write(' '.join(map(str,e))+'\n')
                support_index[(row,p['column'])]={'offset':offset,'count':len(support),'bytes':f.tell()-offset,'source_sha256':p['support_sha256']}
                prep_peak=max(prep_peak,numeric(plans)+numeric(common)+numeric(P)+numeric(support)+len(support)+numeric(support_index)+8192);assert prep_peak<LIMIT
                del support
    caps=sum(p['RAW_literal_polynomial_upper']for p in plans);othercap=RAW-base-caps
    # Binary envelope and primitive codec are no larger than the SOURCE JSON
    # cap: header9; coefficient <=digits(M)+3 bytes; each exponent varint
    # <=decimal_digits(exponent)+1; source per-term48-byte slack/header2048.
    for p in plans:assert all(abs(v)<2**63 for v in p['actual_accumulator_shift'])
    future=200000+100000+20000+300000+200000
    pre={'id':ID,'all_passed':True,'selected_rows':selected,'only_changed_cells':4*len(selected),'source_gate_sha256':sha(gr),
        'whole_owned_limit':LIMIT,'fixed_native_metadata':25000,'whole_RAW_limit':RAW,'actual_READ_input_RAW_bytes':base,'worker_code_header_note_charge':external,
        'all16_binary_output_RAW_cap':caps,'other_artifact_RAW_cap':othercap,'future_other_upper':future,'thin_preparation_owned_peak':prep_peak,
        'streamed_support_offset_index':{str(k):v for k,v in support_index.items()},'prior_completed_layer_report':str(prior.relative_to(ROOT)),'prior_completed_layer_sha256':sha(prior),'binary_output_bound':'Q5CELL1 newline magic plus flag0 or flag1+canonical primitive-int-exponent-v1 payload; proven <=the literal JSON cap for every integer coefficient/shift/exponent in the SOURCE bounds.',
        'input_sha256':{**{str(p.relative_to(ROOT)):h for p,h in active.items()},**{str(p.relative_to(ROOT)):sha(p)for p in extra}},'control_sha256':sha(control),
        'never_parsed_full3872_index_or159249_kernel':True,'source_objects_released_before_native':True,'one_cached_numerator_native_policy':True,'fixed_common_base_direct_addition':True,'all_old_versions_preserved':True}
    dump(out/'preflight.json',pre);del plans,frame,common,records,P,b,g,active,previous;gc.collect();gc.collect();binary=work/'native'
    with(out/'compile.log').open('wb')as log:
        subprocess.run(['bash','-c','source /workspace/cloud-setup/activate.sh; exec c++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint openssl) -lgmp','compile',str(cpp),str(binary)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,59-(time.monotonic()-begin)))
    assert tree(out)+tree(work)+external+future<=othercap
    state={'all_passed':True,'source_only_preparation':True,'no_coefficient_worker_launched':True,'preflight_sha256':sha(out/'preflight.json'),'control_sha256':sha(control),'ignored_executable_sha256':sha(binary),'source_cell_plan':str(fp.relative_to(ROOT)),'coefficient_pair_limit':PAIRS,'compiled_worker_source':str(cpp.relative_to(ROOT)),'driver_source':str(Path(__file__).relative_to(ROOT)),'seconds':time.monotonic()-begin}
    dump(out/'prepared-state.json',state);signal.alarm(0);print(json.dumps(state));print('prepared_sha256',sha(out/'prepared-state.json'))

def execute(a,out,work,begin):
    assert not(out/'native.log').exists()and not(out/'native-summary.json').exists()and not list(out.glob('row-*.bin')),'inspect an ambiguous/in-flight worker; never blind rerun'
    pre=read(out/'preflight.json');statep=out/'prepared-state.json';assert sha(statep)==a.prepared_sha;state=read(statep)
    ap=out/'implementation-storage-addendum.json';assert sha(ap)==a.addendum_sha;addendum=read(ap);assert addendum['all_passed']and addendum['coefficient_job_safe']
    assert addendum['prepared_state_sha256']==a.prepared_sha and addendum['preflight_sha256']==sha(out/'preflight.json')
    ID=a.id;base=pre['actual_READ_input_RAW_bytes'];external=pre['worker_code_header_note_charge'];LIMIT=pre['whole_owned_limit'];RAW=pre['whole_RAW_limit'];PAIRS=state['coefficient_pair_limit'];selected=pre['selected_rows'];control=work/'control.txt';binary=work/'native'
    assert sha(binary)==state['ignored_executable_sha256']and sha(control)==pre['control_sha256']
    for p,h in pre['input_sha256'].items():assert sha(ROOT/p)==h
    assert tree(out)+tree(work)+external+addendum['remaining_noncell_artifact_upper']<=pre['other_artifact_RAW_cap']
    fp=ROOT/state['source_cell_plan'];support_index={tuple(map(int,k.strip('()').split(','))):v for k,v in pre['streamed_support_offset_index'].items()}
    prior=ROOT/pre['prior_completed_layer_report'];prior_report=read(prior);assert sha(prior)==pre['prior_completed_layer_sha256'];prior_api=prior.parent/'evaluated16-cells-and-full177-frame.json';assert sha(prior_api)==prior_report['evidence_sha256'][prior_api.name];previous=read(prior_api)
    history=previous.get('all_evaluated_cells',previous.get('evaluated_cells16',[]));history=[{**z,'literal_cell_file':z.get('literal_cell_file',str((prior.parent/z['file']).relative_to(ROOT)))}for z in history]
    assert not set(selected)&{z['row']for z in history};prep_peak=pre['thin_preparation_owned_peak'];del previous,state,addendum;gc.collect()
    try:
        with(out/'native.log').open('wb')as log:
            subprocess.run([str(binary),str(control),str(out),str(work)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,59-(time.monotonic()-begin)))
        n=read(out/'native-summary.json');assert n['all_passed']and n['cells']==4*len(selected) and n['rows']==len(selected) and n['coefficient_pairs']<=PAIRS
        records=[json.loads(s)for s in(out/'sixteen-evaluated-cell-records.jsonl').read_text().splitlines()];assert len(records)==4*len(selected)
        data=read(fp);plans={(p['row'],p['column']):p for p in data['all16_cell_plans']};del data;post_peak=0
        for z in records:
            p=out/z['file'];assert sha(p)==z['sha256'];blob=p.read_bytes();assert blob[:8]==b'Q5CELL1\n'and blob[8]in(0,1)
            if blob[8]:value=decode(blob[9:]);assert encode(value)==blob[9:]
            else:assert len(blob)==9;value={'content':'0','shift':[0]*5,'terms':[]}
            assert len(value['terms'])==z['terms']and fingerprint(value)==z['actual_QQ_value_sha256']
            pl=plans[(z['row'],z['column'])];sp=ROOT/pl['support_file'];assert sha(sp)==pl['support_sha256'];D=int(pl['exact_integer_scalar_common_denominator']);M=int(pl['integer_numerator_coefficient_abs_upper']);q=Q(value['content']);si=support_index[(z['row'],z['column'])];matched=0
            with control.open('rb')as stream:
                stream.seek(si['offset'])
                for _ in range(si['count']):
                    e=tuple(map(int,stream.readline().split()));assert len(e)==5
                    if matched==len(value['terms']):continue
                    term=value['terms'][matched];actual=tuple(value['shift'][j]+term[j+1]for j in range(5))
                    if e!=actual:assert(sum(e),e)>(sum(actual),actual);continue
                    c=q*int(term[0])*D;assert c.denominator==1 and abs(c.numerator)<=M;matched+=1
                assert stream.tell()==si['offset']+si['bytes']and matched==len(value['terms'])
            post_peak=max(post_peak,numeric(value)+len(value['terms'])+numeric(plans)+numeric(pre)+numeric(support_index)+8192);assert post_peak<LIMIT
            assert p.stat().st_size<=pl['RAW_literal_polynomial_upper']and z['all_source_product_identities']and z['binary_roundtrip']and z['height_support_verified']
        for p,h in pre['input_sha256'].items():assert sha(ROOT/p)==h
        assert sha(control)==pre['control_sha256']and not list(out.rglob('*.partial'))
        all_evaluated=history+[{**z,'literal_cell_file':str((out/z['file']).relative_to(ROOT))}for z in records];assert len({(z['row'],z['column'])for z in all_evaluated})==len(all_evaluated)
        layer={'all_passed':True,'all_evaluated_cells':all_evaluated,'previous_layer_report':str(prior.relative_to(ROOT)),'previous_layer_report_sha256':sha(prior),'global_column_scale_applied_once_not_per_batch':True,'base_full177x24_operator':'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FACTORED_CONSTRAINT_SCHUR_PREFLIGHT_20261105_v64/actual-factored177x24-and-top111x24.json',
            'only_column_scale':{'free_columns':[7,13,27,33],'factor':'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_CACHED_KERNEL_COMMON_SHAPE_20261105_v64/actual-common-denominator-P290.json'},
            'evaluated_cells_local':records,'remaining_top_row_count':111-len({z['row']for z in all_evaluated}),'uncomputed_free_columns':'P290 times corresponding frozen Schur root, once globally; prior evaluated cells are merged without another P factor.',
            'all20_other_columns_immutable_source_references':True,'all34_original_variable_substitutions_kept':True,'full177_rows_and66_formal_zero_proofs_kept':True,
            'cokernel':'full177x24 has coker(top111x24) direct_sum K5^66; no free zero-row summand was removed.',
            'chart':'P290*detA*detB and old parameter-monomial units nonzero; complementary divisors open. No full111 calculation/inverse/exhaustion claim.'}
        dump(out/'evaluated16-cells-and-full177-frame.json',layer)
        report={'id':ID,'all_passed':True,'selected_rows':selected,'changed_QQ_cells':4*len(selected),'cumulative_evaluated_cells':len(all_evaluated),'all_source_polynomial_identities':True,'all_binary_roundtrips':True,
            'source_back_add_checks':n['back_adds'],'coefficient_pairs':n['coefficient_pairs'],'actual_output_terms':sum(z['terms']for z in records),
            'whole_owned_peak':max(prep_peak,post_peak,n['whole_owned_peak'],n['prospective_whole_peak']),'no_full_index_kernel_parse':True,'remaining_top_rows_unevaluated':111-len({z['row']for z in all_evaluated}),
            'all20_columns_and66_zero_row_proofs_retained':True,'input_sha256':{**pre['input_sha256'],str(statep.relative_to(ROOT)):sha(statep),str(ap.relative_to(ROOT)):sha(ap)},'evidence_sha256':{str(p.relative_to(out)):sha(p)for p in out.rglob('*')if p.is_file()},
            'ignored_executable_sha256':sha(binary),'native_seconds':n['seconds'],'total_seconds':time.monotonic()-begin,'actual_fresh_RAW':0,'whole_RAW_upper':RAW}
        for _ in range(10):
            v=base+external+tree(out)+tree(work)+len(enc(report))
            if report['actual_fresh_RAW']==v:break
            report['actual_fresh_RAW']=v
        assert report['whole_owned_peak']<=LIMIT and report['actual_fresh_RAW']<=RAW and report['total_seconds']<60
        rp=out/'report.json';dump(rp,report);assert base+external+tree(out)+tree(work)==report['actual_fresh_RAW']
        print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(rp))
    except Exception as e:
        dump(out/'failure.json',{'all_passed':False,'failure':repr(e),'seconds':time.monotonic()-begin,'all_complete_cells_preserved':True,'no_automatic_retry':True});raise
    finally:signal.alarm(0)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--id',required=True);p.add_argument('--gate',required=True);p.add_argument('--gate-sha',required=True);p.add_argument('--execute',action='store_true');p.add_argument('--prepared-sha');p.add_argument('--addendum-sha');main(p.parse_args())
