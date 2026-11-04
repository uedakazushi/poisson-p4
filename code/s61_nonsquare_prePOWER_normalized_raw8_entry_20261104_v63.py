#!/usr/bin/env python3
"""ONE actual normalized raw8 coefficient job under the frozen nested fusion guard."""
import argparse,hashlib,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_NORMALIZED_RAW8_ENTRY_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
    with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def main(out):
    began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);assert out.is_dir()and not(out/'report.json').exists();phase=read(out/'phase.json');pins={};cache={}
    def source(rel):
        p=ROOT/rel
        if rel not in cache:cache[rel]=read(p);pins[rel]=sha(p)
        return cache[rel]
    gatebase='reports/S61_NONSQUARE_PREPOWER_NESTED_RAW8_SUPPORT_PREFLIGHT_20261104_v63';gate=source(gatebase+'/report.json');plan=source(gatebase+'/nested-plan.json')
    assert sha(ROOT/(gatebase+'/report.json'))=='a3f480f9f1904cade9ffc42af2e746ad07c19dd72cfb3bed669f6c40eeaec230'and gate['resource_gate_passes']and gate['total_live_upper']==131590
    op=source(plan['all_normalized_source_nodes_external']['path']);N=op['scalar_nodes'];assert len(N)==2094;need=set(plan['ordinary_nodes']);assert len(need)==1898
    q=source('reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63/certificate.json');cp=source('reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json');pd=source('reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63/certificate.json')
    def atom(n):
        r=n['reference']
        if 'constant'in r:return [[[0,0,0],r['constant']]]
        if r['path']=='accepted_Q2L_source_atom':
            a=r['definition'];num,den=cp['matrix25x32'][a[1]][a[2]]if a[0]=='checkpoint'else pd['Lframe_dual6x2'][a[1]][1]if a[0]=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][a[1]]
            assert len(den)==1 and den[0][0][1:]==[0,0];return num
        d=source(r['path'])
        if 'old_DAG_root'in r:return next(x['Laurent_polynomial']for x in d['node_values']if x['old_DAG_root']==r['old_DAG_root'])
        if 'h9_target'in r:return d['source_h9'][d['original_target_labels'].index(r['h9_target'])]
        x=d[r['key']]
        for j in r.get('index',[]):x=x[j]
        return x
    atoms=[];ip=out/'native-input.txt'
    with ip.open('w')as f:
        f.write(f"{gate['literal_atom_reserve']} {len(N)}\n")
        for i,n in enumerate(N):
            if i not in need:f.write('x\n');continue
            prefix=f"{n['bound']['support']if n['bound']else 0} {n['QQ_denominator']}"
            if n['op']=='atom':
                p=atom(n);m=min((e[0]for e,c in p),default=0);shift=n['a_shift']+m;f.write(f"p {prefix} {shift} {n['QQ_unit']}\n{len(p)}\n")
                for e,c in p:f.write(f'{e[0]-m} {e[1]} {e[2]} {c}\n')
                atoms.append(dict(node=i,reference=n['reference'],terms=len(p),QQ_unit=n['QQ_unit'],a_shift=shift))
            elif n['op']=='known_a_shift':f.write(f"k {prefix} {n['args'][0]} {n['a_exponent_shift']}\n")
            else:f.write(dict(add='a',sub='s',mul='m')[n['op']]+' '+prefix+' '+' '.join(map(str,n['args']))+'\n')
        f.write('15\n')
        def spec(bound,den):return ' '.join(map(str,[v for p in bound['box']for v in p]+bound['degree']+[den]))
        f.write(spec(N[op['output']]['bound'],plan['outer_QQ_denominator'])+'\n')
        for x in plan['inner15']:f.write(spec(x['bound'],x['QQ_denominator'])+'\n')
        events=[]
        for at,ev in plan['inner_events'].items():
            for x in ev:events.append([int(at),0,x['bottom_index'],x['H2'],x['P2'],x['sign'],x['integer_QQ_scale']])
        for at,ev in plan['outer_events'].items():
            for x in ev:events.append([int(at),1,x['bottom_index'],x['top4'],x['bottom4'],x['sign'],x['integer_QQ_scale']])
        events.sort(key=lambda x:(x[0],x[1]));assert len(events)==105;f.write('105\n')
        for e in events:f.write(' '.join(map(str,e))+'\n')
    identity=dict(actual_input_current15_rows=op['source_current15_rows'],actual_input_original32_columns=op['source_original32_columns'],
      raw_det_equals_normalized_times='(a*F*C^2)^2*gamma*a^(sum row_known_a_factors)',row_known_a_factors=op['row_known_a_factors'],
      normalized_heavy_columns=op['normalized_heavy_columns_divided_by'],known_FULLU_units_ONLY=True,
      all_source_atoms=atoms,all1898_ordinary_roots_pinned=sorted(need),all90_inner_and15_outer_event_definitions=events,
      future_pivot11=op['future_pivot11'],no_p_R_Saux_h22_inverse=True,all8_Q9_other_tails_in_future_complement=True,
      completion_criterion='All105 event contributions and whole exact integer coefficient array complete; partialprefix alone is NOT determinant.')
    dump(out/'source-binding.json',identity);dump(out/'preflight.json',dict(frozen_gate_SHA256=pins[gatebase+'/report.json'],implementation_total_live_upper=131590,all_preloaded_atoms_and15000_work_counted=True,ordinary_node_count=1898,scalar_convolution_pairs_upper=gate['scalar_convolution_pairs_upper'],raw_artifact_bytes_upper=gate['raw_output_and_atom_artifact_bytes_upper'],runtime_uncertain_hard_guard=True))
    phase.update(phase='all literal inputs/events bound; before one native job');dump(out/'phase.json',phase)
    cache.clear();del q,cp,pd,N,op,atoms
    cpp=Path(__file__).with_suffix('.cpp');exe=out/'native-worker';cmd='source /workspace/cloud-setup/activate.sh && c++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'
    run=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,57-(time.monotonic()-began)));(out/'compile.log').write_text(run.stdout);assert run.returncode==0
    cmd='source /workspace/cloud-setup/activate.sh && exec '+str(exe)+' '+str(ip)+' '+str(out)
    run=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,58-(time.monotonic()-began)));(out/'native.log').write_text(run.stdout)
    p=out/'native.json'if(out/'native.json').exists()else out/'frontier.json';n=read(p)if p.exists()else dict(all_passed=False,returncode=run.returncode);ok=run.returncode==0 and n['all_passed']
    phase.update(phase='native complete'if ok else'one native stopped; exactpartial only; no retry');dump(out/'phase.json',phase)
    for pth in [Path(__file__),cpp,ROOT/'code/s61_nonsquare_prePOWER_h22_pure3_fixed_minor_20261104_v63.cpp',ROOT/f'research/notes/{ID}.md']:pins[str(pth.relative_to(ROOT))]=sha(pth)
    size=sum(p.stat().st_size for p in out.iterdir()if p.is_file()and p!=exe);assert size<32*1024*1024
    result=dict(id=ID,all_passed=ok,frozen=True,independent_review_pending=True,normalized_minor_is_zero=n.get('normalized_raw8_is_zero'),normalized_minor_terms=n.get('terms'),
      max_live_coefficient_slots=n.get('max_total_live_coefficient_slots',n.get('max_slots')),actual_scalar_updates=n.get('actual_scalar_addmul_updates',n.get('scalar_updates')),
      completed_inner=n.get('completed_inner90',n.get('completed_inner')),completed_outer=n.get('completed_outer15',n.get('completed_outer')),
      whole_minor_claim_ONLY_if_complete=ok,no_global57_or_remaining4rank_claim=True,
      raw_artifact_bytes_before_report_excluding_binary=size,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},seconds=time.monotonic()-began)
    dump(out/'report.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
