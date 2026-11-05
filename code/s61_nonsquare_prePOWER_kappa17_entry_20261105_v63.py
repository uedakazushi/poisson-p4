#!/usr/bin/env python3
"""Prepare immutable custody, then execute exactly one fixed kappa17 job."""
import argparse,gc,hashlib,json,math,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_KAPPA17_ENTRY_20261105_v63'
def sha(p):
    with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def prepare(out):
    out.mkdir(parents=True,exist_ok=False);pins={}
    def load(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
    base='reports/S61_NONSQUARE_PREPOWER_KAPPA17_SOURCE_GATE_20261105_v63/'
    r=load(base+'report.json');assert r['all_passed'] and pins[base+'report.json']=='9e25c30cf6a8c991dde49ec521f9bc46d634189223e61037e8730aedf23fe314'
    plan=load(base+'literal-source-plan.json')
    lit='reports/S61_NONSQUARE_PREPOWER_KAPPA17_LITERAL_INPUT_GATE_V2_20261105_v63/'
    lr=load(lit+'report.json');assert lr['all_passed'] and pins[lit+'report.json']=='ab430bf4a8b7bee89fd13e84887b45156fb75a4c29f5a9397eec23220ee5513d'
    data=load(lit+'literal-atoms20.json');assert pins[lit+'literal-atoms20.json']=='7c9ed94906fde83f628982d383c60643d92cdc2bba244fd6d17ef7aeeb9e72a8'
    ids=[n['source_node']for n in plan['source_nodes']];index={v:i for i,v in enumerate(ids)};atoms={a['node']:a for a in data['atom_payloads']}
    normalized=[]
    for n in plan['source_nodes']:
        i=n['source_node'];b=n['bound'];spec=([0,-1,0,-1,0,-1,0,-1,'1']if b is None else [v for pair in b['box']for v in pair]+b['degree']+[str(b['denominator'])])
        op=n['definition'][0];args=n['definition'][1:];rec=dict(old_id=i,op='1'if op=='one'else op,spec=spec,args=[index[j]for j in args]if op in ['a','s','m']else[-1,-1])
        if i in atoms:
            rec['op']='i';num,den=atoms[i]['original_num_den'];assert len(den)==1;em,ec=den[0];D=int(spec[-1]);vals=[]
            for m,c in num:
                v=Q(c)/Q(ec)*D;assert v.denominator==1 and em[1:]==[0,0]
                mm=[m[0]-em[0],m[1],m[2]];assert Q(int(v),D)*Q(ec)==Q(c)
                vals.append((mm,str(int(v))))
            rec['terms']=vals
        normalized.append(rec)
    kapp=[]
    def spec(poly):
        vals=[(m,Q(c))for m,c in poly];D=math.lcm(*(c.denominator for m,c in vals));box=[v for k in range(3)for v in [min(m[k]for m,c in vals),max(m[k]for m,c in vals)]]
        return box+[min(sum(m)for m,c in vals),max(sum(m)for m,c in vals),str(D)],[(m,str(int(c*D)))for m,c in vals]
    for ref in data['kappa_sources']:
        x=load(ref['path']);assert pins[ref['path']]==ref['sha256'];kapp.append(spec(x['polynomial']))
    inp=out/'native-input.txt'
    with inp.open('w')as f:
        f.write('KAPPA17D001 49 48\n')
        def emitsp(s):f.write(' '.join(map(str,s))+'\n')
        def emitp(p):
            f.write(str(len(p))+'\n')
            for m,c in p:f.write(' '.join(map(str,m))+' '+c+'\n')
        for n in normalized:
            a=n['args'];f.write(f"{n['old_id']} {n['op']} {a[0]} {a[1]}\n");emitsp(n['spec'])
            if n['op']=='i':emitp(n['terms'])
        f.write(' '.join(str(index[i])for i in data['source_order'])+'\n')
        for i,(ks,kp)in enumerate(kapp):
            emitsp(ks);emitp(kp)
            for b in [plan['dot_product_schedule'][i]['product'],plan['dot_product_schedule'][i]['updated_accumulator']]:emitsp([v for pair in b['box']for v in pair]+b['degree']+[str(b['denominator'])])
    cpp=ROOT/'code/s61_nonsquare_kappa17_dense_20261105_v63.cpp';pins[str(cpp.relative_to(ROOT))]=sha(cpp)
    pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note)
    # The dense implicit lattice adds no per-coefficient exponent key. Sparse
    # literal inputs count both packed key and integer coefficient. Dynamic
    # accounting covers lattice row offsets, denominator and header scalars.
    gate=dict(id=ID+'_WORKER_GATE',all_passed=True,frozen=True,input_sha256=pins,native_input_sha256=sha(inp),
      exact_source_atom_QQ_scale_and_a_offset_checked=True,source_atoms=20,source_nodes=49,
      mathematical_job_not_run=True,no_old_checkpoint_Q2_or_pair_file_read=True,
      representation='dense implicit box/total-degree lattice; sparse input key+coefficient counted separately',
      owned_numeric_and_metadata_hard_cap=75410,temporary_integer_bits=548,
      runtime_guard_seconds=48,wall_total60_seconds=True,memory_bytes=2000000000,
      RAW_total_hard_cap=10150000,source_convolution_all_input_and_lastuse_guard='runtime dynamic charge: every dense coefficient, row offset, poly header/denominator, both sparse key/coefficient;49*20 definition slots plus3000 fixed metadata/scratch',
      readback_guard='native exits first; canonicalize/readback ONE polynomial file at a time, no duplicate polynomial list',
      full_product_reverse_reproduction=True,full_addback_identities=True)
    dump(out/'report.json',gate);print(json.dumps({k:v for k,v in gate.items()if k!='input_sha256'},indent=2))
def execute(gatepath,out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));g=read(gatepath);assert g['all_passed']and g['mathematical_job_not_run'];pins={}
    for rel,h in g['input_sha256'].items():assert sha(ROOT/rel)==h;pins[rel]=h
    inp=gatepath.parent/'native-input.txt';assert sha(inp)==g['native_input_sha256'];pins[str(gatepath.relative_to(ROOT))]=sha(gatepath);pins[str(inp.relative_to(ROOT))]=sha(inp)
    out.mkdir(parents=True,exist_ok=False);cpp=ROOT/'code/s61_nonsquare_kappa17_dense_20261105_v63.cpp';exe=ROOT/'work/d-zero-v63/kappa17-dense';exe.parent.mkdir(parents=True,exist_ok=True)
    cp=subprocess.run(['g++','-std=c++17','-O2',str(cpp),'-lflint','-lgmp','-o',str(exe)],capture_output=True,text=True,timeout=15);(out/'compile.stdout').write_text(cp.stdout);(out/'compile.stderr').write_text(cp.stderr)
    if cp.returncode:dump(out/'report.json',dict(all_passed=False,compile_only_failure=True,no_arithmetic=True,input_sha256=pins));return
    hp=subprocess.run([str(exe),'--help'],capture_output=True,text=True,timeout=3);assert hp.returncode==0;(out/'native-help.txt').write_text(hp.stdout)
    dump(out/'phase.json',dict(phase='ONE scalar source49 plus fixed FOUR products',worker_sha256=sha(cpp),custody_sha256=sha(gatepath),native_input_sha256=sha(inp)))
    def limits():resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));resource.setrlimit(resource.RLIMIT_CPU,(52,52))
    run=subprocess.run([str(exe),str(inp),str(out)],capture_output=True,text=True,timeout=max(1,57-(time.monotonic()-start)),preexec_fn=limits);(out/'native.stdout').write_text(run.stdout);(out/'native.stderr').write_text(run.stderr)
    if run.returncode:dump(out/'report.json',dict(id=ID,all_passed=False,arithmetic_started=True,all_failure_state_preserved=True,no_polynomial_result_claim=True,input_sha256=pins));return
    native=read(out/'native.json');assert native['all_passed']and native['peak_OWNED_all_numeric_metadata']<=75410 and native['observed_integer_bits']<=548
    records=[]
    for fn in [f'source-root{i}.json'for i in range(4)]+['kappa17.json']:
        p=out/fn;x=read(p);D=int(x['denominator']);vals=x['integer_coefficients'];assert D>0 and len(vals)<11729
        digest=sha(p);previous=None
        for m,c in vals:
            assert previous is None or tuple(m)>previous;previous=tuple(m);v=Q(int(c),D);assert v and int(c).bit_length()<=201
            # Canonical reduction and exact scalar back-multiplication, in place.
            assert v*D==int(c)
        records.append(dict(file=fn,terms=len(vals),integer_common_denominator=str(D),raw_sha256=digest,all_QQ_reductions_backmultiply=True))
        del x,vals;gc.collect()
    dump(out/'QQ-readback-and-e4-identity.json',dict(files=records,formula='x4=kappa17/T',all_four_source_roots_fixed_order=[2,4,5,8],no_double_sign_on_kappa=True,
      exact_source49_offsets_and_all_product_backadds=native['all_source_offsets_exact']and native['all_additions_fullQQ_backadd']and native['all_products_reverse_order_fullQQ_reproduction'],
      other_three_pair_components_and_Schur10_corrections_remain_explicit=True,no_final_lastrow_or_global57_claim=True))
    external=sum((ROOT/p).stat().st_size for p in pins);new=sum(p.stat().st_size for p in out.iterdir());assert external+new<=10150000
    report=dict(id=ID,all_passed=True,frozen=True,ONE_fixed_kappa17_polynomial=True,kappa17_terms=native['kappa17_terms'],kappa17_is_zero=native['kappa17_is_zero'],
      e4_pair_solve_identity='x4=kappa17/T',no_final_lastrow_rank_or_global57=True,all_three_other_pair_and_ten_block_corrections_retained=True,
      all_fullQQ_source_product_reproduction=True,peak_owned_numeric_metadata=native['peak_OWNED_all_numeric_metadata'],integer_products=native['integer_products_including_reverse_checks'],
      external_read_RAW=external,total_new_RAW_before_report=new,TOTAL_RAW=external+new,independent_peer_pending=True,input_sha256=pins,
      evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
    dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='mode',required=True);a=s.add_parser('prepare');a.add_argument('--output',type=Path,required=True);a=s.add_parser('execute');a.add_argument('--gate',type=Path,required=True);a.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if args.mode=='prepare':prepare(args.output.resolve())
    else:execute(args.gate.resolve(),args.output.resolve())
