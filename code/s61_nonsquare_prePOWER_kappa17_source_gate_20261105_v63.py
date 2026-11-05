#!/usr/bin/env python3
"""Source-only support/lifetime/height gate for the fixed kappa dot v17."""
import argparse,functools,hashlib,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_KAPPA17_SOURCE_GATE_20261105_v63'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
@functools.lru_cache(None)
def lattice(box,deg):
    return sum(max(0,min(box[2][1],deg[1]-a-e)-max(box[2][0],deg[0]-a-e)+1)
               for a in range(box[0][0],box[0][1]+1) for e in range(box[1][0],box[1][1]+1))
def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    out.mkdir(parents=True,exist_ok=False);pins={}
    def load(rel):
        p=ROOT/rel;pins[rel]=sha(p);return json.loads(p.read_text())
    qr='reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63/certificate.json';q=load(qr)
    cr='reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json';cp=load(cr)
    pr='reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63/certificate.json';pd=load(pr)
    sr='reports/S61_NONSQUARE_PREPOWER_PIVOT14_LASTROW_SOURCE_PROFILE_20261105_v63/eight-tail-source-map.json';src=load(sr)
    assert load('reports/S61_NONSQUARE_PREPOWER_LASTROW17_CORRECTED_PAIR_SOURCE_20261105_v63/report.json')['all_passed']
    staged=load('reports/S61_NONSQUARE_PREPOWER_STAGED443_PREFLIGHT_20261104_v63/operator.json')
    assert staged['nodes'][22]['last_row_equals']=='minus kappa'
    assert load('reports/S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63_REVIEW/report.json')['all_passed']
    pp=next(p for p in src['profiles'] if p['original32_target']==17);ids=pp['pair4_closure'];roots=pp['pair4_roots']
    assert len(ids)==49 and pp['pair4']['numerator_terms']==2838
    atoms={a['node']:a for a in src['atom_inventory']};nd=q['scalar_DAG']['nodes']
    for i in ids:
        if i in atoms:assert atoms[i]['known_FULLU_denominator']['kind']=='QQ_times_a_power'
    N=[]
    def shape(v):
        if not v:return None
        D=math.lcm(*(c.denominator for m,c in v));L=sum(abs(c.numerator)*(D//c.denominator) for m,c in v)
        b=tuple((min(m[j] for m,c in v),max(m[j] for m,c in v)) for j in range(3));d=(min(sum(m) for m,c in v),max(sum(m) for m,c in v))
        return dict(support=len(v),box=b,degree=d,denominator=D,integer_l1=L)
    def combine(a,b,multiply=False):
        if multiply and (a is None or b is None):return None
        if a is None:return b
        if b is None:return a
        if multiply:
            box=tuple((x[0]+y[0],x[1]+y[1]) for x,y in zip(a['box'],b['box']));deg=(a['degree'][0]+b['degree'][0],a['degree'][1]+b['degree'][1])
            D=a['denominator']*b['denominator'];L=a['integer_l1']*b['integer_l1'];s=a['support']*b['support']
        else:
            box=tuple((min(x[0],y[0]),max(x[1],y[1])) for x,y in zip(a['box'],b['box']));deg=(min(a['degree'][0],b['degree'][0]),max(a['degree'][1],b['degree'][1]))
            D=math.lcm(a['denominator'],b['denominator']);L=a['integer_l1']*(D//a['denominator'])+b['integer_l1']*(D//b['denominator']);s=a['support']+b['support']
        return dict(support=min(s,lattice(box,deg)),box=box,degree=deg,denominator=D,integer_l1=L)
    bounds={};atomrefs=[];sourcepolycells=sum(len(n)+len(d) for row in cp['matrix25x32'] for n,d in row)
    assert sourcepolycells==88852
    for i in ids:
        node=nd[i];op=node[0]
        if op in ['a','s','m']:
            a,b=bounds[node[1]],bounds[node[2]]
            bounds[i]=None if op=='s' and node[1]==node[2] else combine(a,b,op=='m')
        elif op=='z':bounds[i]=None
        elif op=='one':bounds[i]=shape([((0,0,0),Q(1))])
        else:
            rec=atoms[i];ref=rec['source'];p=cp if ref['path']==cr else pd if ref['path']==pr else q
            v=p
            for k in ref['json_pointer'].strip('/').split('/'):v=v[int(k)] if isinstance(v,list) else v[k]
            num,den=v;assert len(den)==1;em,ec=den[0];assert em[1:]==[0,0] and Q(ec)
            vals=[((m[0]-em[0],m[1],m[2]),Q(c)/Q(ec)) for m,c in num if Q(c)]
            bounds[i]=shape(vals);atomrefs.append(dict(node=i,source=ref,source_reference_sha256=pins[ref['path']],
              QQ_unit_denominator=str(Q(ec)),known_a_denominator_power=em[0],literal_numerator_terms=len(vals)))
        N.append(dict(source_node=i,definition=node,bound=bounds[i]))
    assert all(j in bounds for j in roots)
    del cp,pd
    kapp=[];krefs=[]
    for j in [179,198,217,237]:
        rel=f'reports/S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63/root{j}.json';p=load(rel)['polynomial']
        kapp.append(shape([(tuple(m),Q(c)) for m,c in p if Q(c)]));krefs.append(dict(path=rel,sha256=pins[rel],json_pointer='/polynomial'))
    assert [x['support'] for x in kapp]==[1257,1257,1232,1174]
    # Source evaluation schedule with four output vectors retained. All input
    # coefficient arrays remain counted independently, even if a worker moves
    # instead of copying their data.
    users={i:0 for i in ids}
    for i in ids:
        if nd[i][0] in ['a','s','m']:
            for j in nd[i][1:]:users[j]+=1
    for i in roots:users[i]+=1
    inputcells=2838+sum(x['support'] for x in kapp);live=0;peak=0;schedule=[]
    def sz(i):return bounds[i]['support'] if bounds[i] else 0
    for i in ids:
        live+=sz(i);scratch=sz(i) if nd[i][0] in ['a','s','m'] else 0
        peak=max(peak,inputcells+live+scratch)
        schedule.append(dict(node=i,source_values_live=live,scratch=scratch,whole_live_input_included=inputcells+live+scratch))
        if nd[i][0] in ['a','s','m']:
            for j in nd[i][1:]:
                users[j]-=1
                if users[j]==0:live-=sz(j)
    outputv=sum(sz(i) for i in roots);assert live==outputv
    acc=None;products=[];dotsteps=[]
    for j in range(4):
        prod=combine(kapp[j],bounds[roots[j]],True);nacc=combine(acc,prod)
        # Old accumulator, new accumulator, and product buffers coexist;
        # source four-vector and all original inputs remain counted.
        owned=inputcells+outputv+(0 if acc is None else acc['support'])+(0 if prod is None else prod['support'])+(0 if nacc is None else nacc['support'])
        peak=max(peak,owned);products.append(prod);dotsteps.append(dict(term=j,product=prod,updated_accumulator=nacc,whole_live=owned));acc=nacc
    maxbits=max([x['integer_l1'].bit_length() for x in bounds.values() if x]+[x['integer_l1'].bit_length() for x in products if x]+[acc['integer_l1'].bit_length()])
    maxden=max([x['denominator'].bit_length() for x in bounds.values() if x]+[x['denominator'].bit_length() for x in products if x]+[acc['denominator'].bit_length()])
    # Canonical QQ pair coefficients: numerator bounded by the cleared L1,
    # denominator divides the tracked D. Count both decimal strings plus all
    # exponent/sign/index overhead; five published polynomials and sourceplan.
    publication_cells=outputv+acc['support'];whole_owned=peak+publication_cells+15000
    perterm=2*((max(maxbits,maxden)*30103+99999)//100000+2)+72
    newRAW=(publication_cells+inputcells)*perterm+2*1024*1024
    readRAW=sum((ROOT/rel).stat().st_size for rel in pins)
    rawgate=readRAW+newRAW<32*1024*1024;slotgate=whole_owned<=150000
    gate=dict(source49_nodes=49,input_source_atom_terms=2838,known_kappa_terms=4920,held_input_cells=inputcells,
      source_four_output_supports=[sz(i) for i in roots],cleared_kappa17_output=acc,
      output_a_exponent_min=acc['box'][0][0],no_B_clearing_needed=True,
      source_and_product_peak_slots=peak,publication_copy_slots=publication_cells,ALL_metadata_and_scratch_reserve=15000,
      whole_worker_owned_upper=whole_owned,tracked_integer_coefficient_bits_upper=maxbits,tracked_denominator_bits_upper=maxden,
      whole_publication_RAW_upper=newRAW,external_read_RAW=readRAW,TOTAL_RAW_upper=readRAW+newRAW,
      source_preparation_CPx25x32_coefficient_cells=sourcepolycells,source_prepare_cp_plus_retained_inputs_upper=sourcepolycells+inputcells+15000,
      slot_gate_passes=slotgate,RAW_gate_passes=rawgate,
      runtime_guard_seconds=60,memory_bytes=2000000000,scalar_job_RUN=False,
      scalar_conclusion_if_nonzero='only x4 numerator nonzero, not final r17 or global57',
      scalar_conclusion_if_zero='only x4=0 for target17, not raw proportionality or other corrections')
    dump(out/'source-support-and-worker-gate.json',gate)
    dump(out/'literal-source-plan.json',dict(variable_order=['a','e','d'],actual_current15_pair_rows=[2,4,5,8],actual_original32_target=17,
      source_DAG_path=qr,source_nodes=N,source_root_order=roots,source_atom_references=atomrefs,
      kappa_coefficients_in_exact_order=krefs,formula='kappa17=sum(kappa[i]*v17[i],i=0..3)',
      A4_solve_identity='x4=kappa17/T; no new adjugate coefficient, signs already present in literal kappa',
      source_evaluation_lastuse=schedule,dot_product_schedule=dotsteps,
      no_output_or_product_coefficients_evaluated=True))
    code=Path(__file__);pins[str(code.relative_to(ROOT))]=sha(code);note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note)
    size=sum(p.stat().st_size for p in out.iterdir());assert size<2*1024*1024
    rep=dict(id=ID,all_passed=True,frozen=True,source_only=True,no_new_scalar_coefficients=True,
      proposed_fixed_target='kappa17 fourth corrected-pair solve coordinate numerator',
      prospective_slot_gate_passes=slotgate,prospective_RAW_gate_passes=rawgate,
      whole_worker_owned_upper=whole_owned,TOTAL_RAW_upper=readRAW+newRAW,
      scalar_job_NOT_run=True,requires_separate_parent_GO=True,total_new_RAW_before_report=size,
      input_sha256=pins,evidence_sha256={p.name:sha(p) for p in out.iterdir()},seconds=time.monotonic()-start)
    dump(out/'report.json',rep);signal.alarm(0);print(json.dumps({k:v for k,v in rep.items() if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
