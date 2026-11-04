#!/usr/bin/env python3
"""Source-only nested dense-accumulator lifetime for the same fixed raw8 minor."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_NESTED_RAW8_SUPPORT_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
    with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def main(out):
    started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
    base=ROOT/'reports/S61_NONSQUARE_PREPOWER_NORMALIZED_RAW8_FUSED_PREFLIGHT_20261104_v63';old=read(base/'report.json');op=read(base/'operator-and-cost.json');proof=read(base/'fused-accumulation-proof.json');assert old['all_passed']and not old['resource_gate_passes']
    N=op['scalar_nodes'];heavy={x['root']for x in op['heavy2_roots']};terms=op['laplace15'];outroot=op['output'];den=lambda i:int(N[i]['QQ_denominator']);size=lambda i:N[i]['bound']['support']if N[i]['bound']else 0
    def lattice(bound):
        lo,hi=bound['degree'];b=bound['box'];return sum(max(0,min(b[2][1],hi-a-e)-max(b[2][0],lo-a-e)+1)for a in range(b[0][0],b[0][1]+1)for e in range(b[1][0],b[1][1]+1))
    def constant(i):
        n=N[i]
        if n['op']=='atom'and 'constant'in n['reference']:return int(n['reference']['constant'])
        return None
    def flatten(i,sign=1):
        n=N[i];a=n['args']
        if n['bound']is None:return []
        if n['op']=='mul'and(a[0]in heavy or a[1]in heavy):
            h,v=(a[0],a[1])if a[0]in heavy else(a[1],a[0]);return [dict(sign=sign,H2=h,P2=v)]
        if n['op']in['add','sub']:return flatten(a[0],sign)+flatten(a[1],sign if n['op']=='add'else-sign)
        if n['op']=='mul':
            for j in [0,1]:
                c=constant(a[j])
                if c in[-1,1]:return flatten(a[1-j],sign*c)
        raise AssertionError(('Unexpected bottom4 source grammar',i,n))
    inner=[];ordinary_roots=[];inner_events={};outer_events={};scalar_inner=0
    for k,t in enumerate(terms):
        b=t['bottom4_root'];top=t['top4_root'];items=flatten(b);assert len(items)==6 and {x['H2']for x in items}==heavy
        for x in items:
            assert den(b)%(den(x['H2'])*den(x['P2']))==0
            x['integer_QQ_scale']=str(den(b)//(den(x['H2'])*den(x['P2'])))
            at=max(x['H2'],x['P2']);inner_events.setdefault(at,[]).append(dict(bottom_index=k,**x))
            ordinary_roots +=[x['H2'],x['P2']];scalar_inner+=size(x['H2'])*size(x['P2'])
        assert den(outroot)%(den(top)*den(b))==0
        outer_events.setdefault(max(top,b),[]).append(dict(bottom_index=k,top4=top,bottom4=b,integer_QQ_scale=str(den(outroot)//(den(top)*den(b))),sign=-1 if sum(t['pure_top4_columns'])%2 else 1))
        ordinary_roots.append(top);inner.append(dict(bottom_index=k,source_root=b,bound=N[b]['bound'],dense_slots=lattice(N[b]['bound']),QQ_denominator=str(den(b)),six_terms=items))
    need=set();todo=ordinary_roots[:]
    while todo:
        i=todo.pop()
        if i in need:continue
        need.add(i);todo+=N[i]['args']
    # Bottom roots and their product/add intermediates are NOT evaluated as polynomials.
    assert all(t['bottom4_root']not in need for t in terms)
    uses={i:0 for i in need}
    for i in need:
        for j in N[i]['args']:uses[j]+=1
    for events in inner_events.values():
        for e in events:
            for j in [e['H2'],e['P2']]:uses[j]+=1
    for events in outer_events.values():
        for e in events:uses[e['top4']]+=1
    assert set(outer_events).issubset(need)
    dense_outer=lattice(N[outroot]['bound']);atom_reserve=sum(size(i)for i in need if N[i]['op']=='atom');live=peak=0;peakid=0;active={};records=[];completed=[]
    def consume(j):
        nonlocal live
        uses[j]-=1;assert uses[j]>=0
        if not uses[j]:live-=size(j)
    for i in sorted(need):
        n=N[i];live+=size(i);scratch=size(i)if n['op']in['mul','add','sub']else 0
        total=live+sum(active.values())+scratch+dense_outer+atom_reserve+15000
        if total>peak:peak=total;peakid=i
        for j in n['args']:consume(j)
        for e in inner_events.get(i,[]):
            k=e['bottom_index']
            if k not in active:
                assert not active,('Two bottom buffers would overlap',i,active,k)
                active[k]=inner[k]['dense_slots']
            now=live+sum(active.values())+dense_outer+atom_reserve+15000
            if now>peak:peak=now;peakid=i
            consume(e['H2']);consume(e['P2'])
        for e in outer_events.get(i,[]):
            k=e['bottom_index'];assert k in active;consume(e['top4']);del active[k];completed.append(k)
        records.append(dict(node=i,ordinary_live_upper=live,active_bottom_buffers=dict(active),output_outer_slots=dense_outer,inner_event_count=len(inner_events.get(i,[])),outer_event_count=len(outer_events.get(i,[]))))
    assert completed==list(range(15))and not active and live==0
    total_scalar=scalar_inner+proof['scalar_convolution_pairs_upper'];artifact=old['raw_output_and_atom_artifact_bytes_upper']
    plan=dict(inner15=inner,inner_events={str(i):v for i,v in inner_events.items()},outer_events={str(i):v for i,v in outer_events.items()},ordinary_nodes=sorted(need),all_normalized_source_nodes_external=dict(path=str((base/'operator-and-cost.json').relative_to(ROOT)),sha256=sha(base/'operator-and-cost.json')),
      outer_output_root=outroot,outer_dense_slots=dense_outer,outer_QQ_denominator=str(den(outroot)),
      last_use_records=records,one_bottom_buffer_at_a_time_proved=True,all90_inner_scales_and15_outer_scales_integer=True,
      exactly_same_15term_determinant_no_new_candidate=True,
      no_p_R_Saux_unknown_inverse=True,scalar_coefficients_or_determinants_NOT_evaluated=True)
    dump(out/'nested-plan.json',plan)
    pins={str(p.relative_to(ROOT)):sha(p)for p in [base/'report.json',base/'operator-and-cost.json',base/'fused-accumulation-proof.json',Path(__file__),ROOT/f'research/notes/{ID}.md']}
    result=dict(id=ID,all_passed=True,frozen=True,source_only=True,resource_gate_passes=peak<=150000 and artifact<32*1024*1024,
      total_live_upper=peak,peak_node=peakid,outer_dense_coefficient_slots=dense_outer,max_bottom_dense_slots=max(x['dense_slots']for x in inner),
      literal_atom_reserve=atom_reserve,unchanged_work_reserve=15000,ordinary_required_nodes=len(need),
      scalar_convolution_pairs_upper=total_scalar,raw_output_and_atom_artifact_bytes_upper=artifact,
      coefficient_integer_bits=old['output_integer_l1_bits'],QQ_denominator_bits=old['QQ_denominator_bits'],
      runtime_60s_NOT_guaranteed=True,no_new_scalar_job_or_value=True,input_sha256=pins,evidence_sha256={'nested-plan.json':sha(out/'nested-plan.json')},seconds=time.monotonic()-started)
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;dump(out/'report.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
