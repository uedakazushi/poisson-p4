#!/usr/bin/env python3
"""Read-only exact source roots/atom inventories for all eight remaining tails."""
import argparse, hashlib, json, resource, signal, time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_PIVOT14_LASTROW_SOURCE_PROFILE_20261105_v63'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    out.mkdir(parents=True,exist_ok=False);pins={}
    def load(rel):
        p=ROOT/rel;pins[rel]=sha(p);return json.loads(p.read_text())
    qr='reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63/certificate.json'
    cr='reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json'
    pr='reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63/certificate.json'
    q=load(qr);cp=load(cr);pd=load(pr);nd=q['scalar_DAG']['nodes']
    p14='reports/S61_NONSQUARE_PREPOWER_K92_PIVOT14_SCHUR_20261105_v63/'
    assert load(p14+'report.json')['all_passed'];op=load(p14+'operator.json')
    assert op['source_remaining_current15']==[14] and op['target_remaining_original32']==[1,3,4,6,7,8,9,17]
    BP={(1,0,1):Q(50),(0,2,1):Q(18),(0,1,1):Q(-9),(0,0,1):Q(-2),(0,1,0):Q(60),(0,0,0):Q(-40)}
    def closure(roots):
        seen=set();todo=list(roots)
        while todo:
            i=todo.pop()
            if i in seen:continue
            seen.add(i);a=nd[i]
            if a[0] in ['a','s','m']:todo.extend(a[1:])
        return sorted(seen)
    atoms={}
    def atom(i):
        if i in atoms:return atoms[i]
        a=nd[i];kind=a[0]
        if kind in ['a','s','m','one','z']:return None
        if kind=='checkpoint':
            num,den=cp['matrix25x32'][a[1]][a[2]];ref=dict(path=cr,json_pointer=f'/matrix25x32/{a[1]}/{a[2]}')
        elif kind=='accepted_pair_second':
            num,den=pd['Lframe_dual6x2'][a[1]][1];ref=dict(path=pr,json_pointer=f'/Lframe_dual6x2/{a[1]}/1')
        else:
            assert kind=='new_unit_atom';num,den=q['scalar_DAG']['new_unit_atoms_num_den'][a[1]]
            ref=dict(path=qr,json_pointer=f'/scalar_DAG/new_unit_atoms_num_den/{a[1]}')
        dd={tuple(m):Q(c) for m,c in den if Q(c)};assert dd
        if len(dd)==1:
            m,c=next(iter(dd.items()));assert m[1:]==(0,0);unit=dict(kind='QQ_times_a_power',QQ=str(c),a_power=m[0])
        else:
            assert set(dd)==set(BP);k=dd[next(iter(BP))]/BP[next(iter(BP))]
            assert k and all(dd[m]==k*BP[m] for m in BP);unit=dict(kind='QQ_times_FULLU_B',QQ=str(k))
        vals=[(tuple(m),Q(c)) for m,c in num if Q(c)]
        height=max([max(abs(v.numerator).bit_length(),v.denominator.bit_length()) for m,v in vals]+[0])
        box=[[min(m[k] for m,v in vals),max(m[k] for m,v in vals)] for k in range(3)] if vals else None
        rec=dict(node=i,definition=a,source=ref,numerator_terms=len(vals),denominator_terms=len(dd),known_FULLU_denominator=unit,
                 literal_numerator_box=box,max_literal_QQ_bits=height)
        atoms[i]=rec;return rec
    def stats(ids):
        aa=[atom(i) for i in ids];aa=[v for v in aa if v is not None]
        return dict(nodes=len(ids),atoms=len(aa),numerator_terms=sum(v['numerator_terms'] for v in aa),
                    max_numerator_terms=max([v['numerator_terms'] for v in aa]+[0]),all_atom_denominators_known_FULLU=True)
    profiles=[];union=set()
    for t in op['target_remaining_original32']:
        j=q['new22_target_original32_columns'].index(t);roots=[r[j] for r in q['remaining15x22_scalar_roots']]
        allids=closure(roots);pairids=closure([roots[r] for r in [2,4,5,8]]);lastids=closure([roots[14]])
        union.update(allids)
        profiles.append(dict(original32_target=t,current22_column=j,all15_source_roots=roots,
                             pair4_roots=[roots[r] for r in [2,4,5,8]],last_G134_root=roots[14],
                             all15_closure=allids,pair4_closure=pairids,lastrow_closure=lastids,
                             all15=stats(allids),pair4=stats(pairids),lastrow=stats(lastids),
                             actual_final_entry_reference=dict(path=p14+'operator.json',node_index=op['exports']['residual1x8'],component=[0,op['target_remaining_original32'].index(t)]),
                             cleared_entry_reference=dict(path=p14+'operator.json',node_index=op['exports']['cleared_residual1x8'],component=[0,op['target_remaining_original32'].index(t)])))
    assert profiles[-1]['all15']==dict(nodes=159,atoms=64,numerator_terms=10406,max_numerator_terms=394,all_atom_denominators_known_FULLU=True)
    dump(out/'eight-tail-source-map.json',dict(profiles=profiles,source_DAG_reference=qr,source_union_nodes=sorted(union),
          atom_inventory=[atoms[i] for i in sorted(atoms)],union=stats(sorted(union)),
          fixed_rawpair_current15=[2,4,5,8],known_rawpair_rank1_columns=[22,11,13,25,15,19,14,18],
          all_eight_remaining_columns_OUTSIDE_known_rawpair_columns=True,
          warning='Atom counts are NOT expanded-output/lifetime bounds. No pair correction or Q9 tail may be discarded.'))
    code=Path(__file__);pins[str(code.relative_to(ROOT))]=sha(code)
    note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note)
    total=sum(p.stat().st_size for p in out.iterdir());assert total<1024*1024
    rep=dict(id=ID,all_passed=True,frozen=True,read_only_source_profile=True,no_new_coefficient_minor_rank_point_prime=True,
             cheapest_literal_source_target=17,target17_all15_source_terms=10406,
             no_scalar_job_or_expanded_worker_budget=True,all_pair_corrections_and_tail_effects_retained=True,
             total_new_RAW_before_report=total,input_sha256=pins,evidence_sha256={p.name:sha(p) for p in out.iterdir()},seconds=time.monotonic()-start)
    dump(out/'report.json',rep);signal.alarm(0);print(json.dumps({k:v for k,v in rep.items() if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())
