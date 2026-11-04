#!/usr/bin/env python3
"""ONE cheap m*full-first38 test of the saved15 Tor2 quadratic targets."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_T3_zero_structural_20261104_v63 import MON,POS
from ambient_B18_Tor63_quadratic_stop_gate_20261104_v63 import read,sha,enc,subtract
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_TOR63_FULL_FIRST_MULTIPLE_STOP_20261104_v63'
def basis_with_certificates(columns):
    b=[];ind=[]
    for k,w in enumerate(columns):
        v={i:a for i,a in enumerate(w) if a};c={k:F(1)}
        for p,bv,bc in b:
            if p in v:
                a=v[p];subtract(v,bv,a);subtract(c,bc,a)
        if not v:continue
        p=min(v);a=v[p]
        b.append((p,{i:q/a for i,q in v.items()},{i:q/a for i,q in c.items()}));ind.append(k)
    return b,ind
def reduce_with_certificate(vector,basis):
    v={i:a for i,a in enumerate(vector) if a};c={}
    for p,b,bc in basis:
        if p in v:
            a=v[p];subtract(v,b,a)
            for k,x in bc.items():c[k]=c.get(k,F(0))+a*x
    return v,{k:x for k,x in c.items() if x}
def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    tp=ROOT/'reports/AMBIENT_B18_TOR63_QUADRATIC_STOP_GATE_20261104_v63/Tor63-quadratic-targets-and-intrinsic-stop.json'
    fp=ROOT/'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json'
    gates=[tp.parent/'report.json',fp.parent/'report.json',
      ROOT/'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63_REVIEW/report.json']
    assert all(read(p)['all_passed'] for p in gates)
    T,C=map(read,[tp,fp]);first=[[F(v) for v in col] for col in C['all640_exact_D0_firstjet81_columns']]
    assert len(first)==640 and all(len(v)==81 for v in first)
    first_basis,selected=basis_with_certificates(first);assert len(first_basis)==38
    all640expressions=[]
    for col in first:
        rem,co=reduce_with_certificate(col,first_basis);assert not rem
        assert [sum((v*first[k][i] for k,v in co.items()),F(0)) for i in range(81)]==col
        all640expressions.append(sorted(co.items()))
    multiples=[];labels=[]
    for original in selected:
        col=first[original]
        for k in range(9):
            w=[F(0)]*405
            for generator in range(9):
                for j in range(9):
                    a=col[9*generator+j]
                    if a:w[45*generator+POS[tuple(sorted((j,k)))]-10]+=a
            multiples.append(w);labels.append({'actual_D0_original_relation_index':original,'F9_multiplier_index':k})
    assert len(multiples)==342
    basis,ind=basis_with_certificates(multiples);results=[]
    for old in T['all15_Tor2_quadratic_targets_and_reductions']:
        target=list(map(F,old['actual_C0_quadratic_target405']))
        residual,coeff=reduce_with_certificate(target,basis)
        assert [sum((c*multiples[k][i] for k,c in coeff.items()),F(0))+residual.get(i,F(0)) for i in range(405)]==target
        results.append({'q_cycle_pair':old['q_cycle_pair'],'saved_exact_quadratic_target405':target,
          'full_first342_membership_coefficients':sorted(coeff.items()),'exact_residual':sorted(residual.items()),
          'in_maxideal_times_full_first_relation_span':not bool(residual)})
    allin=all(r['in_maxideal_times_full_first_relation_span'] for r in results)
    proof={'actual_CA_F9_order':['zero0..zero6','r','s'],'quadratic_monomials9':list(map(list,MON[10:])),
      'all640_actual_D0_first_columns81':first,'first_span_rank':38,'independent38_original_relation_indices':selected,
      'all640_first_span_membership_coefficients':all640expressions,
      'all342_maxideal_times_independent38_columns405':multiples,'all342_column_provenance':labels,
      'multiple_span_rank':len(basis),'independent_multiple_columns':ind,
      'all15_saved_Tor2_targets_and_exact_reductions':results,'all15_in_multiple_span':allin,
      'finite63_image_contained_in_m_cubed_D0_if_allin':allin,
      'quadratic_route_stops_as_inconclusive_if_allin':allin,
      'second_carries_unknown_sensitivities_new_jets_or_literal_transform_used':False,
      'actual_finite_image_zero_nonzero_socle_or_depth_claimed':False}
    ep=out/'full-first-multiples-and-Tor63-targets.json';text=json.dumps(enc(proof),separators=(',',':'))+'\n'
    assert len(text.encode())<32*1024*1024;ep.write_text(text)
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',tp,fp,*gates,
      ROOT/'code/ambient_B18_Tor63_quadratic_stop_gate_20261104_v63.py',
      ROOT/'code/ambient_B18_T3_zero_structural_20261104_v63.py']
    rep={'id':STEM,'all_passed':True,'status':'cheap_full_first_multiple_stop_gate_only',
      'actual_first_span_rank':38,'multiple_columns':342,'quadratic_QQ_rows':405,
      'multiple_span_rank':len(basis),'saved_finite_Tor2_targets':15,'all15_in_multiple_span':allin,
      'finite63_image_zero_nonzero_or_full_D_depth_decided':False,
      'no_second_carries_opaque_sensitivities_new_jets_or_higher_job':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{ep.name:sha(ep)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps({k:rep[k] for k in ['all_passed','multiple_span_rank','all15_in_multiple_span','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
