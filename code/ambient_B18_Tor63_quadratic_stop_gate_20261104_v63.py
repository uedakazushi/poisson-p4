#!/usr/bin/env python3
"""ONE cheap Tor2 quadratic-target test against intrinsic C0 Koszul only."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_T3_zero_structural_20261104_v63 import MON,POS
from ambient_B18_zero_normalizer_quadratic_20261104_v63 import compose
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_TOR63_QUADRATIC_STOP_GATE_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def enc(a):
    if isinstance(a,F):return str(a)
    if isinstance(a,(list,tuple)):return [enc(x) for x in a]
    if isinstance(a,dict):return {str(k):enc(v) for k,v in a.items()}
    return a
def subtract(v,b,c):
    for i,a in b.items():
        q=v.get(i,F(0))-c*a
        if q:v[i]=q
        else:v.pop(i,None)
def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    tp=ROOT/'reports/AMBIENT_B18_FINITE_TOR63_IMAGE_PREFLIGHT_20261104_v63/finite63-proof-and-input-preflight.json'
    ep=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json'
    fp=ROOT/'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json'
    jp=ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json'
    gates=[tp.parent/'report.json',ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63_REVIEW/report.json']
    assert all(read(p)['all_passed'] for p in gates)
    T,E,C,J=map(read,[tp,ep,fp,jp]);q0=[[F(v) for v in a] for a in T['all6_zero_column_cycle_origins16']]
    assert len(q0)==6 and all(len(q)==16 for q in q0)
    rows=T['actual_ideal_generator_physical_rows8']
    X=[tuple(map(F,a)) for a in E['families']['CA']['physical_R0_readback_jets_in_own_F9']]
    pairs={f:{r:[[F(v) for v in p] for p in pp] for r,pp in E['families'][f]['all12_trace_pairs_in_common_R0_through2']} for f in ['CA','CB']}
    delta=[]
    for row in rows:
        diff=tuple(a-b for a,b in zip(pairs['CA'][row][0],pairs['CB'][row][0]))
        assert not any(diff[:10])
        d=compose(diff,X);assert not any(d[:10]);delta.append(d)
    # Constant C0 free9 images of the eight actual CA odd target generators.
    fa=J['both57_source2']['CA'];f63=[F(v) for v in fa['pure_function2jets'][str(fa['f63'])]]
    assert not any(f63[8:10])
    Xi=[]
    for i in range(7):Xi.append([F(int(k==i+1)) for k in range(9)])
    Xi.append([F(0)]+f63[1:8]+[F(1)])
    # Rebuild the intrinsic28 first columns directly, then compare the frozen actual matrix.
    first=[];labels=[]
    for i,j in itertools.combinations(range(8),2):
        v=[F(0)]*81
        if j<7:v[9*(1+i)+j]=-1
        if i<7:v[9*(1+j)+i]=1
        first.append(v);labels.append([i,j])
    assert first==[[F(v) for v in c] for c in C['all640_exact_D0_firstjet81_columns'][:28]]
    columns=[];column_labels=[]
    for k,v in enumerate(first):
        for variable in range(9):
            w=[F(0)]*405
            for generator in range(9):
                for j in range(9):
                    a=v[9*generator+j]
                    if a:w[45*generator+POS[tuple(sorted((j,variable)))]-10]+=a
            columns.append(w);column_labels.append({'intrinsic_relation':k,'Koszul_indices':labels[k],'multiplier_F9_variable':variable})
    assert len(columns)==252
    basis=[];independent=[]
    for k,w in enumerate(columns):
        v={i:a for i,a in enumerate(w) if a};c={k:F(1)}
        for p,b,bc in basis:
            if p in v:
                a=v[p];subtract(v,b,a);subtract(c,bc,a)
        if not v:continue
        p=min(v);a=v[p]
        basis.append((p,{i:q/a for i,q in v.items()},{i:q/a for i,q in c.items()}));independent.append(k)
    targets=[]
    for alpha,beta in itertools.combinations(range(6),2):
        h=[[q0[alpha][i]*q0[beta][8+j]-q0[beta][i]*q0[alpha][8+j] for j in range(8)] for i in range(8)]
        vector=[F(0)]*405
        for i in range(8):
            for j in range(8):
                if not h[i][j]:continue
                for generator in range(9):
                    c=h[i][j]*Xi[i][generator]
                    if c:
                        for m,a in enumerate(delta[j][10:]):vector[45*generator+m]+=c*a
        residual={i:a for i,a in enumerate(vector) if a};coeff={}
        for p,b,bc in basis:
            if p in residual:
                a=residual[p];subtract(residual,b,a)
                # This is the coefficient of the ORIGINAL252 columns being subtracted.
                for k,c in bc.items():coeff[k]=coeff.get(k,F(0))+a*c
        coeff={k:c for k,c in coeff.items() if c}
        rebuilt=[sum((c*columns[k][i] for k,c in coeff.items()),F(0))+residual.get(i,F(0)) for i in range(405)]
        assert rebuilt==vector
        targets.append({'q_cycle_pair':[alpha,beta],'constant_h8_by8':h,
          'actual_C0_quadratic_target405':vector,'intrinsic252_membership_coefficients':sorted(coeff.items()),
          'exact_residual':sorted(residual.items()),'in_intrinsic_m_times_Koszul_first_span':not bool(residual)})
    allin=all(a['in_intrinsic_m_times_Koszul_first_span'] for a in targets)
    proof={'actual_CA_F9_order':['zero0..zero6','r','s'],'quadratic_monomials9':list(map(list,MON[10:])),
      'target_row_order':'C0free9 generator-major, then45 ordinary quadratic monomials',
      'actual_C0_basis_origin_images_of_CA8_odd_generators':Xi,
      'all8_even_differences_in_actual_CA_F9_through2':delta,
      'all28_intrinsic_C0_Koszul_first_columns81':first,
      'all252_intrinsic_maxideal_times_first_columns405':columns,'all252_column_provenance':column_labels,
      'intrinsic_quadratic_span_rank':len(basis),'independent_intrinsic_columns':independent,
      'all15_Tor2_quadratic_targets_and_reductions':targets,
      'all15_in_intrinsic_span':allin,
      'entire_finite63_image_is_in_m_cubed_D0_if_allin':allin,
      'no_actual_full_D0_second_relation_upper_bound_computed':True,
      'outside_intrinsic_span_would_not_alone_prove_actual_nonzero':True,
      'no_full_transform_or_new_source_jet_or_opaque_second_matrix':True}
    ep=out/'Tor63-quadratic-targets-and-intrinsic-stop.json';text=json.dumps(enc(proof),separators=(',',':'))+'\n'
    assert len(text.encode())<32*1024*1024;ep.write_text(text)
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',tp,
      ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json',fp,jp,*gates,
      ROOT/'code/ambient_B18_zero_normalizer_quadratic_20261104_v63.py']
    rep={'id':STEM,'all_passed':True,'status':'cheap_intrinsic_quadratic_stop_gate_only',
      'actual_finite_Tor2_candidate_targets':15,'target_QQ_rows':405,'intrinsic_relation_columns':252,
      'intrinsic_quadratic_span_rank':len(basis),'all15_targets_in_intrinsic_span':allin,
      'finite63_image_nonzero_or_zero_proved':False,'full_D_depth_or_socle_solved':False,
      'no_new_source_jets_higher_order_or_full640_second_carries':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{ep.name:sha(ep)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n')
    print(json.dumps({k:rep[k] for k in ['all_passed','intrinsic_quadratic_span_rank','all15_targets_in_intrinsic_span','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
