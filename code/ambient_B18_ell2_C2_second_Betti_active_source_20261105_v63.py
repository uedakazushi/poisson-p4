#!/usr/bin/env python3
"""Active-five-variable d3 source split, then one guarded known-matrix rank."""
import argparse,hashlib,itertools,json,resource,signal,time,traceback
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
ROOT=Path(__file__).resolve().parents[1]
PREF='AMBIENT_B18_ELL2_C2_SECOND_BETTI_ACTIVE_PREFLIGHT_20261105_v63'
STEM='AMBIENT_B18_ELL2_C2_SECOND_BETTI_ACTIVE_RANK_20261105_v63'
ACTIVE=[0,1,4,5,7];ZERO=[2,3,6]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(v):return {i:F(c) for i,c in v}
def blocks(n,variables):return [(n-q,w) for q in range(min(n,len(variables))+1) for w in itertools.combinations(variables,q)]
def differential(n,variables,A):
    b=blocks(n,variables);t=blocks(n-1,variables);idx={w:i for i,(p,w) in enumerate(t)};out=[]
    for p,w in b:
        for j in range(14):
            v={}
            if p:v=add(v,{14*idx[w]+i:(-1)**p*c for i,c in A[8][j].items()})
            for s,i in enumerate(w):v=add(v,{14*idx[w[:s]+w[s+1:]]+k:(-1)**(p+s+1)*c for k,c in A[i][j].items()})
            out.append(v)
    return out
def split_map(n):
    out=[]
    indices={k:{w:i for i,(p,w) in enumerate(blocks(n-k,ACTIVE))} for k in range(min(n,3)+1)}
    for p,w in blocks(n,list(range(8))):
        a=tuple(i for i in w if i in ACTIVE);z=tuple(i for i in w if i in ZERO)
        sign=(-1)**sum(i>j for i in a for j in z)
        for j in range(14):out.append((z,14*indices[len(z)][a]+j,sign))
    return out
def preflight(out):
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);hashes={}
    def read(p):hashes[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
    c=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json')
    first=read('reports/AMBIENT_B18_ELL2_C2_FIRST_BETTI_SOURCE_PROFILE_20261105_v63/actual-d2-origin-rank-and-source-factorization.json')
    A=[[parse(v) for v in M] for M in c['C2_coordinate_actions9_sparse_columns']]
    assert all(not v for i in ZERO for v in A[i])
    D={n:differential(n,ACTIVE,A) for n in [1,2,3]}
    assert [len(D[n]) for n in [1,2,3]]==[84,224,364]
    assert [sum(map(len,D[n])) for n in [1,2,3]]==[37,193,447]
    assert all(not mv(D[1],v) for v in D[2]) and all(not mv(D[2],v) for v in D[3])
    split={n:split_map(n) for n in [1,2,3]}
    for n in [2,3]:
        full=differential(n,list(range(8)),A)
        if n==2:assert full==list(map(parse,first['all_d2_constant_columns']))
        for j,col in enumerate(full):
            z,cj,sj=split[n][j];expected={}
            if len(z)<n:
                for i,a in D[n-len(z)][cj].items():expected[(z,i)]=sj*a
            got={}
            for i,a in col.items():
                zi,ci,si=split[n-1][i];got[(zi,ci)]=si*a
            assert got==expected,'literal signed zero-wedge split'
    assert first['rank']==89
    r1=11;r2=89-3*r1;assert r2==56
    bound={'basis224_times168':224*168,'all364factor168':364*168,'pivot_operations_max':168*167//2,
      'input447':447,'fixed_original_minor_sparse_max':447,'source_metadata_reserve':10000}
    assert sum(bound.values())<150000
    evidence={'active_variables':ACTIVE,'zero_actions':ZERO,
      'active_d1_columns':[srec(v) for v in D[1]],'active_d2_columns':[srec(v) for v in D[2]],'active_d3_columns':[srec(v) for v in D[3]],
      'signed_zero_wedge_decomposition':{str(n):[[list(z),i,s] for z,i,s in split[n]] for n in [1,2,3]},
      'decomposition_basis_convention':'active exterior factors first, zero exterior factors last; sign=(-1)^#(activeindex>zeroindex)',
      'all304_and1137_source_entries_reconstructed_exactly':True,
      'active_r1':r1,'active_r2':r2,'active_r3_upper':224-r2,
      'global_rank_d3_formula':'active_r3+3*56+3*11=active_r3+201',
      'global_beta2_formula':'429-globalrankd3=228-active_r3',
      'source_complete_complex':'periodic14 plusKoszul8; not a coefficient or finitekernel query',
      'whole_guard_bound':bound,'whole_guard_slots':sum(bound.values()),
      'limits':{'seconds':60,'bytes':2_000_000_000,'artifacts':32*1024*1024,'total_slots':150000,'coefficient_bits':384}}
    ep=out/'active-five-d3-source-split-and-budget.json';ep.write_text(json.dumps(evidence,separators=(',',':'))+'\n')
    report={'id':PREF,'all_passed':True,'one_new_rank_shape':[224,364],'input_nonzeros':447,'rank_upper':168,
      'whole_guard_slots':sum(bound.values()),'numeric_rank_started':False,'input_sha256':hashes,'code_sha256':sha(Path(__file__)),
      'note_sha256':sha(ROOT/'research/notes'/f'{PREF}.md'),'evidence_sha256':{ep.name:sha(ep)}}
    (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps(report))
def rank(out,gate):
    start=time.monotonic();out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
    def alarm(signum,frame):raise RuntimeError('60second guard')
    signal.signal(signal.SIGALRM,alarm);signal.alarm(60)
    gp=ROOT/gate;g=json.loads(gp.read_text());assert g['all_passed'] and not g['numeric_rank_started']
    ip=gp.parent/'active-five-d3-source-split-and-budget.json';inp=json.loads(ip.read_text());assert sha(ip)==g['evidence_sha256'][ip.name]
    cols=list(map(parse,inp['active_d3_columns']));d2=list(map(parse,inp['active_d2_columns']));assert len(cols)==364
    assert all(not mv(d2,v) for v in cols)
    piv={};receipts=[];factors=[];current_column=0
    def cguard(v):
        assert all(max(abs(c.numerator).bit_length(),c.denominator.bit_length())<=384 for c in v.values()),'384bit guard'
        return v
    def liveguard():
        slots=sum(map(len,piv.values()))+sum(len(v['subtract']) for v in receipts)+sum(map(len,factors))+2*447+10000
        assert slots<150000,'150k whole support guard';return slots
    try:
        for j,orig in enumerate(cols):
            current_column=j;v=dict(orig);steps=[]
            while v:
                p=min(v)
                if p not in piv:
                    c=v[p];piv[p]=cguard(scale(v,1/c));receipts.append({'source_column':j,'pivot':p,'divide_by':str(c),'subtract':[[i,str(a)] for i,a in steps]});break
                c=v[p];steps.append((p,c));v=cguard(add(v,piv[p],-c))
            assert len(piv)<=168;liveguard()
        porder=sorted(piv);pi={p:i for i,p in enumerate(porder)};B=[piv[p] for p in porder]
        for j,orig in enumerate(cols):
            current_column=j;v=dict(orig);c={}
            while v:
                p=min(v);assert p in piv;a=v[p];c[pi[p]]=a;v=cguard(add(v,piv[p],-a))
            factors.append(cguard(c));assert mv(B,c)==orig;liveguard()
        selected=[x['source_column'] for x in receipts];sequence=[x['pivot'] for x in receipts]
        inversions=sum(sequence[i]>sequence[j] for i in range(len(sequence)) for j in range(i+1,len(sequence)))
        rowpos={p:i for i,p in enumerate(porder)}
        minor=[{rowpos[i]:v for i,v in cols[j].items() if i in rowpos} for j in selected]
        # Keep the exact determinant factored; never manufacture a >384-bit QQ.
        det_product={'sign':(-1)**inversions,'nonzero_QQ_factors':[x['divide_by'] for x in receipts]}
        r=len(piv);beta2=228-r
        evidence={'shape':[224,364],'rank':r,'normalized_independent_columns':[srec(v) for v in B],
          'fixed_pivot_receipts':receipts,'all364_column_factorization':[srec(v) for v in factors],
          'fixed_original_minor_rows':porder,'fixed_original_minor_columns':selected,'fixed_original_minor_sparse_columns':[srec(v) for v in minor],
          'fixed_minor_exact_nonzero_determinant_product':det_product,
          'minor_proof':'normalized pivot-row matrix triangular diagonal1; exact source-column operations give the stated signed product',
          'full_rank_d3':r+201,'actual_beta2_C2_over_B':beta2,
          'new_honest_T2_generator_count_if_unit_transport_applied':beta2+30,
          'all_new_source_coefficient_queries':False,'T2_length_or_H0_computed':False}
        ep=out/'active-d3-rank-and-all-source-factorizations.json';ep.write_text(json.dumps(evidence,separators=(',',':'))+'\n')
        assert ep.stat().st_size<32*1024*1024
        report={'id':STEM,'all_passed':True,'active_d3_rank':r,'full_d3_rank':r+201,'actual_beta2_C2_over_B':beta2,
          'fixed_original_minor_nonzero':True,'all364_column_factorizations_exact':True,'whole_named_slots':liveguard(),
          'new_physical_coefficients_queried':False,'input_sha256':{str(gp.relative_to(ROOT)):sha(gp),str(ip.relative_to(ROOT)):sha(ip)},
          'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
          'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start}
        (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps(report))
    except Exception as e:
        cp={'error':str(e),'column':current_column,'pivots':[[p,srec(v)] for p,v in sorted(piv.items())],
          'receipts':receipts,'completed_factorizations':[srec(v) for v in factors],'traceback':traceback.format_exc(),'no_automatic_retry':True}
        (out/'guarded-frontier.json').write_text(json.dumps(cp,separators=(',',':'))+'\n');raise
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--rank-source-gate');a=p.parse_args()
    if a.rank_source_gate:rank(a.output,a.rank_source_gate)
    else:preflight(a.output)
