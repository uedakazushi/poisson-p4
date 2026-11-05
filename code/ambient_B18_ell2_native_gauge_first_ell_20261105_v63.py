#!/usr/bin/env python3
"""Selected cached-Jacobian first-ell gauge extension; no physical solve."""
import argparse,hashlib,json,re,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
ROOT=Path(__file__).resolve().parents[1]
PREF='AMBIENT_B18_ELL2_NATIVE_GAUGE_FIRST_ELL_PREFLIGHT_20261105_v63'
STEM='AMBIENT_B18_ELL2_NATIVE_GAUGE_FIRST_ELL_20261105_v63'
JPATH='reports/AMBIENT_B18_T_JOINT135_SELECTED_FUNCTIONALS_20261104_v63/phase1-actual-J135-blocks-and-target-Jacobian.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(v):return {i:F(c) for i,c in v}
def dense_mv(A,v):return [sum((a*v[j] for j,a in row.items()),F(0)) for row in A]
def guard(v):assert all(max(abs(c.numerator).bit_length(),c.denominator.bit_length())<=384 for c in v);return v
def selected_top_field(path,key):
    # Parse one UNIQUE top-level JSON field; do not retain unused N/T arrays.
    s=path.read_text();matches=list(re.finditer(re.escape(json.dumps(key))+r'\s*:',s));assert len(matches)==1
    pos=matches[0].start();depth=0;quoted=False;escaped=False
    for ch in s[:pos]:
        if quoted:
            if escaped:escaped=False
            elif ch=='\\':escaped=True
            elif ch=='"':quoted=False
        elif ch=='"':quoted=True
        elif ch in '[{':depth+=1
        elif ch in ']}':depth-=1
    assert depth==1 and not quoted
    start=matches[0].end()
    while s[start].isspace():start+=1
    return json.JSONDecoder().raw_decode(s,start)[0]
def preflight(out):
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);hs={}
    def read(p):hs[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
    J=read(JPATH)
    peer=read('reports/AMBIENT_B18_T_JOINT135_SELECTED_FUNCTIONALS_20261104_v63_REVIEW/report.json');assert peer['all_passed']
    phi=read('reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63/phase1-actual105-B56-cache-composition.json')
    oldJ=read('reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase2-fixed30-Jacobian.json')
    oldg=read('reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase3-gauge30-firsts.json')
    gp='reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase5-native-adapter-firsts.json';hs[gp]=sha(ROOT/gp);oldG=selected_top_field(ROOT/gp,'G8')
    np='reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase1-source-firsts.json';hs[np]=sha(ROOT/np);oldnu=selected_top_field(ROOT/np,'nu_first9')
    assert J['all_base_G8_and_gauge_constraints_checked'] and J['source_dependent_monic_carry_included']
    assert J['canonical_gauge30_order']==['even0..14','odd0..14'] and J['physical105_order']==phi['unknown105_order']
    P=list(map(parse,J['Jg_phi30_by105']));Q=list(map(parse,J['Jg_gamma30']));Qi=list(map(parse,J['Jg_gamma_inverse']));T=list(map(parse,J['properG32_and_nu2_target_Jacobian34_by135']))
    expected=[]
    for i in range(30):
        row={}
        for j in range(30):
            c=F(oldJ['A0'][i%15][j%15]) if (i<15)==(j<15) else (F(oldJ['A1_eta_delta'][i-15][j]) if i>=15 and j<15 else F(0))
            if c:row[j]=c
        expected.append(row)
    assert Q==expected and Qi==[{j:F(c) for j,c in enumerate(row) if F(c)} for row in oldJ['J30_inverse']]
    for j in range(30):
        v=[F(i==j) for i in range(30)];assert dense_mv(Q,dense_mv(Qi,v))==v and dense_mv(Qi,dense_mv(Q,v))==v
    basis=[tuple(m) for m in phi['basis56']];idx={m:i for i,m in enumerate(basis)};u=list(map(parse,phi['all105_cached_values_B56']))
    first=[]
    for k in range(9):
        m=[0]*9;m[k]=1;first.append([v.get(idx[tuple(m)],F(0)) for v in u])
    allold=[]
    for k in range(8):
        rhs=dense_mv(P,first[k]);gamma=guard([-c for c in dense_mv(Qi,rhs)])
        expectedg=[F(row[k+1]) for row in oldg['even15_first9']+oldg['odd15_first9']]
        assert gamma==expectedg and [a+b for a,b in zip(dense_mv(Q,gamma),rhs)]==[F(0)]*30
        target=guard(dense_mv(T,first[k]+gamma))
        expectedG=[F(oldG[g][4*(k+1)+mask]) for g in range(8) for mask in range(4)]
        assert target[:32]==expectedG and target[32:]==[F(row[k+1]) for row in oldnu]
        allold.append({'direction':k,'all30constraints_zero':True,'gauge30':list(map(str,gamma)),'G32_nu2':list(map(str,target))})
    rpath='reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json'
    S=read(rpath)['native_S12_by32']
    slots={'all_sparse_J135_selected_matrix_records':sum(len(r) for A in [P,Q,Qi,T] for r in A),
      'all_old30_Jacobian_input':15*15*2+30*30,'B56_cached105_dense_upper':105*56,
      'all9_phi_firsts':105*9,'old_gauge_input_arrays':30*8*2+30*9,'old_G8_selected_field':8*36,'old_nu_selected_field':2*9,
      'old8_response_arrays':8*(30+34),'new_ell_response_and_RHS':105+30+30+34,
      'new_S1_ell_sparse_upper':504,'source_frame_metadata_reserve':1000}
    assert sum(slots.values())<20000
    evidence={'physical105_order':J['physical105_order'],'canonical_gauge30_order':J['canonical_gauge30_order'],
      'Jg_phi_sparse_rows':J['Jg_phi30_by105'],'Jg_gamma_sparse_rows':J['Jg_gamma30'],
      'Jg_gamma_inverse_sparse_rows':J['Jg_gamma_inverse'],'G32_nu2_target_Jacobian_sparse_rows':J['properG32_and_nu2_target_Jacobian34_by135'],
      'actual_phi_first_ell':list(map(str,first[8])),'all8_old_direction_crosschecks':allold,
      'old_J30_and_inverse_bothproducts_exact':True,'source_dependent_monic_carry_included':True,
      'source_branch':'same canonical group15 positions/even15+odd15 order and saved odd center as peered JOINT135',
      'ell_direction':'sA=rA+ell/5 atfixedz7,rA, hencephi_ell=(partial_s phi)/5',
      'trace_caution':'actual p=tracefree_tau+s_trace(phi)/2; s_trace is the original source trace, never the F9coordinate sA; coherent J135alreadyincludes this',
      'native_S12_by32':S,'new_ell_gauge_and_G8_values_evaluated':False,
      'new_job_formula':'gammaell=-Jgamma_inverse*Jphi*phiell; Gnuell=Jtarget*(phiell,gammaell)',
      'new_native_output':'S1ell only; effectiveX1=L0*S1ell onlyafterseparateS0*cN=0, no individualnativeX1 or N1 claim',
      'whole_slot_bound':sum(slots.values()),'slot_breakdown':slots,
      'limits':{'seconds':60,'bytes':2_000_000_000,'artifact_bytes':32*1024*1024,'slots':20000,'coefficient_bits':384},
      'physical_implicit_higherjet_or_N1_queries':False}
    ep=out/'actual-first-ell-input-and-old8-source-crosschecks.json';ep.write_text(json.dumps(evidence,separators=(',',':'))+'\n')
    report={'id':PREF,'all_passed':True,'old8_all30_gauge_and34_target_crosschecks_exact':True,
      'new_ell_gauge_values_computed':False,'whole_slot_bound':sum(slots.values()),'input_sha256':hs,'code_sha256':sha(Path(__file__)),
      'note_sha256':sha(ROOT/'research/notes'/f'{PREF}.md'),'evidence_sha256':{ep.name:sha(ep)}}
    (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'whole_slot_bound':sum(slots.values()),'new_ell_values':False}))
def values(out,gate_path,peer_path):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);gp=ROOT/gate_path;g=json.loads(gp.read_text());assert g['all_passed']
    pp=ROOT/peer_path;p=json.loads(pp.read_text());assert p['all_passed']
    pins={}
    for k,v in p.get('input_sha256',{}).items():
        path,h=(v['path'],v['sha256']) if isinstance(v,dict) else (k,v);assert sha(ROOT/path)==h;pins[path]=h
    assert pins.get(gate_path)==sha(gp),'peer must bind exactfirstellsourcegate'
    ip=gp.parent/'actual-first-ell-input-and-old8-source-crosschecks.json';assert sha(ip)==g['evidence_sha256'][ip.name];d=json.loads(ip.read_text())
    P=list(map(parse,d['Jg_phi_sparse_rows']));Q=list(map(parse,d['Jg_gamma_sparse_rows']));Qi=list(map(parse,d['Jg_gamma_inverse_sparse_rows']));T=list(map(parse,d['G32_nu2_target_Jacobian_sparse_rows']))
    phi=list(map(F,d['actual_phi_first_ell']));rhs=guard(dense_mv(P,phi));gamma=guard([-c for c in dense_mv(Qi,rhs)])
    assert [a+b for a,b in zip(dense_mv(Q,gamma),rhs)]==[F(0)]*30
    target=guard(dense_mv(T,phi+gamma));assert target[32:]==[F(0),F(0)],'actualmonicfirstellzero'
    G=[target[4*i:4*i+4] for i in range(8)];cols=[{} for _ in range(128)]
    for row,col,k,sign in d['native_S12_by32']['koszul28_in_first8_rows']:
        for source_mask in range(4):
            for mask,c in enumerate(G[k]):
                if c and not(mask&source_mask):
                    i=4*row+(mask|source_mask);j=4*col+source_mask
                    cols[j]=add(cols[j],{i:F(sign)*c})
    assert sum(map(len,cols))<=504
    evidence={'actual_phi_first_ell':list(map(str,phi)),'RHS30':list(map(str,rhs)),
      'actual_gauge30_first_ell':list(map(str,gamma)),'all30canonical_constraints_firstell_zero':True,
      'actual_G8_C4_first_ell':[list(map(str,v)) for v in G],'actual_nuA_nuB_first_ell':list(map(str,target[32:])),
      'actual_native_S1_ell_48_by128_sparse_columns':[srec(v) for v in cols],
      'sourceframe':'literal complete normalizer native_S12_by32, existing signs andmonics; free4 columns derivativezero',
      'no_native_L1_or_N1_evaluated':True,'no_individual_actual_native_X1_claim':True,
      'effective_use_gate':'onlyonoriginalsourcecombinations withnativeS0c0=0, useL0*S1ell*cN',
      'all_new_physical_source_or_gauge_inverse_solves':False}
    ep=out/'actual30-firstell-and-native-S1ell.json';ep.write_text(json.dumps(evidence,separators=(',',':'))+'\n')
    report={'id':STEM,'all_passed':True,'new_canonical_firstell_coefficients':30,
      'actual_G8_firstell_components':32,'actual_native_S1ell_nonzeros':sum(map(len,cols)),
      'all30constraints_exact':True,'new_physical_source_or_higherjet_or_N1':False,
      'actual_T2_origin_or_length_computed':False,'input_sha256':{gate_path:sha(gp),str(ip.relative_to(ROOT)):sha(ip),peer_path:sha(pp)},
      'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
      'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start}
    (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps(report))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--value-source-gate');p.add_argument('--source-peer');a=p.parse_args()
    if a.value_source_gate:values(a.output,a.value_source_gate,a.source_peer)
    else:preflight(a.output)
