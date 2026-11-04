#!/usr/bin/env python3
"""Exact small mu(T) assembly from actual X0/X1; no coefficient evaluation.

Stage0 persists hashes/schema while native128 firsts are pending. Stage1 reads
both frozen first arrays, uses the saved27 constant cancellation combinations,
and computes one sparse84x741 QQ rank with complete receipts/factorization.
"""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_COMPARISON_T_CONSTANT_RANK_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rec(xs):return {i:F(c) for i,c in xs if F(c)}
def addscaled(v,w,c):
    for i,a in w.items():
        b=v.get(i,F(0))+c*a
        if b:v[i]=b
        else:v.pop(i,None)
def encode(v):return [[i,str(c)] for i,c in sorted(v.items()) if c]
def sparse_rank(columns,nrows):
    pivots=[];atrow={};selected=[];ops=[];factors=[];divides=[]
    for j,col in enumerate(columns):
        v=dict(col);coeff={};steps=[]
        while v and min(v) in atrow:
            p=atrow[min(v)];c=v[min(v)];steps.append([p,str(c)])
            coeff[p]=coeff.get(p,F(0))+c;addscaled(v,pivots[p],-c)
        r={'column':j,'subtract':steps}
        if v:
            i=min(v);c=v[i];p=len(pivots);v={k:a/c for k,a in v.items()}
            atrow[i]=p;pivots.append(v);selected.append(j);coeff[p]=c;divides.append(c)
            r.update(pivot_row=i,divide=str(c),normalized=encode(v))
        else:r['zero']=True
        factors.append(coeff);ops.append(r)
    # Both lower rank (distinct leading rows) and all-column upper span.
    for col,coeff in zip(columns,factors):
        value={}
        for p,c in coeff.items():addscaled(value,pivots[p],c)
        assert value==col
    determinant=F(1)
    for c in divides:determinant*=c
    return {'rank':len(pivots),'basis_normalized_columns':list(map(encode,pivots)),
      'selected_original_columns':selected,'selected_minor_rows':[r['pivot_row'] for r in ops if 'pivot_row' in r],
      'selected_minor_determinant':str(determinant),'all_column_factorization_in_normalized_basis':list(map(encode,factors)),
      'fixed_elimination_receipts':ops}

def main(out,native_report,native_data):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=True)
    assert not(out/'report.json').exists(),'frozen result already exists'
    paths={
      'X0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',
      'first640':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',
      'resolution':'reports/AMBIENT_B18_CBAR_LITERAL_RESOLUTION_20261104_v63/actual-Cbar-minimal-resolution-and-three-targets.json',
      'fullV':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
      'genericK':'reports/AMBIENT_B18_CRITICAL_CONTACT_AND_COMPARISON_RANK4_20261104_v63/report.json'}
    data={k:read(ROOT/p) for k,p in paths.items()};X0,A,R,V,G=[data[k] for k in paths]
    source_peers={
      'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63_REVIEW/report.json':'baaec4550a00c25ee24b2ab5d03e9eebaf67943caa216a01abe624bf19f3dd17',
      'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63_REVIEW/report.json':'0000a1aa2790891a329f7ff708ef6dcfcaea8b04571de07b79ee5b503caabcb1'}
    for p,h in source_peers.items():assert sha(ROOT/p)==h and read(ROOT/p)['all_passed']
    assert G['all_passed'] and G['rank_D_K']==4 and G['dim_k_C_ell_torsion']==4
    assert G['length_k_T_lower_bound']==8 and not G['origin_depth_decided']
    for p,h in G['input_sha256'].items():assert sha(ROOT/p)==h,p
    assert X0['rank_X0']==27
    combinations=[r for r in X0['all768_fixed_elimination_receipts_and_constant_combinations'] if r.get('zero')]
    assert len(combinations)==741
    # Persist this fixed source-only injection check while the derivative
    # producer runs. Reuse it at execution rather than rerunning the check.
    anchors_path=out/'literal-d2-linear-anchors.json'
    if anchors_path.exists():
        anchors_data=read(anchors_path)
        assert anchors_data['literal_resolution_sha256']==sha(ROOT/paths['resolution'])
        basis=anchors_data['literal_wedge_basis84'];d2linear=list(map(rec,anchors_data['literal_d2_linear_columns84']))
    else:
        basis=R['resolution_basis_sector_and_wedge'][2];d2linear=[]
        for column in R['differentials8_column_sparse'][1]:
            entry={}
            for row,p in column:
                for m,c in p:
                    if sum(m)==1:entry[8*row+m.index(1)]=F(c)
            d2linear.append(entry)
        anchors=[8*(8*s+j)+i for s,(i,j) in basis]
        assert len(d2linear)==len(basis)==len(set(anchors))==84
        assert [[v.get(k,F(0)) for v in d2linear] for k in anchors]==[[F(i==j) for j in range(84)] for i in range(84)]
        anchors_data={'literal_resolution_sha256':sha(ROOT/paths['resolution']),
          'literal_wedge_basis84':basis,'literal_d2_linear_columns84':list(map(encode,d2linear)),
          'selected_derivative_coordinates84':anchors,'selected84_matrix_equals_identity_exact':True}
        anchors_path.write_text(json.dumps(anchors_data,separators=(',',':'))+'\n')
    phase={'phase':'waiting_for_frozen_native128_X1_PASS','source_hashes':{p:sha(ROOT/p) for p in paths.values()}|source_peers,
      'expected_native_schema':'source_columns128_rows48_sparse_dual9; each column [[row,[[dualindex,QQ],...]],...]; dual0constant,1..7zA0..6,8rA',
      'expected_assembly':[84,741],'constant_rank_r':27,'post_constant_target_count':89,
      'no_348_matrix_or_new_coefficient_query':True,'literal84_d2_injective_anchors_verified':True,
      'literal_anchor_cache_sha256':sha(anchors_path)}
    (out/'phase.json').write_text(json.dumps(phase,indent=2)+'\n')
    if native_report is None:
        print(json.dumps({'status':'PREPARED_waiting_native128_PASS','shape':[84,741],'r':27}));return
    assert native_data is not None
    native_report=native_report.resolve();native_data=native_data.resolve();nr=read(native_report);N=read(native_data)
    assert nr['all_passed']
    # The caller must bind the exact native data through its frozen report.
    def contains_hash(obj,h):
        if isinstance(obj,dict):return any(contains_hash(v,h) for v in obj.values())
        if isinstance(obj,list):return any(contains_hash(v,h) for v in obj)
        return obj==h
    assert contains_hash(nr,sha(native_data)),'native data is not pinned by its report'
    ns=N['source_columns128_rows48_sparse_dual9']
    assert len(ns)==128
    first=A['source_columns640_rows56_sparse_dual9'];assert len(first)==640
    # Strict actual sparse-source budget before constructing Fraction arrays.
    source_entries=sum(len(p) for col in first+ns for row,p in col)
    combo_entries=sum(len(r['constant_original_source_combination']) for r in combinations)
    guard={'native_and_first640_source_nonzeros':source_entries,'constant_combination_entries':combo_entries,
      'single_linear_scratch_bound':24*8,'max_w_matrix_slots':84*741,
      'sparse_only':True,'expanded_coefficient_cap':150000}
    assert source_entries+combo_entries+84*741+24*8<150000
    (out/'preflight.json').write_text(json.dumps(guard,indent=2)+'\n')
    X=[]
    for j,col in enumerate(first+ns):
        entries={i:rec(p) for i,p in col};assert len(entries)==len(col)
        assert all(0<=i<(56 if j<640 else 48) for i in entries)
        for row in range(56):assert entries.get(row,{}).get(0,F(0))==F(X0['full_X0_56_by768'][row][j])
        X.append(entries)
    m0=[list(map(F,row)) for row in V['all24_by56_m0']]
    # Literal d2 signs and linear terms; actual CB quadratic shifts have no
    # effect on this first derivative, and r^n has no linear term.
    assert len(d2linear)==len(basis)==84 and anchors_data['selected84_matrix_equals_identity_exact']
    wcols=[];a1_columns=[];source_ids=[];constraint_entries=0
    for relation in combinations:
        combo=rec(relation['constant_original_source_combination']);xprime={}
        for j,c in combo.items():
            for row,p in X[j].items():
                target=xprime.setdefault(row,{})
                for k,a in p.items():
                    b=target.get(k,F(0))+c*a
                    if b:target[k]=b
                    else:target.pop(k,None)
        assert all(not p.get(0,F(0)) for p in xprime.values())
        v={}
        for row,coeffs in enumerate(m0):
            for i,c in enumerate(coeffs):
                if c and i in xprime:
                    for k,a in xprime[i].items():
                        if k:
                            key=8*row+k-1;v[key]=v.get(key,F(0))+c*a
        v={i:c for i,c in v.items() if c}
        w={}
        for k,(s,pair) in enumerate(basis):
            i,j=pair
            q=v.get(8*(8*s+7)+i,F(0)) if j==7 else (v.get(8*(8*s+j)+i,F(0))-v.get(8*(8*s+i)+j,F(0)))/2
            if q:w[k]=q
        check={}
        for k,c in w.items():addscaled(check,d2linear[k],c)
        assert check==v,'failed actual linear cycle constraint (includes all diag/skew/r entries)'
        constraint_entries+=24*8;wcols.append(w);a1_columns.append(encode(v));source_ids.append(relation['source_column'])
    phase['phase']='all741_actual_linear_cycle_constraints_checked';phase['checked_linear_entries']=constraint_entries
    (out/'phase.json').write_text(json.dumps(phase,indent=2)+'\n')
    rank=sparse_rank(wcols,84);rho=rank['rank'];mu=89-rho
    # Dynamic total support includes the coefficient matrix, factorization,
    # normalized pivots, all saved source-linear readbacks and inputs.
    actual_support=source_entries+combo_entries+sum(map(len,wcols))+sum(map(len,rank['basis_normalized_columns']))+sum(map(len,rank['all_column_factorization_in_normalized_basis']))+sum(map(len,a1_columns))
    assert actual_support<150000
    ev={'actual_remaining_source_column_labels':source_ids,
      'actual_constant_source_combinations':[r['constant_original_source_combination'] for r in combinations],
      'actual24_by8_linear_cycles_sparse':a1_columns,'literal_d2_linear_columns84':list(map(encode,d2linear)),
      'literal_wedge_basis84':basis,'actual84_by741_w0_columns_sparse':list(map(encode,wcols)),
      'all142272_linear_constraint_entries_checked':constraint_entries==741*24*8,
      'Xprime_constant_terms_all_zero':True,'m1_cancels_exactly':True,
      'complete_rank_certificate':rank,'r':27,'rank_w0':rho,'mu_Rbar_T_exact':mu,
      'length_T_lower_bound':mu,'mu_D_K_lower_bound':max(0,mu-4),
      'depth19_proved_by_lower_bound':mu>8,'depth20_not_inferred_if_lower_bound_insufficient':True,
      'no_T_length_or_348_matrix_computed':True,'expanded_live_support':actual_support}
    ep=out/'actual-linear-cycle-muT-certificate.json';ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
    phase['phase']='FINAL_actual_mu_T_computed';phase['mu_T']=mu;(out/'phase.json').write_text(json.dumps(phase,indent=2)+'\n')
    report={'id':STEM,'all_passed':True,'r_X0':27,'actual_rank_w0':rho,'actual_mu_Rbar_T':mu,
      'QQ_length_T_lower_bound':mu,'mu_D_K_lower_bound':max(0,mu-4),'origin_depth19_proved_by_this_bound':mu>8,
      'inconclusive_for_T_length_if_mu_at_most8':mu<=8,
      'all768_source_origins_and741_linear_cycle_constraints_bound':True,
      'no_new_coefficients_source_jet_point_GB_348_matrix':True,
      'input_sha256':{p:sha(ROOT/p) for p in paths.values()}|source_peers|{str(native_report.relative_to(ROOT)):sha(native_report),str(native_data.relative_to(ROOT)):sha(native_data)},
      'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['id','all_passed','actual_rank_w0','actual_mu_Rbar_T','origin_depth19_proved_by_this_bound','seconds']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--native-report',type=Path);p.add_argument('--native-data',type=Path)
    a=p.parse_args();main(a.output,a.native_report,a.native_data)
