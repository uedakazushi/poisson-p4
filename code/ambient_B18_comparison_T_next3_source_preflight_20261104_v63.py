#!/usr/bin/env python3
"""Fixed rank-receipt extraction for the three unresolved quadratic T classes.

No source coefficient, source jet, scalar inverse, or fresh rank search.  The
620 combinations kill the ACTUAL minimal linear relations before any second
coefficient functional is considered.
"""
import argparse,itertools,json,resource,signal,sys,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_COMPARISON_T_NEXT3_SOURCE_PREFLIGHT_20261104_v63'

def compose(c,columns):
    v={}
    for j,a in c.items():addscaled(v,columns[j],a)
    return v

def receipt_combos(cert,n):
    normalized=[];zero=[];labels=[]
    for j,op in enumerate(cert['fixed_elimination_receipts']):
        assert j==op['column'];v={j:F(1)}
        for p,a in op['subtract']:addscaled(v,normalized[p],-F(a))
        if 'pivot_row' in op:normalized.append({i:a/F(op['divide']) for i,a in v.items()})
        else:assert op['zero'];zero.append(v);labels.append(j)
    assert len(normalized)==cert['rank'] and len(zero)==n-cert['rank']
    return normalized,zero,labels

def write(out,name,data):
    p=out/name;p.write_text(json.dumps(data,separators=(',',':'))+'\n');return p

def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));sys.set_int_max_str_digits(0)
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'LQ':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/actual-linear-quadratic-T-certificate.json',
      'L1':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase1-actual-five-target-linear-relations.json',
      'Q2':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase2-certified-quadratic-input.json',
      'M':'reports/AMBIENT_B18_COMPARISON_T_CONSTANT_RANK_20261104_v63/actual-linear-cycle-muT-certificate.json',
      'X0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',
      'X1a':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',
      'X1b':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json',
      'V':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json'}
    hashes={p:sha(ROOT/p) for p in paths.values()}
    write(out,'phase.json',{'phase':'fixed_sources_bound_before_arithmetic','input_sha256':hashes,'no_new_source_coefficient':True})
    D={k:read(ROOT/p) for k,p in paths.items()};LQ,L1,Q2,M,X0,A,B,V=[D[k] for k in paths]
    assert LQ['actual_linear_relation_rank']==37 and LQ['quadratic_certificate']['rank']==177
    mons=[tuple(m) for m in Q2['quadratic_monomials36']];assert mons==list(itertools.combinations_with_replacement(range(8),2))
    mi={m:i for i,m in enumerate(mons)};qcert=LQ['quadratic_certificate'];qb=list(map(rec,qcert['basis_normalized_columns']))
    lead={min(v):i for i,v in enumerate(qb)};missing=[i for i in range(180) if i not in lead];assert missing==[143,176,179]
    dual=[{} for _ in missing]
    for j in range(180):
        v={j:F(1)}
        for i in sorted(lead):
            if i in v:addscaled(v,qb[lead[i]],-v[i])
        assert set(v)<=set(missing)
        for k,i in enumerate(missing):
            if v.get(i):dual[k][j]=v[i]
    for col in qb+list(map(rec,Q2['all464_certified_quadratic_columns_sparse'])):
        assert all(sum((a*col.get(j,F(0)) for j,a in d.items()),F(0))==0 for d in dual)
    assert [[d.get(i,F(0)) for i in missing] for d in dual]==[[F(i==j) for j in range(3)] for i in range(3)]
    ds=[]
    for d in dual:
        targets=[]
        for t in range(5):
            rows=[{'variables':list(mons[j%36]),'coefficient':str(a)} for j,a in sorted(d.items()) if j//36==t]
            targets.append({'target':t,'fullV_source_row':L1['five_remaining_free_source_rows'][t],'terms':rows})
        ds.append({'targets':targets,'variable_support':sorted({k for j in d for k in mons[j%36]}),'nonzero_count':len(d)})
    write(out,'three-quadratic-duals.json',{'monomials36':mons,'five_targets':L1['five_remaining_free_source_rows'],
      'quotient_representative_rows':missing,'projection3_by180_rows_sparse':list(map(encode,dual)),
      'support_by_dual_and_target':ds,'all177basis_and464columns_annihilated':True,'quotient_representatives_map_to_I3':True})
    write(out,'phase.json',{'phase':'exact_three_duals_complete','dual_support':[x['variable_support'] for x in ds]})
    # Replay both frozen column-elimination receipts, carrying exact source
    # coefficient vectors. No new rank computation is performed here.
    wn,wzero,wlabels=receipt_combos(M['complete_rank_certificate'],741)
    ln,lzero,llabels=receipt_combos(LQ['complete_linear_rank_certificate'],657)
    linear=list(map(rec,L1['actual40_by657_linear_relation_columns']))
    assert len(wzero)==657 and len(lzero)==620
    assert wlabels==[M['actual_remaining_source_column_labels'].index(x) for x in L1['actual657_minimal_relation_source_labels']]
    for c in lzero:assert not compose(c,linear)
    c741=[compose(c,wzero) for c in lzero]
    base768=list(map(rec,M['actual_constant_source_combinations']))
    c768=[compose(c,base768) for c in c741]
    xcols=A['source_columns640_rows56_sparse_dual9']+B['source_columns128_rows48_sparse_dual9']
    X=[{i:rec(p) for i,p in col} for col in xcols]
    pivot_records=[r for r in X0['all768_fixed_elimination_receipts_and_constant_combinations'] if 'pivot' in r]
    xpcombs=[rec(r['constant_original_source_combination']) for r in pivot_records]
    fullx0=[{i:F(X0['full_X0_56_by768'][i][j]) for i in range(56) if F(X0['full_X0_56_by768'][i][j])} for j in range(768)]
    xpfull=[compose(c,fullx0) for c in xpcombs]
    U=X0['free_kernel_rows32'];pb={r['pivot']:rec(r['normalized']) for r in pivot_records}
    pi_at={r['pivot']:i for i,r in enumerate(pivot_records)}
    assert all({i:v.get(row,F(0)) for i,row in enumerate(U) if v.get(row)}==rec(r['normalized']) for r,v in zip(pivot_records,xpfull))
    zero_records=[];total_support=0
    m0=[list(map(F,row)) for row in V['all24_by56_m0']]
    for number,(c657,c1,c) in enumerate(zip(lzero,c741,c768)):
        assert not compose(c,fullx0)
        first=[{} for _ in range(8)]
        for j,a in c.items():
            for row,p in X[j].items():
                for k,b in p.items():
                    if k:first[k-1][row]=first[k-1].get(row,F(0))+a*b
        alpha={}
        for k,v0 in enumerate(first):
            v0={i:a for i,a in v0.items() if a};first[k]=v0
            v={i:v0[row] for i,row in enumerate(U) if v0.get(row)};ak={}
            for pivot in sorted(pb):
                a=v.get(pivot,F(0))
                if a:ak[pi_at[pivot]]=a;addscaled(v,pb[pivot],-a)
            assert not v and compose(ak,xpfull)==v0
            assert all(sum((a*v0.get(i,F(0)) for i,a in enumerate(row)),F(0))==0 for row in m0)
            for j,a in ak.items():alpha[8*j+k]=a
        total_support+=len(c657)+len(c1)+len(c)+len(alpha)+sum(map(len,first))
        zero_records.append({'linear_zero_receipt':llabels[number],'combination_in657':encode(c657),
          'combination_in741':encode(c1),'combination_in_actual768':encode(c),
          'actual56_firsts_by_variable':list(map(encode,first)),
          'linear_Xp27_correction_row8j_plus_variable':encode(alpha)})
    assert total_support<150000
    write(out,'zero-linear620-source-combinations.json',{'coordinate_order':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA'],
      'all620':zero_records,'normalized_Xp27_original768_combinations':list(map(encode,xpcombs)),
      'all_Xc0_zero_all_Xc1_in_Xp0_span':True,'total_saved_QQ_support':total_support,
      'source_quadratic_coefficients_not_evaluated':True})
    write(out,'phase.json',{'phase':'exact_three_duals_and620_source_combinations_complete','source_support':total_support})
    # Explicit conservative adjoint sensitivity. Euler is extended linearly to
    # all homogeneous quadratics; only its restriction to genuine cycles is a
    # lift. Pure r² obstruction terms are retained as a separate source gate.
    qnorm=list(map(rec,L1['normalized84_unitcycle_Q1_columns_sparse']))
    wb=list(map(rec,M['complete_rank_certificate']['basis_normalized_columns']));wl={min(v):i for i,v in enumerate(wb)}
    wedge=M['literal_wedge_basis84'];wi={(s,tuple(pair)):i for i,(s,pair) in enumerate(wedge)}
    # Coefficients pairing dual with minus Qstar1*beta1.
    pairbeta=[]
    for d in dual:
        row={}
        for p,q in enumerate(qnorm):
            for k in range(8):
                a=-sum((b*d.get(36*(i//8)+mi[tuple(sorted((i%8,k)))],F(0)) for i,b in q.items()),F(0))
                if a:row[8*p+k]=a
        pairbeta.append(row)
    # Pi5 on full56 coefficients is the U-coordinate constant quotient by Xp.
    pi=[];tail=L1['constant_projection_tail5']
    for row in range(56):
        v={U.index(row):F(1)} if row in U else {}
        for pivot in sorted(pb):
            if pivot in v:addscaled(v,pb[pivot],-v[pivot])
        pi.append({j:v[i] for j,i in enumerate(tail) if v.get(i)})
    wbcoords=[]
    for j in range(84):
        v={j:F(1)};b={}
        while v:
            i=min(v);p=wl[i];a=v[i];b[p]=b.get(p,F(0))+a;addscaled(v,wb[p],-a)
        wbcoords.append(b)
    sensitivity=[{} for _ in range(3)];cycle_obstruction=[]
    for src in range(56):
        for mindex,(u,v) in enumerate(mons):
            # d2-linear right inverse on admissible homogeneous quadratic
            # cycles, separated by z-degree to avoid an incorrect uniform1/3.
            wc={};k=int(u<7)+int(v<7)
            for s in range(3):
                for i in range(7):
                    for j in range(i+1,7):
                        # (partial_i a_j-partial_j a_i)/(zdegree+1)
                        for derivative,outrow,sign in [(i,8*s+j,1),(j,8*s+i,-1)]:
                            multiplicity=(u==derivative)+(v==derivative)
                            if multiplicity:
                                remain=v if u==derivative else u;a=F(sign*multiplicity,k+1)*m0[outrow][src]
                                if a:idx=8*wi[(s,(i,j))]+remain;wc[idx]=wc.get(idx,F(0))+a
                    if k:
                        mult=(u==i)+(v==i)
                        if mult:
                            remain=v if u==i else u;a=F(mult,k)*m0[8*s+7][src]
                            if a:idx=8*wi[(s,(i,7))]+remain;wc[idx]=wc.get(idx,F(0))+a
            beta={}
            for index,a in wc.items():
                wedgeidx,var=divmod(index,8)
                for p,b in wbcoords[wedgeidx].items():idx=8*p+var;beta[idx]=beta.get(idx,F(0))+a*b
            for t,d in enumerate(dual):
                a=sum((b*d.get(36*target+mindex,F(0)) for target,b in pi[src].items()),F(0))
                a+=sum((b*pairbeta[t].get(i,F(0)) for i,b in beta.items()),F(0))
                if a:sensitivity[t][36*src+mindex]=a
    sinfo=[]
    for d in sensitivity:
        sinfo.append({'nonzeros':len(d),'variable_support':sorted({k for i in d for k in mons[i%36]}),
          'source_rows56':sorted({i//36 for i in d}),
          'variable_pairs':sorted({mons[i%36] for i in d})})
    write(out,'projected-source-X2-sensitivity.json',{'source_rows56_monomials36_order':True,
      'three_rows_sparse':list(map(encode,sensitivity)),'support':sinfo,
      'Euler_extension_requires_genuine_cycle_condition':True,
      'no_source_second_coefficients_evaluated':True,'m1_removed_only_after620_linear_cancellation_and_mod177':True,
      'source_correction':'Z2=sum(c_j Xj2)-sum(alpha_p1 Xp1); unknown Xp0 alpha2 annihilated',
      'Koszul_lift_ambiguity_is_d3_linear_and_annihilated_by_dual3':True})
    write(out,'phase.json',{'phase':'FINAL','residual_dual_support':[x['variable_support'] for x in ds],
      'source_X2_sensitivity_support':[x['variable_support'] for x in sinfo]})
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'fixed_dual_shape':[3,180],'zero_linear_combinations':620,
      'residual_dual_support':[x['variable_support'] for x in ds],
      'source_X2_sensitivity_support':[x['variable_support'] for x in sinfo],
      'source_X2_sensitivity_QQ_nonzeros':sum(map(len,sensitivity)),
      'total_saved_source_combinations_QQ_nonzeros':total_support,'no_new_source_jet_rank_GB_point':True,
      'T_length_or_depth_not_decided':True,'input_sha256':hashes,'code_sha256':sha(Path(__file__)),
      'note_sha256':sha(note),'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    write(out,'report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','residual_dual_support','source_X2_sensitivity_support','source_X2_sensitivity_QQ_nonzeros','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
