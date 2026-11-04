#!/usr/bin/env python3
"""Independent univariate coefficient/sensitivity audit of T1/T2 phi.

Uses previously independently checked source2 coefficients. No producer replay,
Newton solve, new source jet, point, prime, or full honest2 matrix.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
from ambient_B18_top_graph_telescoping_review_20261104_v63 import solve
from ambient_B18_T3_zero_structural_20261104_v63 import insert_basis

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T1_T2_QUADRATIC_PHI_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def pa(a,b):
    c=dict(a)
    for m,v in b.items():
        q=c.get(m,F(0))+v
        if q:c[m]=q
        else:c.pop(m,None)
    return c
def ps(a,c):return {m:v*c for m,v in a.items() if v*c}
def pm(a,b):
    c={}
    for m,v in a.items():
        for n,w in b.items():
            p=tuple(sorted(m+n));q=c.get(p,F(0))+v*w
            if q:c[p]=q
            else:c.pop(p,None)
    return c
def pc(c):return {():F(c)} if c else {}
def ja(a,b):return (a[0]+b[0],pa(a[1],b[1]),pa(a[2],b[2]))
def jn(a):return (-a[0],ps(a[1],-1),ps(a[2],-1))
def jm(a,b):return (a[0]*b[0],pa(ps(a[1],b[0]),ps(b[1],a[0])),pa(pa(ps(a[2],b[0]),ps(b[2],a[0])),pm(a[1],b[1])))
def ji(a):
    assert a[0]
    return (1/a[0],ps(a[1],-1/a[0]**2),pa(ps(pm(a[1],a[1]),1/a[0]**3),ps(a[2],-1/a[0]**2)))

def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    base=ROOT/'reports/AMBIENT_B18_T1_T2_QUADRATIC_PHI_20261104_v63'
    pp=base/'quadratic-phi-certificate.json';pr=base/'report.json'
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    fp=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json'
    hp=ROOT/'reports/AMBIENT_B18_T3_HONEST10_PRESENTATION_20261104_v63/honest10-presentation-DAG.json'
    up=ROOT/'reports/AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63/unit-reduced-presentation-DAG.json'
    jp=ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json'
    qp=ROOT/'reports/AMBIENT_B18_T1_T2_FIRSTJET_20261104_v63/Q1-Q2-input-and-firstjets.json'
    zp=ROOT/'reports/AMBIENT_B18_T3_ZERO_STRUCTURAL_20261104_v63/report.json'
    ze=zp.parent/'independent-firstjet-and-structural-proof.json'
    np=ROOT/'reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63/actual-N-R0-binding-DAG.json'
    nr=ROOT/'reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63_REVIEW/report.json'
    cert,old,full,honest,unit,jet,qdata,indep,N=map(read,[pp,sp,fp,hp,up,jp,qp,ze,np])
    assert all(read(p)['all_passed'] for p in [pr,zp,nr])
    nodes=unit['DAG_scalar_nodes'];oldnodes=old['DAG_scalar_nodes'];nold=len(oldnodes)
    mon=list(map(tuple,jet['source2_monomials']))
    variables={n[1]:i for i,n in enumerate(oldnodes) if n[0]=='v'}
    xids=[variables[f'x{r}'] for r in old['physical_zero_base_rows9']]
    xi=old['physical_zero_base_rows9'].index(61)
    bindings={};lines={}
    for name in ['CA','CB']:
        src=jet['both57_source2'][name]
        sol={i:list(map(F,v)) for i,v in src['all57_source2jets']}
        L=[[sol[i][1+j] for j in range(9)] for i in xids]
        y,_=solve(L,[[F(i==xi)] for i in range(9)]);d=[r[0] for r in y]
        def qvalue(v):return sum((v[k]*d[a]*d[b] for k,(a,b) in enumerate(mon[10:],10)),F(0))
        y,_=solve(L,[[-qvalue(sol[i])] for i in xids]);a=[r[0] for r in y]
        source={i:[F(0),sum(v[1+j]*d[j] for j in range(9)),sum(v[1+j]*a[j] for j in range(9))+qvalue(v)] for i,v in sol.items()}
        assert all(source[i]==[F(0),F(j==xi),F(0)] for j,i in enumerate(xids))
        saved=cert['both_family_exact_second_coordinate_contractions'][name]
        assert d==list(map(F,saved['physical_x61_line_Y_first']))
        assert a==list(map(F,saved['physical_x61_line_Y_second']))
        assert source=={i:list(map(F,v)) for i,v in saved['all57_source_coefficient_line_jets']}
        for i,v in source.items():
            if i in bindings:assert bindings[i]==v
            bindings[i]=v
        lines[name]={'first':list(map(str,d)),'second':list(map(str,a))}
    YA=list(map(F,lines['CA']['first']))
    original=[]
    for n in oldnodes:
        o=n[0]
        if o=='c':v=(F(n[1],n[2]),{},{})
        elif o=='v':a=bindings[variables[n[1]]];v=(a[0],pc(a[1]),pc(a[2]))
        elif o=='a':v=ja(original[n[1]],original[n[2]])
        elif o=='n':v=jn(original[n[1]])
        elif o=='m':v=jm(original[n[1]],original[n[2]])
        elif o=='i':v=ji(original[n[1]])
        else:raise ValueError(o)
        original.append(v)
    # Independently derive full-H correction: b_H,lin=b_phys,lin-a_CA,x61*b61,lin.
    grad=[list(map(F,r)) for r in N['all150_actual_target_odd_linear_on_Cohen8']]
    selected=N['chosen8_odd_target_physical_rows'];M=[grad[i] for i in selected]
    Mi,detM=solve(M,[[F(i==j) for j in range(8)] for i in range(8)])
    assert detM==F(-4,9)
    rows=[None,18,19,28,29,67,72,46,95];cfirst=[F(0)];correction_proof=[]
    ca=old['families']['CA']
    for row in rows[1:]:
        roots=[ca['complete_normalized_physical_w_coefficients'][d][row][0] for d in range(4)]
        a1=sum((original[i][1].get((),F(0)) for i in roots),F(0))
        bh=[grad[row][k]-a1*grad[61][k] for k in range(8)]
        coef=sum((bh[k]*Mi[k][0] for k in range(8)),F(0));cfirst.append(coef)
        correction_proof.append({'row':row,'CA_even_x61_derivative':str(a1),
          'full_H_odd_linear_on_R8':list(map(str,bh)),'z1_Hadamard_coefficient_origin':str(coef)})
    assert cfirst==[F(0)]*7+[F(1),F(-1,2)]
    assert cfirst==list(map(F,cert['BN_firstrow_linear_v_coefficients_from_FULL_composed_H_image']))

    # Separate t-degree from coefficient polynomials in formal gamma/theta.
    # This is a different arithmetic representation from the author's weighted
    # monomials involving t, gamma, theta in a single polynomial.
    known={};old_r=[q['node'] for q in full['all_coordinate_divisions']]
    gnodes=[q['node'] for q in full['all4_known_regular_g_divisions']]
    for node in old_r:
        a=list(map(F,jet['all784_CB_h_firstjets_in_YA'][str(node)]))
        first=sum((a[1+j]*YA[j] for j in range(9)),F(0))
        known[node]=(a[0],pc(first),{(1000000+node,):F(1)})
    for j,node in enumerate(gnodes):known[node]=(F(0),{(j,):F(1)},{(1000000+node,):F(1)})
    for key,a in jet['all792_CA_h_firstjets'].items():
        a=list(map(F,a));known[int(key)]=(a[0],pc(sum((a[1+j]*YA[j] for j in range(9)),F(0))),{})
    known[full['CB_g_restriction_node']]=(F(0),{},{})
    active=unit['remaining_D3_relation_original_indices'];A=unit['all_reduced_D3_relation_columns']
    mapping={}
    for n,fam in qdata['families'].items():
        raw=full['actual_Qn_R0_squared_generator_matrices'][n]
        maps=fam['all_raw_pair_to_Q3_column_indices']
        assert len(raw)==len(maps)
        assert all(pair==full['actual_Qn_R0_squared_generator_matrices']['3'][j] for pair,j in zip(raw,maps))
        ids=fam['honest9_column_source_positions_in_Q3_reduced']
        assert len(ids)==len(raw)+35
        pivot=fam['unit_graph52_column_index_in_raw_Q']
        expected_indices=list(range(36))+[36+j for k,j in enumerate(maps) if k!=pivot]
        active_inverse={j:i for i,j in enumerate(active)}
        assert ids==[active_inverse[j] for j in expected_indices]
        assert all(c==A[i] for c,i in zip(fam['honest9_relation_columns'],ids))
        assert all(list(map(F,v))==list(map(F,indep['all2069_rebuilt_exact_firstjet90_columns'][i])) for v,i in zip(fam['all_exact_firstjet90_columns'],ids))
        for i in ids:mapping[i]=active[i]
        assert raw[pivot]==full['actual_Qn_R0_squared_generator_matrices']['3'][1353]
    rawroots={};c5roots={}
    for i,j in mapping.items():
        rawroots[i]=honest['all36_Koszul_kernel_columns_R0_10'][j][0] if j<36 else full['actual_Qn_R0_squared_generator_matrices']['3'][j-36][0]
        c5roots[i]=honest['all36_Koszul_kernel_columns_R0_10'][j][5] if j<36 else honest['all2034_Q3_columns_lifted_R0_10'][j-36][5]
    step=unit['all_constant_unit_pivot_steps'][0];pe=step['relation_before_projection'][0];unode=step['pivot_node']
    assert all(i<len(full['DAG_scalar_nodes']) for i in list(rawroots.values())+[pe])
    assert all(int(i)>=len(full['DAG_scalar_nodes']) for i in jet['all792_CA_h_firstjets'])
    roots=list(rawroots.values())+list(c5roots.values())+[pe,unode]
    needed=set();stack=roots[:]
    while stack:
        i=stack.pop()
        if i in needed:continue
        needed.add(i)
        if i<nold or i in known:continue
        n=nodes[i];o=n[0]
        if o in ['a','m']:stack+=n[1:]
        elif o in ['n','i']:stack.append(n[1])
        elif o!='c':raise ValueError((i,n))
    V={}
    for i in sorted(needed):
        if i<nold:v=original[i]
        elif i in known:v=known[i]
        else:
            n=nodes[i];o=n[0]
            if o=='c':v=(F(n[1],n[2]),{},{})
            elif o=='a':v=ja(V[n[1]],V[n[2]])
            elif o=='n':v=jn(V[n[1]])
            elif o=='m':v=jm(V[n[1]],V[n[2]])
            elif o=='i':v=ji(V[n[1]])
            else:raise ValueError(o)
        V[i]=v
    p0,p1,p2=V[pe];u0,u1,_=V[unode]
    assert p0==0 and not p1 and u0==F(10,3)
    records={};phis={}
    for i in sorted(mapping):
        r0,r1,r2=V[rawroots[i]];s0,s1,_=V[c5roots[i]]
        assert not r0
        red=pa(pa(r2,ps(p2,-s0/u0)),pa(ps(pm(p1,s1),-1/u0),ps(pm(p1,u1),s0/u0**2)))
        first=list(map(F,indep['all2069_rebuilt_exact_firstjet90_columns'][i]))
        vel=[sum((first[10*k+1+j]*YA[j] for j in range(9)),F(0)) for k in range(9)]
        assert vel[0]==0
        corr=sum((cfirst[j]*vel[j] for j in range(1,9)),F(0))
        phi=pa(red,pc(corr));assert all(m==() for m in phi)
        phis[i]=phi.get((),F(0));records[i]={'old_E0_second':str(red.get((),F(0))),
          'full_H_BN_correction':str(corr),'exact_phi':str(phis[i]),'all_gamma_theta_coefficients_zero':True}
    expected={r['source_position_in_Q3_reduced']:r for r in cert['all_relation_phi_records']}
    assert set(expected)==set(mapping)
    assert all(phis[i]==F(expected[i]['exact_phi']) for i in mapping)
    results={}
    for n,fam in qdata['families'].items():
        ids=fam['honest9_column_source_positions_in_Q3_reduced']
        augmented=[list(map(F,indep['all2069_rebuilt_exact_firstjet90_columns'][i]))+[phis[i]] for i in ids]
        saved=cert['families'][n]
        assert augmented==[list(map(F,r)) for r in saved['all_augmented91_columns']]
        basis=[]
        for a in augmented:insert_basis(basis,a)
        assert len(basis)==44
        columns=saved['independent_column_indices'];coords=saved['minor_coordinate_indices']
        minor=[[augmented[i][j] for j in coords] for i in columns]
        _,det=solve(minor,[[F(i==j) for j in range(44)] for i in range(44)])
        assert det==F(saved['minor_determinant'])==F(-390625,279936)
        results[n]={'columns':len(ids),'rank':44,'minor_determinant':str(det),
          'all_exact_phi':list(map(str,[phis[i] for i in ids])),
          'minor_coordinate_indices':coords,'minor_original_column_indices':columns,
          'minor_matrix':[[str(x) for x in r] for r in minor]}
    ev={'both_source2_physical_line_inverse_coefficients':lines,
        'independent_full_H_BN_firstrow_derivation':correction_proof,
        'full_H_BN_firstrow_v_coefficients':list(map(str,cfirst)),
        'all1961_independent_phi_records':{str(k):v for k,v in records.items()},'families':results,
        'formal_unknowns_encoding':'univariate t jets with coefficient polynomials in gamma0..3 and theta_(1000000+node)',
        'all788_theta_and_four_gamma_sensitivities_cancel_exactly':True,
        'mathematical_application':'rank44 spans Jinf/mJinf, Nakayama gives Q1=Q2=Jinf, hence T1=T2=0',
        'only_T0_remaining_and_full_depth_open':True}
    ep=out/'independent-univariate-phi-and-minors.json';ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',pp,pr,sp,fp,hp,up,jp,qp,zp,ze,np,nr,
      ROOT/'code/ambient_B18_top_graph_telescoping_review_20261104_v63.py',ROOT/'code/ambient_B18_T3_zero_structural_20261104_v63.py']
    report={'id':STEM,'all_passed':True,'status':'independent_actual_T1_T2_zero_by_univariate_exact_phi_and_complete_minimal44',
      'own_second_physical_coordinate_inverse_and_all57_series_contractions':True,
      'BN_correction_independently_derived_by_full_H_chain_rule_not_copied':True,
      'all788_theta_and_four_gamma_sensitivities_cancel':True,
      'Q1_augmented_rank':44,'Q2_augmented_rank':44,'both_exact44_minor':'-390625/279936',
      'actual_Q1_Q2_equal_Jinf_T1_T2_zero':True,'w_torsion_only_degree0_hence_annihilated_by_w':True,
      'T0_nonzero_full_depth_or_nil_solution_claimed':False,
      'no_new_source_IFT_jets_GB_points_primes_or_author_rerun':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
      'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start,
      'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'certificate_bytes':ep.stat().st_size}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['all_passed','Q1_augmented_rank','Q2_augmented_rank','both_exact44_minor','seconds','max_rss_kib','certificate_bytes']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
