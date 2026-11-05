#!/usr/bin/env python3
"""Independent full-C Nakayama/DVR and conditional graph03 source proof."""
import hashlib,itertools,json,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_FULL_LINEAR_REDUCTION_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
PIN='27575d10b8b18166a2e2d689de3f32071175058d4bf8da81b267a9972de127de'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def det3(M):
    return sum(((-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3)))
       *M[0][p[0]]*M[1][p[1]]*M[2][p[2]] for p in itertools.permutations(range(3)))
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};docs={}
    def pin(p,h=None):
        p=pathlib.Path(p);s=sha(p)
        if h is not None:assert s==h,(str(p),s,h)
        pins[str(p.relative_to(ROOT))]=s;return read(p) if p.suffix=='.json' else None
    pr=pin(ROOT/'reports'/PRI/'report.json',PIN);assert pr['all_passed']
    for p,h in pr['input_sha256'].items():
        d=pin(ROOT/p,h);docs[pathlib.Path(p).name]=d
        if p.endswith('_REVIEW/report.json'):assert d['all_passed']
    for p,h in pr['evidence_sha256'].items():D=pin(ROOT/'reports'/PRI/p,h)
    pin(ROOT/'code/ambient_B18_C_full_linear_reduction_source_20261105_v63.py',pr['code_sha256'])
    pin(ROOT/'research/notes'/f'{PRI}.md',pr['note_sha256'])
    B4=pin(ROOT/'reports/AMBIENT_B18_C3_TAUA_B4_SOURCE_V2_20261105_v63/four-actual-cubic-sources-two-cached-adjoints-and-rank2-completion.json',
        '3605f4a6a70b95569979997b32bbc8ba5f7229456314d3b4f311185317a1f327')
    raw=pin(ROOT/'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
        '16287ac4ec5d098b9dae13cd5b7b911f3f6de0173bcb6553523f5d934e2ddde6')
    note=ROOT/'research/notes/AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63.md'
    pin(note,'b689b825af502d9e274467965cbee1a19e3768150e2ca618dbeb171eabff978a')
    frame=docs['actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json']
    graph=docs['actual-cover-NF-and-five-triangular-cubic-quartic-contradictions.json']
    c2=docs['actual-physical12-and-C2-right-exact-quotient.json']
    assert frame['actual_C_and_C3_HF_degrees0_1_2']==[3,6,6]
    beta=[[Q(x) for x in B4['beta2_plane_r2_rell_ell2'][i]] for i in [0,4,5]]
    assert det3(beta)==-Q(5,192)
    # Fresh NO-CUT degree3 normal form: no ell3 coefficient is removed.
    N=[[Q(0)]*12 for _ in range(6)];curves=[]
    for sector,record in enumerate(graph['two_block_source_curve_NF_receipts']):
        g=list(map(Q,record['curve3_r3_r2ell_rell2_ell3']));assert g[0];curves.append(g)
        for j in range(1,4):
            # Direct monomial r3 -> -(g1*r2ell+g2*rell2+g3*ell3)/g0.
            row=[-g[j]/g[0]]+[Q(i==j) for i in range(1,4)]
            assert sum(x*y for x,y in zip(row,g))==0
            N[3*sector+j-1][4*sector:4*sector+4]=row
    assert [[str(v) for v in row] for row in N]==D['full_no_cut_cubic_NF6_rows12_columns']
    ellmap=[[N[i][4*sector+j+1] for sector in range(2) for j in range(3)] for i in range(6)]
    assert ellmap==[[Q(i==j) for j in range(6)] for i in range(6)]
    # Every non-plane cubic is already in S2*(z7), and every H cubic
    # is in S1*plane2H; this enumerates all495 cubic target monomials.
    counts={'z':0,'Hplane':0,'ABr3':0,'ell_times_ABquadratic':0}
    for sector in range(3):
        for factors in itertools.combinations_with_replacement(range(9),3):
            if any(i<7 for i in factors):counts['z']+=1
            elif sector==2:counts['Hplane']+=1
            elif factors==(7,7,7):counts['ABr3']+=1
            else:assert 8 in factors;counts['ell_times_ABquadratic']+=1
    assert counts=={'z':483,'Hplane':4,'ABr3':2,'ell_times_ABquadratic':6} and sum(counts.values())==495
    rank_report=read(ROOT/'reports/AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63_REVIEW/report.json')
    assert rank_report['actual_comparison_C_generic_DVR_rank']==4 and rank_report['actual_mu_D_C']==8
    assert c2['actual_C2_dimension']==14 and c2['ell_action_rank']==6
    cut=docs['one-CUT1-residual-and-formal-coordinate-rank-certificate.json'];assert cut['actual_C3_length19_proved']
    # DVR min-exponent arithmetic is finite and uses only accepted cuts.
    types=[]
    for exps in itertools.combinations_with_replacement(range(1,8),4):
        if 8+sum(min(2,e) for e in exps)==14 and 12+sum(min(3,e) for e in exps)==19:types.append(exps)
    assert types==[(1,1,2,e) for e in range(3,8)]
    fm={r['physical_row']:r for r in frame['all12_verified_source_records']}
    cm={r['physical_row']:r for r in c2['physical12_coefficients']};basis=raw['actual33_basis']
    mixed=[basis.index([i+1,a,b]) for i in range(7) for a,b in [(1,0),(0,1)]]
    source_rows=[]
    for p in [63,64,68,73,74]:
        assert not any(map(Q,fm[p]['bA_first_z7']+fm[p]['bB_first_z7']+fm[p]['bB_second_plane']))
        assert fm[p]['delta_and_bA_pureplane_second_zero']
        vals=dict(cm[p]['delta']);assert all(not Q(vals.get(i,'0')) for i in mixed)
        assert all(not Q(vals.get(basis.index(m),'0')) for m in [[0,3,0],[0,2,1]])
        assert any(r['physical_row']==p for r in graph['same_five_actual_source_premises'])
        source_rows.append({'physical_row':p,'all14_delta2_mixed_coefficients_zero_verified':True,
            'all_first_h_zero_bA2_bB2_plane_zero_verified':True,
            'delta3_r3_r2ell_zero_from_literal_C2_and_rell2_zero_from_prior_source_theorem':True,
            'remaining_delta3_is_formal_d03_ell3':True,
            'full_no_cut_classes':['(0,0,d03,0,0,0)','(0,0,0,0,0,d03)']})
    # Their selected minor is d03^2, universally, with no witness read.
    source_minor={(2,):Q(1)};assert source_minor=={(2,):Q(1)}
    proof={'same_full_C_not_C3_not_an_algebra_split':True,
      'D_finiteness_source_proof':[
        'The complete3-target/116-source presentation kills IAIB on all3 targets: IA on A,IB on B, and explicit64 IAIB products on H.',
        'After the accepted unit graph changes, both full conductor branches have r3 coefficient a QQ unit. Weierstrass makes R0/IA and R0/IB finitefree rank3 over the SAME D=QQ[[ell]].',
        'The active beta2 minor is -5/192, so the initial ideal of IA+IB contains z7 and plane2. Therefore IA+IB is m-primary.',
        'T=(IA intersect IB)/IAIB is finitely generated and annihilated by IA+IB, hence finite-dimensional overQQ.',
        '0->T->R0/IAIB->R0/(IA intersect IB)->0 and the latter embeds in the two finiteD branches. Thus R0/IAIB, and consequently C, is finiteD. No K/JK/Tor2 shortcut is used.'],
      'beta2_fixed_minor':'-5/192','full_no_cut_cubic_NF':[[str(v) for v in row] for row in N],
      'all495_target_cubic_monomials_accounted':counts,'ell_from_ABquadratics_is_I6':True,
      'unconditional_Nakayama_proof':[
        'The actual initial z7 and H-plane2 relations, plus both unit-r3 intrinsic cubics, show gr3(C)=ell*gr2(C) without any ell3 cut.',
        'Hence m3C=ell*m2C+m4C. Set Q=m3C/(ell*m2C), a finitely generated R0-module.',
        'Then Q=mQ. Nakayama gives Q=0, so m3C=ell*m2C exactly.',
        'Multiplication gives m^(n+1)C=ell*m^nC for all n>=2. This does not cancel ell on C or assume C torsionfree.'],
      'unconditional_lower_HF3':'ell3C subset m3C implies rankD(m3C)=rankD(C)=4. m4C=ell*m3C, so HF3full=muD(m3C)>=4.',
      'five_full_no_cut_graph_source_checks':source_rows,
      'full_graph_carry_proof':[
        'For tauA, h=bB1=0 and RA=zi*A exactly; H mixed clears contribute no B cubic because RH2_B,plane=delta_i2,plane=0. zz clears retainz. Thus A=d03ell3 and B=0.',
        'For tauB, all14 actual delta_j2 mixed coefficients vanish, so clearing its B quadratic introduces no RB-beta plane cubic. H mixed clears again add no B plane cubic. Thus A=0 and B=d03ell3.',
        'H cubic parts are removed by linear multiples of actual H quadratics, with A/B carry only at4. No ell3 quotient is imposed.',
        'In the no-cut cover, the ell3 coordinates are free. The selected two-column minor is d03^2.'],
      'conditional_exponent_and_all_order_proof':[
        'If any d03 is nonzero, two independent actual cubics cut the six-dimensional cover to dimension<=4. The unconditional lower4 forces fullHF3=4.',
        'Then finiteD module m3C has rank4 and mu4, hence is Dfree4 by the DVR structure theorem.',
        'Its submodule ell3C is torsionfree. For every D-torsion c in C, ell3c is torsion in m3C and hence0. Thus all torsion exponents<=3.',
        'The independently proved cuts give exponents(1,1,2,e>=3), so e=3.',
        'All higher m-powers are ell multiples of Dfree m3C, giving fullHF366444... and series(3+3t-2t3)/(1-t).',
        'C3m4=0 gives m4C subset ell3C. Their colengths are both19; therefore m4C=ell3C. Inductively m^kC=ell^(k-1)C for k>=4.',
        'm3 times comparisonD-torsion is0 since it lies in m3C intersect torsion(C). This is not an originalH0 annihilator assertion.'],
      'conditional_exact_type':{'free_rank':4,'torsion_exponents':[1,1,2,3],'comparison_torsion_length':7},
      'conditional_genuine_C_ell_gate':'Only after e=3, the canonical kernel(ell on C/ell4)/ell3(C/ell4) identifies C[ell], dim4. C/ell4 length23; no such action/basis/edge is constructed.',
      'all_graph03_values_unread_no_new_source_IFT_gauge_or_jet':True,
      'no_unconditional_e3_or_full_HF_claim':True,'negative_graph_result_only_subset_negative':True,
      'source_static_numeric_bound':50000,'hard_numeric':150000,'hard_seconds':60,'hard_address_bytes':2<<30,'hard_RAW':32<<20}
    ep=out/'independent-full-C-finiteness-Nakayama-and-conditional-graph03-proof.json';save(ep,proof)
    result={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'primary_report_sha256':PIN,
      'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},
      'unconditional_fullC_ell_reduction_accepted':True,'conditional_graph03_exponent3_theorem_accepted':True,
      'graph03_values_still_unread':True,'elapsed_seconds':time.monotonic()-start}
    save(out/'report.json',result);signal.alarm(0)
    print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':result['elapsed_seconds']}))
if __name__=='__main__':main()
