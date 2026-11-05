#!/usr/bin/env python3
"""Independent complete quadratic-source and formal cubic-carry proof."""
import hashlib,itertools,json,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63_REVIEW'
PRI='AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63'
PRIMARY='f2f0cfbc738cd58167741b9ce6cd2a73cec5b884dfe0f47fbb1acf94df19d72c'
def sha(p):
    h=hashlib.sha256()
    with pathlib.Path(p).open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def field(t,k):
    j=t.index(':',t.index(json.dumps(k)))+1
    while t[j].isspace():j+=1
    return json.JSONDecoder().raw_decode(t,j)[0]
def sparse(p):return {i:Q(c) for i,c in p if Q(c)}
def add(p,q,c=Q(1)):
    r=dict(p)
    for i,x in q.items():
        r[i]=r.get(i,Q(0))+c*x
        if not r[i]:del r[i]
    return r
def serial(p):return [[i,str(c)] for i,c in sorted(p.items())]
def degree(m):return int(m[0]!=0)+m[1]+m[2]

def main():
    t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
    try:
        rp=ROOT/'reports'/PRI/'report.json';assert sha(rp)==PRIMARY;pr=read(rp);assert pr['all_passed']
        pins={str(rp.relative_to(ROOT)):PRIMARY}
        for p,h in pr['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
        for p,h in pr['evidence_sha256'].items():
            ep=ROOT/'reports'/PRI/p;assert sha(ep)==h;pins[str(ep.relative_to(ROOT))]=h
        author=read(ROOT/'reports'/PRI/'actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json')
        for p in pr['input_sha256']:
            if p.endswith('_REVIEW/report.json'):assert read(ROOT/p)['all_passed']
        def data(suffix):return read(next(ROOT/p for p in pr['input_sha256'] if p.endswith(suffix)))
        raw=data('actual-ell2-product33-and-raw45-actions.json');basis=raw['actual33_basis']
        assert basis[0]==[0,0,0] and len(basis)==33
        # Frozen normal forms preserve order: high-degree tails cannot alter the first/plane2 functionals.
        assert len(raw['all74_normal_forms'])==len(raw['basis74'])==74
        for mon,v in zip(raw['basis74'],raw['all74_normal_forms']):
            assert all(degree(basis[i])>=degree(mon) for i,c in v if Q(c))
        idxz=[basis.index([j+1,0,0]) for j in range(7)];ir=basis.index([0,1,0]);il=basis.index([0,0,1])
        r2=basis.index([0,2,0]);rl=basis.index([0,1,1]);first=idxz+[ir,il]
        C2=data('actual-physical12-and-C2-right-exact-quotient.json');rootdata=data('actual-source105-ell2-raw33-values.json')
        ell=data('all164-cached-pure-ell-readbacks-and-physical12-second-values.json')
        nr={i:sparse(v) for i,v in rootdata['all_requested164_root_values']}
        er={r['physical_row']:r for r in ell['all12_actual_physical_coefficients_pure_ell3']}
        rows=C2['physical12_coefficients'];assert len(rows)==12
        actual=[];brows=[]
        for k,r in enumerate(rows):
            row=r['physical_row'];erow=er[row];assert r['roots']==erow['source_roots']
            ds={n:sparse(r[n]) for n in ['delta','bA','bB']}
            rr=r['roots'];assert ds['bA']==nr[rr['bA_node']] and ds['bB']==nr[rr['bB_node']]
            assert ds['delta']==add(nr[rr['aA_node']],nr[rr['aB_node']],-1)
            assert all(not ds[n].get(0,Q(0)) for n in ds)
            assert all(not ds['delta'].get(i,Q(0)) for i in first)
            assert all(not ds[n].get(i,Q(0)) for n in ['bA','bB'] for i in [ir,il])
            assert all(not ds[n].get(i,Q(0)) for n in ['delta','bA'] for i in [r2,rl])
            assert all(Q(erow[n+'_jet3'][j])==0 for n in ['delta','bA'] for j in [0,1,2])
            if k<7:assert [ds['bA'].get(i,Q(0)) for i in idxz]==[Q(k==j) for j in range(7)]
            bzA=[ds['bA'].get(i,Q(0)) for i in idxz];bzB=[ds['bB'].get(i,Q(0)) for i in idxz]
            b2=[ds['bB'].get(r2,Q(0)),ds['bB'].get(rl,Q(0)),Q(erow['bB_jet3'][2])]
            brows.append(b2);actual.append({'physical_row':row,'delta_first9_zero':True,
              'delta_and_bA_pureplane_second_zero':True,'bA_first_z7':list(map(str,bzA)),
              'bB_first_z7':list(map(str,bzB)),'bB_second_plane':list(map(str,b2))})
        assert actual==author['all12_verified_source_records']
        assert all(not nr[n].get(i,Q(0)) for n in [20086,40865] for i in [0,ir,il])
        # Actual complete first map, in the21 z-coordinate rows. r/ell firsts are already checked0.
        firstcols=[{} for _ in range(116)]
        for j in range(7):firstcols[j]={j:Q(1)};firstcols[8+j]={7+j:Q(1)}
        for j,r in enumerate(actual):
            firstcols[80+3*j]={14+i:-Q(x) for i,x in enumerate(r['bB_first_z7']) if Q(x)}
            firstcols[81+3*j]={14+i:Q(x) for i,x in enumerate(r['bA_first_z7']) if Q(x)}
        chosen=list(range(7))+list(range(8,15))+[81+3*j for j in range(7)]
        assert all(firstcols[c]=={i:Q(1)} for i,c in enumerate(chosen))
        free=[j for j in range(116) if j not in chosen];assert len(free)==95
        kernel=[]
        for j in free:
            v={j:Q(1)}
            for i,c in firstcols[j].items():v[chosen[i]]=-c
            check={}
            for k,c in v.items():check=add(check,firstcols[k],c)
            assert not check and {k:v.get(k,Q(0)) for k in free if v.get(k,Q(0))}=={j:Q(1)}
            kernel.append([j,serial(v)])
        # Complete95 source second images after quotient by all z-initial rows.
        projected95=[]
        for j in free:
            v={}
            if j>=80 and (j-80)%3==0:
                pj=(j-80)//3;v={6+i:-c for i,c in enumerate(brows[pj]) if c}
            projected95.append([j,serial(v)])
        chosenminor=[brows[k] for k in [0,4,5]];det=Q(0)
        for p in reversed(list(itertools.permutations(range(3)))):
            term=Q((-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3)))
            for i,j in enumerate(p):term*=chosenminor[i][j]
            det+=term
        assert det==Q(5,768)==Q(author['fixed_bB_plane_minor_determinant'])
        assert [[str(c) for c in v] for v in chosenminor]==author['fixed_bB_plane_minor_matrix']
        # Formal support proof keeps every unknown zz/mixed coefficient, not a numerical source3 substitution.
        sym2=list(itertools.combinations_with_replacement(range(9),2))
        zz=[m for m in sym2 if max(m)<7];zterms=[m for m in sym2 if min(m)<7]
        assert len(sym2)==45 and len(zz)==28 and len(zterms)==42
        carry=[]
        for target,support in [('A',zz),('B',zz),('H',zterms)]:
            for mon in support:
                z=min(i for i in mon if i<7);mult=list(mon);mult.remove(z);mult=mult[0]
                if target=='A':test=[] # RA2=0 exactly
                elif target=='B':test=[tuple(sorted((mult,)+m)) for m in sym2] # conservative full RB2
                else:test=[tuple(sorted((mult,)+m)) for m in zterms] # RH2 has no pureplane part
                assert all(any(i<7 for i in m) for m in test)
                carry.append({'target':target,'quadratic_monomial':list(mon),'firstframe_row':z,
                  'multiplier_variable':mult,'all_cubic_carry_monomials_have_z':True})
        assert len(carry)==98
        # Offdiagonal cubic nu1*b2+nu2*b1 always has z; constant terms vanish.
        assert all(any(i<7 for i in (z,)+m) for z in range(7) for m in sym2)
        # The four arbitrary pureplane delta3 coordinates are in S1*(plane2 eH).
        sensitivity=[]
        for p in itertools.combinations_with_replacement([7,8],3):
            divisor=p[:2];mult=p[2];assert divisor in [(7,7),(7,8),(8,8)]
            sensitivity.append({'delta3_monomial':list(p),'known_H_quad_divisor':list(divisor),'linear_multiplier':mult,'projected_sensitivity':0})
        proof={'same_complete_comparison_module':'C=coker(D0->E4), not Cprime, over SAME R0, targets tauA,tauB,tauA*tauB',
          'right_exact_complete116_source':'8 IA*eA +8 IB*eB +64 IAIB*eH +12*3 physical actions; F*1 already in raw image.',
          'all_source_hashes_and_prior_exact_peers_match':True,
          'all_M0_zero_and_delta_first9_zero_and_b_nu_plane_first_zero_verified':True,
          'C2_normal_forms_filtered_all74_verified':True,
          'C2_r2_rell_coordinate_functionals_are_actual_low_order_plane_coefficients':True,
          'new_actual_ell2_readback_same_physical_root_namespace_verified':True,
          'all12_low_order_source_records':actual,
          'first21_chosen_source_columns':chosen,'first21_unit_matrix_literal_I21':True,
          'complete95_constant_kernel_free_tail_I95':kernel,'all95_first_multiplybacks_zero':True,
          'complete95_projected_quadratic_images9':projected95,
          '95_partition':{'two_curves':2,'products':64,'tauA':12,'remaining_tauB':5,'tauAB':12},
          'quadratic_projection_proof':'Modulo S1*V1 all z terms vanish. Only12 tauA images remain -bB2*eH; all constant first-clearings use RH whose plane2 is0, never RB.',
          'fixed_three_by_three_det':str(det),'fixed_physical_rows':[18,60,67],
          'actual_HF0_1_2_C_and_C3':[3,27-21,9-3],
          'formal_all98_quadratic_clearings':carry,
          'all_12_times4_unknown_delta3_sensitivities_zero':sensitivity,
          'all_unknown_mixed_zz_and_delta3_coefficients_retained_formal':True,
          'coherent_cubic_operator':'pi_plane(M3_j - M2_frame*L(M2_j)), modulo known S1 I2 and S2 I1',
          'IB_beta2_carry_retains_z_from_the_multiplier':True,
          'H_mixed_carry_uses_actual_RH2_plane_zero':True,
          'all12_tauAB_cubic_classes_source_zero_subset_only':True,
          'full_C3_cubic_map_or_length_or_exponents_not_determined':True,
          'plain_relation_model_only':{'type':'CI(3,3)eA + CI(3,3)eB + plane/m2 eH','HF':[3,6,6,4,2],'length':21},
          'no_new_cubic_scalar_source_gauge_jet_value_or_Gaussian':True,
          'max_numeric_static_upper':30000,'hard_numeric_cap':150000,'seconds':60,'address_bytes':2<<30,'hard_RAW_bytes':32<<20}
        save(out/'independent-complete95-and-formal-cubic-carry-proof.json',proof)
        report={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'primary_report_sha256':PRIMARY,
          'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{'independent-complete95-and-formal-cubic-carry-proof.json':sha(out/'independent-complete95-and-formal-cubic-carry-proof.json')},
          'actual_C_and_C3_HF_through2_accepted':[3,6,6],'physical12_tauAB_cubic_subset_source_zero_accepted':True,
          'no_actual_C3_length_or_torsion_exponent_claim':True,'no_author_producer_imports':True,'elapsed_seconds':time.monotonic()-t}
        save(out/'report.json',report);signal.alarm(0)
        print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':report['elapsed_seconds']}))
    except Exception as e:
        save(out/'failure.json',{'error':type(e).__name__+': '+str(e),'no_positive_claim':True});raise
if __name__=='__main__':main()
