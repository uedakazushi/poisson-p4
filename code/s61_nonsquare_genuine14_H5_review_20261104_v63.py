#!/usr/bin/env python3
"""Independent70-entry genuine14/H5 binding, fixed units, no Gaussian.

Only sparse QQ[a,e,d] polynomial products are used; all rational checkpoint
denominators in these70 entries are constants. No POWER node is evaluated.
"""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
AUTHOR='S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63'
ID=AUTHOR+'_REVIEW'
ZERO=(0,0,0)
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def decode(items):
    out={tuple(m):Q(c) for m,c in items}
    assert len(out)==len(items) and all(len(m)==3 and min(m)>=0 and c for m,c in out.items())
    return out
def add(a,b):
    out=dict(a)
    for m,c in b.items():
        out[m]=out.get(m,Q(0))+c
        if not out[m]:del out[m]
    return out
def scale(a,c):return {m:x*c for m,x in a.items() if x*c}
MAXTERMS=0
def mul(a,b):
    global MAXTERMS
    out={}
    for m,x in a.items():
        for n,y in b.items():
            k=tuple(u+v for u,v in zip(m,n));out[k]=out.get(k,Q(0))+x*y
    out={m:x for m,x in out.items() if x};MAXTERMS=max(MAXTERMS,len(out))
    assert len(out)<=13854
    return out
def power(a,n):
    out={ZERO:Q(1)}
    for _ in range(n):out=mul(out,a)
    return out
def unit_product(a,b,divide=False):
    return {'constant':str(Q(a['constant'])/Q(b['constant']) if divide else Q(a['constant'])*Q(b['constant'])),
      'a_d_F_B_C_J_powers':[x-y if divide else x+y for x,y in zip(a['a_d_F_B_C_J_powers'],b['a_d_F_B_C_J_powers'])]}

def main(out):
    began=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/AUTHOR;report=read(folder/'report.json')
    assert report['all70_literal_heavyroots_bind_unchanged_old_checkpoint'] and report['all70_exact_knownunit_scaling_identities_passed']
    hashes={}
    for p,h in report['input_sha256'].items():assert sha(ROOT/p)==h;hashes[p]=h
    certificate=folder/'certificate.json';assert sha(certificate)==report['certificate_sha256']
    hashes[str(certificate.relative_to(ROOT))]=sha(certificate)
    hashes[str((folder/'report.json').relative_to(ROOT))]=sha(folder/'report.json')
    gates=[
      'S61_NONSQUARE_CUBIC_LAST_Q9_UNIT_COMPLETION_20261104_v63_REVIEW',
      'S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63_REVIEW',
      'S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63_REVIEW',
      'S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63_REVIEW',
      'S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63_REVIEW',
      'S61_NONSQUARE_RAW_H_MINOR_OPEN_20261104_v63_REVIEW']
    for name in gates:
        p=ROOT/'reports'/name/'report.json';g=read(p);assert g['all_passed']
        hashes[str(p.relative_to(ROOT))]=sha(p)
    cert=read(certificate)
    final=read(ROOT/'reports/S61_NONSQUARE_CUBIC_LAST_Q9_UNIT_COMPLETION_20261104_v63/certificate.json')
    parent=read(ROOT/'reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63/certificate.json')
    old=read(ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json')
    raw=read(ROOT/'reports/S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63/checkpoint.json')
    pre=read(ROOT/'reports/S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63/preflight.json')
    data=read(ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63/input.json')
    rows=final['new14_source_labels'];columns=cert['final5_heavy_target_labels']
    assert rows==cert['final14_source_labels']==[10,12,13,14,15,16,17,18,19,20,21,22,23,24]
    assert columns==[27,28,29,30,31]
    positions=[final['new21_target_labels'].index(i) for i in columns]
    assert positions==cert['final_target_positions']==[16,17,18,19,20]
    assert final['scalar_DAG']['nodes'][0]==parent['scalar_DAG']['nodes'][0]==['z']
    traces=[];literal=[]
    for i,row in enumerate(rows):
        this=[]
        for j,pos in zip(columns,positions):
            ref=final['actual14x21_scalar_roots'][i][pos]
            if ref==0:
                assert not old['matrix25x32'][row][j][0];trace=[0]
            else:
                op,ref2=final['scalar_DAG']['nodes'][ref]
                assert op=='matrix_parent'
                leaf=parent['scalar_DAG']['nodes'][ref2]
                assert leaf==['checkpoint',row,j]
                trace=[ref,ref2,leaf]
            traces.append([row,j,trace]);this.append(old['matrix25x32'][row][j])
        literal.append(this)
    assert traces==cert['all70_literal_external_root_traces']
    assert literal==cert['old_checkpoint70_rationals']
    # Every path inspected above has at most two direct external hops. It
    # therefore has no POWER ancestor; no pure target tail is expanded.
    a={(1,0,0):Q(1)};e={(0,1,0):Q(1)};d={(0,0,1):Q(1)};one={ZERO:Q(1)}
    F=add(add(scale(a,500),scale(power(e,2),90)),add(scale(e,3),scale(one,-24)))
    C=add(add(scale(a,100),scale(power(e,2),18)),add(scale(e,-9),scale(one,-2)))
    B=add(mul(add(add(scale(a,50),scale(power(e,2),18)),add(scale(e,-9),scale(one,-2))),d),add(scale(e,60),scale(one,-40)))
    J=add(scale(mul(a,d),5),add(scale(e,-6),scale(one,-4)))
    factors=[a,d,F,B,C,J]
    def unit_pair(u):
        c=Q(u['constant']);ex=u['a_d_F_B_C_J_powers']
        assert c and len(ex)==6 and all(isinstance(x,int) for x in ex)
        numerator={ZERO:c};denominator=one
        for f,k in zip(factors,ex):
            if k>=0:numerator=mul(numerator,power(f,k))
            else:denominator=mul(denominator,power(f,-k))
        return numerator,denominator
    sunit={'constant':'1/12500000','a_d_F_B_C_J_powers':[0,0,1,2,1,0]}
    gamma={'constant':'-4/625','a_d_F_B_C_J_powers':[0,0,1,1,2,0]}
    s_num,s_den=unit_pair(sunit)
    assert s_den==one and decode(pre['matrix15x6'][1][1])==s_num
    assert cert['known_column_units']==['1','1','1','1',gamma]
    gamma_poly,gamma_den=unit_pair(gamma);assert gamma_den==one
    keep=pre['remaining15_genuine_indices']
    assert keep==[0,1,2,3,4,5,7,8,9,10,11,12,13,14,15]
    step,=raw['steps'];assert step['pivot_row']==step['pivot_column']==1
    eliminated={q['row']:q for q in step['eliminations']}
    assert raw['remaining_rows']==[i for i in range(15) if i!=1]
    assert raw['remaining_columns']==[0,2,3,4,5]
    source_map=[];compact=[];numterms=0;nonzero=0
    for i,row in enumerate(rows):
        label=data['source25_labels'][row];assert label[0]=='ZW_genuine'
        genuine=label[1];rawindex=keep.index(genuine)
        assert rawindex in raw['remaining_rows'];source_map.append([row,genuine,rawindex])
        init=raw['initial_row_divisors']
        if rawindex in eliminated:
            tau=unit_product(unit_product(init[1],init[rawindex]),eliminated[rawindex]['row_divisor_after'])
        else:
            assert not pre['matrix15x6'][rawindex][1]
            tau=unit_product(sunit,init[rawindex])
        rho=cert['known_row_units_rho14'][i]
        # Recover the claimed original clearing unit from rho*s/tau and check
        # its complete source polynomial. No factor search/division scan.
        clearing_unit=unit_product(unit_product(rho,sunit),tau,divide=True)
        cn,cd=unit_pair(clearing_unit)
        clearing=decode(data['genuine16_source_clearing_units'][genuine])
        assert mul(clearing,cd)==cn
        rhon,rhod=unit_pair(rho)
        assert rhod==one
        assert tuple(rho['a_d_F_B_C_J_powers']) in ((0,0,0,0,1,0),(1,0,1,0,2,0),(0,0,1,0,2,0))
        compactrow=[]
        for k,column in enumerate(raw['remaining_columns']):
            p=decode(raw['matrix15x6'][rawindex][column])
            assert p==decode(cert['raw_compact14x5'][i][k])
            numerator=decode(literal[i][k][0]);denominator=decode(literal[i][k][1])
            assert set(denominator)=={ZERO} and denominator[ZERO]
            expected=mul(rhon,p)
            if k==4:expected=mul(expected,gamma_poly)
            assert numerator==mul(denominator,expected),(row,columns[k])
            numterms+=len(numerator);nonzero+=bool(numerator);compactrow.append(p)
        compact.append(compactrow)
    assert source_map==cert['actual_genuine_raw_checkpoint_index_map']
    assert numterms==13854 and nonzero==43
    assert sum(len(p) for r in compact for p in r)==2407
    assert max(len(p) for r in compact for p in r)==148
    colterms=[sum(len(row[j]) for row in compact) for j in range(5)]
    assert colterms==[415,304,302,302,1084]
    zero_rows=[r for r,row in zip(rows,compact) if not any(row)]
    assert zero_rows==[10,13,16]
    # The L-coordinate cancellation is a universal identity: after multiplying
    # by Lx1, its two correction terms have equal opposite monomials.
    # Symbols stand for (s,Xih,Xi1,X1h,Lx1,Lh).
    def mon(*v):return {tuple(sorted(v)):Q(1)}
    left={}
    for coefficient,word in ((1,(0,4,1)),(-1,(0,5,2)),(-1,(2,4,3)),(1,(2,5,0))):
        for m,v in mon(*word).items():left[m]=left.get(m,Q(0))+coefficient*v
    left={m:v for m,v in left.items() if v}
    assert left=={tuple(sorted((4,0,1))):Q(1),tuple(sorted((4,2,3))):Q(-1)}
    proof=dict(all70_literal_final_parent_checkpoint_bindings=True,
      no_POWER_ancestor_or_pure_degree207_expansion=True,
      all14_source_clearing_polynomials_and_fixed_tau_s_rho_identities=True,
      all70_QQ_polynomial_cross_multiplications_exact=True,
      all_checkpoint_denominators_constant=True,old_numerator_terms=numterms,
      compact_terms=2407,compact_nonzeros=43,compact_max_entry_terms=148,
      compact_column_term_counts=colterms,zero_heavy_source_rows=zero_rows,
      rho14_and_gamma_units_only_on_FULLU_adFBCJ=True,
      gamma='-4*F*B*C^2/625',universal_L_correction_cancellation=True,
      generic_rank5_and_open_D_pQ_only_by_independent_raw_rank_gates=True,
      p_and_Q_not_global_units=True,all_pure_tails_retained_lazily=True,
      no_Gaussian_scan_Bblock_minor_GB_newpoint_or_rank_probe=True,
      no_full_H7_full57_quartic239_or_exhaustion_claim=True,
      maximum_polynomial_intermediate_terms=MAXTERMS)
    (out/'independent-70-binding-and-unit-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in (Path(__file__),ROOT/f'research/notes/{ID}.md'):hashes[str(p.relative_to(ROOT))]=sha(p)
    review=dict(id=ID,all_passed=True,input_sha256=hashes,
      evidence_sha256={'independent-70-binding-and-unit-proof.json':sha(out/'independent-70-binding-and-unit-proof.json')},
      all70_heavy_bindings_and_knownunit_identities_independently_accepted=True,
      no_global_heavy_rank5_or_full57_claim=True,
      seconds=time.monotonic()-began,max_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(review,indent=2)+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(review))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
