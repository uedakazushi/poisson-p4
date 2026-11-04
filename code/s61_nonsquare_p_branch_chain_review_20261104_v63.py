#!/usr/bin/env python3
"""Independent reversed-variable Fraction Laurent audit of five frozen packets."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_RAW_H_P_BRANCH_CHAIN_20261104_v63_REVIEW'
NAMES=['RAW_H_P_BRANCH_PROFILES','RAW_H_G1_UNIT_PAIR','RAW_H_THETA_UNIT_PROFILE','RAW_H_THETA_ONE_DEFECT','THETA_RATIONAL_COORDINATES']
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def dec(a):
    p={}
    for m,c in a:
        m=tuple(reversed(m));p[m]=p.get(m,F(0))+F(c)
    return {m:c for m,c in p.items() if c}
def add(a,b,c=F(1)):
    p=dict(a)
    for m,v in b.items():
        q=p.get(m,F(0))+c*v
        if q:p[m]=q
        else:p.pop(m,None)
    return p
def sc(a,c):return {m:c*v for m,v in a.items() if c*v}
def mul(a,b):
    p={}
    for m,c in a.items():
        for n,d in b.items():
            k=tuple(x+y for x,y in zip(m,n));p[k]=p.get(k,F(0))+c*d
    return {m:c for m,c in p.items() if c}
def pw(a,k,n=3):
    p={(0,)*n:F(1)}
    for i in range(k):p=mul(p,a)
    return p
def shift(a,m,c=F(1)):return {tuple(x+y for x,y in zip(k,m)):c*v for k,v in a.items() if c*v}
def det3(A):
    p={}
    for perm in itertools.permutations(range(3)):
        q={(0,0,0):F(1)}
        for i,j in enumerate(perm):q=mul(q,A[i][j])
        p=add(p,q,F((-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))))
    return p
def monic_reduce(p,H,var,degree):
    q={};r=dict(p);step=0
    while r:
        bad=[m for m in r if m[var]>=degree]
        if not bad:break
        m=max(bad,key=lambda a:(a[var],a));e=list(m);e[var]-=degree
        term={tuple(e):r[m]};q=add(q,term);r=add(r,mul(term,H),F(-1));step+=1
        assert step<10000 and len(r)<100000
    assert add(mul(q,H),r)==p
    return q,r
def substitute(p,images,n):
    out={};cache={}
    for m,c in p.items():
        term={(0,)*n:c}
        for j,k in enumerate(m):
            if not k:continue
            if k<0:
                assert j==0 and images[j]=={(1,)+(0,)*(n-1):F(1)}
                factor={(k,)+(0,)*(n-1):F(1)}
            else:
                if (j,k) not in cache:cache[j,k]=pw(images[j],k,n)
                factor=cache[j,k]
            term=mul(term,factor)
        out=add(out,term)
    return out
def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    folders=[ROOT/'reports'/f'S61_NONSQUARE_{s}_20261104_v63' for s in NAMES]
    reports=[read(p/'report.json') for p in folders];cert=[read(p/'certificate.json') for p in folders]
    checked={}
    for folder,r in zip(folders,reports):
        assert sha(folder/'certificate.json')==r['certificate_sha256'] and r['frozen']
        for p,h in r['input_sha256'].items():
            path=ROOT/p;assert sha(path)==h,p;checked[p]=h
    P,G,T,O,R=cert
    rawp=ROOT/'reports/S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63/checkpoint.json'
    genp=ROOT/'reports/S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63/certificate.json'
    bp=ROOT/'reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63/certificate.json'
    br=ROOT/'reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63_REVIEW/report.json'
    assert read(br)['all_passed'];raw,generic,bridge=map(read,[rawp,genp,bp])
    M=[[dec(p) for p in row] for row in raw['matrix15x6']];remaining=raw['remaining_rows'];heavy=[0,2,3,4,5]
    assert remaining==generic['input14_remaining_rows']==[0,2,3,4,5,6,7,8,9,10,11,12,13,14]
    assert [[M[i][j] for j in heavy] for i in remaining]==[[dec(p) for p in row] for row in bridge['raw_compact14x5']]
    for q in bridge['known_row_units_rho14']:
        assert F(q['constant']) and len(q['a_d_F_B_C_J_powers'])==6 and all(isinstance(x,int) for x in q['a_d_F_B_C_J_powers'])
    for q in bridge['known_column_units']:
        if isinstance(q,str):assert F(q)
        else:assert F(q['constant']) and len(q['a_d_F_B_C_J_powers'])==6
    one={(0,0,0):F(1)};d={(1,0,0):F(1)};e={(0,1,0):F(1)};a={(0,0,1):F(1)}
    C=add(add(sc(a,100),sc(pw(e,2),18)),add(sc(e,-9),sc(one,-2)))
    B=add(mul(add(sc(a,50),add(sc(pw(e,2),18),add(sc(e,-9),sc(one,-2)))),d),add(sc(e,60),sc(one,-40)))
    FF=add(sc(a,500),add(sc(pw(e,2),90),add(sc(e,3),sc(one,-24))))
    JJ=add(sc(mul(a,d),5),add(sc(e,-6),sc(one,-4)))
    g1=add(sc(mul(e,d),3),add(sc(d,-2),sc(one,20)))
    g2=add(sc(mul(e,d),6),add(d,sc(one,40)))
    theta=add(pw(B,2),mul(C,mul(g1,g2)),F(-1))
    assert dec(P['g1'])==g1 and dec(P['g2'])==g2 and dec(P['Theta'])==theta
    assert M[2][0]==M[4][0]==dec(P['source_p'])==dec(generic['p'])
    assert sc(M[2][0],750000)==mul(g1,theta)
    assert add(theta,mul(C,mul(g1,g2)))==pw(B,2)
    assert {m:c for m,c in theta.items() if m[2]==2}=={(2,0,2):F(2500)}
    profiles={};count=0
    for branch in P['branches']:
        name=branch['branch'];H=g1 if name=='g1' else theta;var=1 if name=='g1' else 2;degree=1 if name=='g1' else 2
        rm={}
        for row in branch['full14x2_records']:
            i,j=row['raw_checkpoint_row'],row['column'];q,r=dec(row['quotient_by_original_branch_polynomial']),dec(row['remainder'])
            assert add(mul(q,H),r)==M[i][j]
            assert all(m[var]<degree for m in r);rm[i,j]=r;count+=1
        assert len(rm)==28
        assert not any(M[i][j] for i in branch['sparse_other3_zero_rows'] for j in [2,3,4])
        profiles[name]=rm
    assert count==56
    TT=add(mul(a,pw(d,2)),add(sc(d,-2),sc(one,8)))
    expected=[[shift(pw(TT,3),(-3,0,0),F(1,100)),shift(pw(TT,4),(-4,0,0),F(1,10))],
      [sc(pw(TT,2),F(1,300)),shift(pw(TT,3),(-1,0,0),F(1,6))]]
    assert [[dec(v) for v in row] for row in G['pair_matrix']]==expected
    for rec in G['selected_full_source_quotient_identities']:
        i,j=rec['row'],rec['column'];assert dec(rec['original'])==M[i][j]
        assert add(mul(dec(rec['quotient_by_g1']),g1),dec(rec['remainder']))==M[i][j]
    assert all(not M[i][j] for i in [5,8] for j in [2,3,4])
    bq=add(sc(e,6),add(one,shift(one,(-1,0,0),F(-20))))
    assert add(B,shift(TT,(-1,0,0),F(-50)))==mul(g1,bq)==mul(g1,dec(G['B_minus_50T_over_d_quotient_by_g1']))
    determinant=add(mul(expected[0][0],expected[1][1]),mul(expected[0][1],expected[1][0]),F(-1))
    assert determinant==dec(G['determinant'])==shift(pw(TT,6),(-4,0,0),F(1,750))
    unitdet=shift(pw(B,6),(2,0,0),F(1,750*50**6))
    assert unitdet==dec(G['unit_determinant_representative'])
    assert add(determinant,unitdet,F(-1))==mul(g1,dec(G['det_minus_unitrepresentative_quotient_by_g1']))
    adj=[[expected[1][1],sc(expected[0][1],-1)],[sc(expected[1][0],-1),expected[0][0]]]
    for A0,A1 in [(expected,adj),(adj,expected)]:
        assert all(add(mul(A0[i][0],A1[0][j]),mul(A0[i][1],A1[1][j]))==(determinant if i==j else {}) for i in range(2) for j in range(2))
    pivot=sc(mul(g1,pw(g2,2)),F(1,37500));n=sc(mul(pw(C,2),g1),37500);B4=pw(B,4)
    assert pivot==dec(T['pivot_remainder']) and n==dec(T['inverse_numerator']) and B4==dec(T['inverse_denominator'])
    assert M[8][0]==add(pivot,mul(theta,dec(T['pivot_original_quotient_by_Theta'])))
    assert add(mul(pivot,n),B4,F(-1))==mul(theta,add(theta,pw(B,2),F(-2)))
    assert add(mul(M[8][0],n),B4,F(-1))==mul(theta,dec(T['full_pivot_times_inverse_minus_one_cleared_quotient_by_Theta']))
    assert not any(M[8][j] for j in [2,3,4])
    bRows=[7,9,10,11,12,13,14];assert T['full_source_rows']==bRows==generic['all_other3_column_support_exact_rows']
    blocks=[[dec(p) for p in row] for row in T['fixed_7x4_matrix']]
    actualB=[[M[i][j] for j in [2,3,4]] for i in bRows]
    assert actualB==[[dec(p) for p in row] for row in generic['actual_block7x3']]
    for rec in T['all28_full_original_quotient_identities']:
        i,j=rec['raw_checkpoint_row'],rec['new_column']
        value=M[i][j] if j in [2,3,4] else add(mul(B4,M[i][5]),mul(mul(M[i][0],n),M[8][5]),F(-1))
        rem=dec(rec['remainder']);assert add(mul(theta,dec(rec['quotient_by_Theta'])),rem)==value
        jj=[2,3,4,5].index(j)
        assert rem==blocks[bRows.index(i)][jj] and all(m[2]<2 for m in rem)
    for i in T['zero_rows_after_specialization']:
        assert not profiles['Theta'][i,0] and not profiles['Theta'][i,5] and not any(M[i][j] for j in [2,3,4])
    Delta=dec(generic['nonzero_selected_minor_Delta']);selected=generic['selected_rows7'];Lamb=[[dec(p) for p in row] for row in generic['cleared_upper_row_coordinates7xr']]
    assert Delta and Delta==det3([[actualB[i][j] for j in generic['selected_cols3']] for i in selected])
    for i in range(7):
        for j in range(3):
            rhs={}
            for k,b in enumerate(selected):rhs=add(rhs,mul(Lamb[i][k],actualB[b][j]))
            assert mul(Delta,actualB[i][j])==rhs
    U,V=actualB[3][0],actualB[0][0];assert U==dec(O['U']) and V==dec(O['V'])
    assert all(not add(mul(V,actualB[3][j]),mul(U,actualB[0][j]),F(-1)) for j in range(3))
    assert dec(O['Delta'])==Delta and dec(O['old_lambda_first'])==Lamb[3][0]
    assert not any(Lamb[3][j] for j in [1,2])
    assert mul(Delta,U)==mul(Lamb[3][0],V)
    monic=shift(theta,(-2,0,0),F(1,2500))
    _,u=monic_reduce(U,monic,2,2);_,v=monic_reduce(V,monic,2,2)
    assert u==dec(O['Theta_reduced_U']) and v==dec(O['Theta_reduced_V'])
    for j,q in enumerate(O['all3_middle_cleared_relations_Theta_quotients']):
        assert add(mul(v,blocks[3][j]),mul(u,blocks[0][j]),F(-1))==mul(theta,dec(q))
    defect=add(mul(v,blocks[3][3]),mul(u,blocks[0][3]),F(-1))
    assert defect==dec(O['defect_before_final_Theta_reduction'])==mul(theta,dec(O['quotient_by_Theta'])) and not dec(O['defect_remainder'])
    NN=add(mul(d,B),mul(g1,g2),F(-1));assert NN==dec(R['squarecompletion']['N'])
    assert add(pw(NN,2),sc(mul(g1,pw(g2,2)),20),F(-1))==mul(pw(d,2),theta)
    # Target order is (d,z), the reverse of the producer's order.
    O2={(0,0):F(1)};D={(1,0):F(1)};z={(0,1):F(1)};zp=add(z,sc(O2,20))
    LL=add(pw(z,2,2),sc(D,50));KK=add(mul(zp,add(z,sc(O2,-10))),sc(D,50))
    AA=shift(mul(zp,KK),(-2,0),F(1,500));EE=add(sc(O2,F(2,3)),shift(add(pw(z,2,2),sc(O2,-400)),(-1,0),F(1,60)))
    BB=shift(mul(mul(LL,z),zp),(-1,0),F(1,200));CC=shift(mul(LL,pw(zp,2,2)),(-2,0),F(1,200))
    tc={'L':LL,'K':KK,'a':AA,'e':EE,'B':BB,'C':CC,'g1':sc(pw(z,2,2),F(1,20)),'g2':sc(LL,F(1,10))}
    Fnum=add(mul(pw(z,2,2),pw(zp,2,2)),add(mul(D,add(sc(pw(z,2,2),82),add(sc(z,2000),sc(O2,7200)))),sc(pw(D,2,2),720)))
    Jnum=add(mul(pw(z,2,2),zp),sc(mul(D,add(z,sc(O2,4))),50));tc.update(Fnum=Fnum,Jnum=Jnum)
    assert all(p==dec(R['target_coordinates'][k]) for k,p in tc.items())
    source={'a':a,'e':e,'d':d,'g1':g1,'g2':g2,'B':B,'C':C,'F':FF,'J':JJ,'Theta':theta,'N':NN}
    mapped={k:substitute(p,[D,EE,AA],2) for k,p in source.items()}
    assert all(p==dec(R['forward_all_images'][k]) for k,p in mapped.items())
    assert not mapped['Theta'] and mapped['N']==mul(z,tc['g2'])
    assert shift(mapped['F'],(2,0),40)==Fnum and shift(mapped['J'],(1,0),100)==Jnum
    zrep=R['source_Z_regular_representative'];assert dec(zrep['numerator'])==mul(mul(NN,C),g1) and dec(zrep['denominator'])==pw(B,2)
    assert add(mul(dec(zrep['numerator']),g2),mul(NN,pw(B,2)),F(-1))==sc(mul(NN,theta),-1)
    def rational_sub_z(poly):
        degree=max([m[1] for m in poly]+[0]);num={}
        for (dp,zp0),c in poly.items():num=add(num,shift(mul(pw(NN,zp0),pw(g2,degree-zp0)),(dp,0,0),c))
        return num,pw(g2,degree)
    for rec in R['all_six_reverse_composition_quotients']:
        k=rec['coordinate'];num,den=rational_sub_z(tc[k]);an=dec(rec['substituted_numerator']);ad=dec(rec['substituted_denominator'])
        assert mul(num,ad)==mul(an,den)
        assert add(an,mul(ad,source[k]),F(-1))==mul(theta,dec(rec['difference_quotient_by_Theta']))
    # All target unit witnesses hold on the entire scheme; no reducedness premise.
    unit_checks=[add(mul(g1,mul(C,g2)),pw(B,2),F(-1)),
      add(pw(add(NN,sc(g2,20)),2),sc(mul(mul(pw(d,2),C),g2),20),F(-1))]
    for value in unit_checks:assert not monic_reduce(value,monic,2,2)[1]
    for name,scaleq,dp in [('Fnum',40,2),('Jnum',100,1)]:
        num,den=rational_sub_z(tc[name]);assert not monic_reduce(add(num,shift(mul(den,source['F' if name=='Fnum' else 'J']),(dp,0,0),F(scaleq)),F(-1)),monic,2,2)[1]
    evidence={'variable_order':'independent source(d,e,a), target(d,z)',
      'all_input_hashes_current':True,'raw70_checkpoint_bridge_coefficients_exact':True,'known_row_column_scalings_use_only_FULLU_units':True,
      'p_factor_and_comaximality':True,'profile_full_quotient_identities':56,
      'g1_fixed_pair_four_formulas_source_quotients_inverse_and_B_unit':True,
      'Theta_original_unit_row_and_full28_Schur_transport':True,
      'all21_generic_B_Cramer_identities_and_fixed_Delta_rebuilt':True,
      'one_fixed_raw11_defect_only_exact_zero':True,
      'Theta_two_way_coordinate_maps_and_all_FULLU_unit_witnesses':True,
      'Theta_target_Fnum_terms':len(Fnum),'Theta_target_Jnum_terms':len(Jnum),
      'whole_Theta_heavy_rank_or_full_cubic_rank_not_inferred_from_these5_packets':True}
    op=out/'independent-QQ-p-branch-chain.json';op.write_text(json.dumps(evidence,indent=2)+'\n')
    inputs=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',rawp,genp,bp,br]
    for p in folders:inputs += [p/'report.json',p/'certificate.json']
    report={'id':STEM,'all_passed':True,'independent_fraction_Laurent_all_five_packet_arithmetic':True,
      'source_checkpoint_and_actual_heavy_unit_transport_bound':True,'all56_plus28_original_quotient_identities_verified':True,
      'rational_Theta_scheme_isomorphism_and_all_localizers_accepted':True,
      'no_new_ranksearch_points_primes_GB_or_full_heavy_surjectivity_claim':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},'checked_primary_source_hashes':checked,
      'evidence_sha256':{op.name:sha(op)},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['all_passed','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
