#!/usr/bin/env python3
"""Actual N2/L2 source functionals from previously solved B11 caches only.

No physical source coefficient or gauge equation is solved. Full ordered
telescopes retain off-intersection terms before N*F=d*(G,0) is checked.
"""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_comparison_X_origin_20261104_v63 as H
from ambient_B18_native128_first_20261104_v63 import PATHS as OLD_PATHS
from ambient_B18_T_native_S2_B11_20261104_v63 import Finite,dense,counts,encode,QUAD,LINEAR
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63'
def main(out):
    st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    keys=['source','interface','proper','equations','first5','second5','membership','preparation','normalizer']
    paths={k:OLD_PATHS[k] for k in keys}
    paths.update({
      'B11':'reports/AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT_20261104_v63/actual-B11-source-and-gauge-preflight.json',
      'phase1':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase1-source-B11-values.json',
      'phase2':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase2-gauge30-B11-correction.json',
      'phase3':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase3-actual-proper-G8-B11.json',
      'cache_report':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/report.json',
      'cache_peer':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63_REVIEW/report.json',
      'first_source':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase4-native-source-and-R12.json',
      'first_adapter':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase5-native-adapter-firsts.json',
      'first_preflight':'reports/AMBIENT_B18_NATIVE128_FIRST_SOURCE_PREFLIGHT_20261104_v63/actual-native128-first-input-and-mixed-index-profile.json',
      'adjoint':'reports/AMBIENT_B18_T_NATIVE_L2_ADJOINT_PREFLIGHT_20261104_v63/all619-full-c1-native-L2-adjoints.json',
      'adjoint_report':'reports/AMBIENT_B18_T_NATIVE_L2_ADJOINT_PREFLIGHT_20261104_v63/report.json'})
    saved=0;livepeak=0
    def save(name,data,live=()):
        nonlocal saved,livepeak
        q=counts(data);lv=counts(live);livepeak=max(livepeak,saved+q+lv)
        assert livepeak<150000,('TOTAL_QQ_GUARD',name,saved,q,lv)
        p=out/name;p.write_text(json.dumps(encode(data),separators=(',',':'))+'\n');saved+=q
        assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
        (out/'phase.json').write_text(json.dumps({'phase':name,'saved_QQ_slots':saved,'max_saved_plus_live_QQ':livepeak,'seconds':time.monotonic()-st})+'\n')
    hashes={p:H.sha(ROOT/p) for p in paths.values()}
    save('source-manifest.json',{'input_sha256':hashes,'no_source_or_gauge_solve':True})
    D={k:H.read(ROOT/p) for k,p in paths.items()}
    assert D['cache_report']['all_passed'] and D['cache_peer']['all_passed'] and D['adjoint_report']['all_passed']
    S,I,P,Q,R,K,C,prep,E=[D[k] for k in keys]
    source=D['phase3'];actual=[tuple(map(F,p)) for p in source['actual14_C_B44']];kappa=[tuple(map(F,p)) for p in source['kappa9_C_B44']]
    F12=[tuple(map(F,p)) for p in source['physicalF12_C_B44']];gs=[tuple(map(F,p)) for p in source['actualG8_C_B44']]
    nu=[tuple(map(F,p)) for p in D['phase1']['nu_A_B_B11']]
    first={int(i):dense(v) for i,v in D['B11']['all105_actual_B11_root_values']}
    ev=[tuple(map(F,p)) for p in D['phase2']['even15_B11']];od=[tuple(map(F,p)) for p in D['phase2']['odd15_B11']]
    requested={tuple(k) for row in D['adjoint']['all619'] for k,c in row['actual_L2_source_adjoint']}
    assert all(len(p)==44 for p in actual+kappa+F12+gs)
    inputslots=counts([actual,kappa,F12,gs,nu,first,ev,od]);matrixslots=44*(17*12+17*5+12*5+5*5+4*12*12)
    bound=inputslots+matrixslots+2*44*12*12+len(requested)+12000
    assert bound<150000
    save('preflight.json',{'actual_cached_input_QQ_slots':inputslots,'matrix_and_scratch_forecast':matrixslots+2*44*12*12+12000,
      'whole_saved_live_forecast':bound,'cap':150000,'requested_L2_scalar_entries':len(requested),
      'all105_source_values_already_cached':True,'all30_gauge_values_already_cached':True,
      'new_physical_source_or_gauge_solve':False,'source_adapter':'ordered full17 first5 telescopes, then full12 core5 telescope'},live=[actual,kappa,F12,gs,first])
    a=Finite(2,nu);a.install();GL=D['first_preflight']['GL18_source_parameter_positions']
    def frame_family(family,tauindex):
        f=S['families'][family];ss=a.dual(first[f['p_quadratic_relation']['s']]);p=H.add(a.tau(tauindex),H.scale(ss,F(1,2)));G=H.ident(5)
        for name,(aa,bb) in f['coefficient_unknown_pairs'].items():
            if name.startswith('gamma') and name!='gamma0':
                index=int(name[5:]);i,j=divmod(GL[index-1],5);G[i][j]=H.add(G[i][j],a.dual(first[aa]),H.mul(a.dual(first[bb]),p))
        core=H.ident(5)
        if family=='CA':core[0][3]=H.scale(p,F(-2));core[1][0]=p;core[1][3]=H.neg(H.mul(p,p))
        else:core[0][3]=H.scale(p,F(1,2));core[1][0]=p;core[1][3]=H.scale(H.mul(p,p),F(1,4))
        return H.mm(G,core)
    G=H.ident(5)
    for (i,j),x,y in zip(I['actual_group15_positions'],ev,od):G[i][j]=H.add(G[i][j],a.dual(x),H.mul(a.dual(y),a.tau(0)))
    HA=frame_family('CA',0);assert H.frame(H.mm(G,HA))==H.ident(5)
    FB=H.frame(H.mm(G,frame_family('CB',1)))
    assert [FB[1][0],FB[1][3],FB[1][4],FB[0][3],FB[0][4]]==actual[9:]
    rows=P['zero_coordinate_rows36'];targets=[18,19,28,29,60,67,72,63,64,68,73,74]
    transform=H.mm(FB,G);itransform=H.mi(transform);ann=H.constmat(P['CB_annihilator17_by36']);R17=[[a.Z]*12 for _ in range(17)]
    for j,target in enumerate(targets):
        inp=[a.Z]*150;inp[target]=a.ONE
        for k,row in enumerate(E['eliminated_zero_rows4']):inp[row]=H.c(E['eliminated_four_coordinate_reconstruction'][k][rows.index(target)])
        pushed=H.physical_push(inp,transform,itransform);v=H.mm(ann,[[pushed[i]] for i in rows])
        for i in range(17):R17[i][j]=v[i][0]
    f17=[H.evalpoly(p,actual) for p in Q['equations14variable_sparse']]
    assert H.mm(R17,[[p] for p in F12])==[[p] for p in f17]
    du=H.add(a.ONE,H.neg(H.mul(kappa[1],kappa[1])));pu=H.add(H.c(5),H.mul(kappa[1],kappa[1]))
    roots=[]
    for name in ['v0','v1','b','c','e']:
        p=R['five_explicit_rational_roots'][name];roots.append(H.mul(H.evalpoly(p['numerator'],kappa),H.power(du,-p['D_denominator_power'])))
    coeff=[H.telescope(p,actual,[0,1,10,11,13],roots) for p in Q['equations14variable_sparse']]
    asinv=H.mi(coeff[:5]);BA=H.mm(coeff[5:],asinv)
    R12=[[H.add(R17[i+5][j],H.neg(H.add(*(H.mul(BA[i][k],R17[k][j]) for k in range(5))))) for j in range(12)] for i in range(12)]
    for i,p in enumerate(R['all12_actual_remaining_and_divided_records']):
        unit=H.scale(H.power(du,p['D_denominator_power']),F(p['constant_denominator'],p['primitive_content_unit']));R12[i]=[H.mul(unit,q) for q in R12[i]]
    hpolys=[p['divided_primitive_polynomial'] for p in R['all12_actual_remaining_and_divided_records']];hv=[[H.evalpoly(p,kappa)] for p in hpolys]
    assert H.mm(R12,[[p] for p in F12])==[[H.mul(kappa[1],p[0])] for p in hv]
    save('phase1-actual-offintersection-R12-B11.json',{'R17':R17,'R12':R12,'all17_offintersection_identities':True,'all12_dH12_R12F12':True,
      'full_ordered_first5_telescopes_used':True},live=[actual,kappa,F12,gs,coeff,BA,hv,first])
    del R17,coeff,BA,asinv,transform,itransform,inp,pushed,v,HA,FB,G
    core=[kappa[j] for j in [0,1,2,6]];inds=[K['source9_variable_order'].index(name) for name in K['selected5_pivot_order']];rr=[]
    for name in K['selected5_pivot_order']:
        p=K['five_actual_rational_roots'][name];rr.append(H.mul(H.evalpoly(p['numerator'],core),H.power(du,-p['D_power']),H.power(pu,-p['P_power'])))
    co=[H.telescope(p,kappa,inds,rr) for p in hpolys];perm=K['selected5_source_H_rows']+[i for i in range(12) if i not in K['selected5_source_H_rows']]
    A5=[co[i] for i in perm[:5]];B7=[co[i] for i in perm[5:]];Ai=H.mi(A5);BA=H.mm(B7,Ai);T=[[a.Z]*12 for _ in range(12)]
    for i in range(5):
        for j in range(5):T[i][perm[j]]=Ai[i][j]
    for i in range(7):
        T[i+5][perm[i+5]]=a.ONE
        for j in range(5):T[i+5][perm[j]]=H.neg(BA[i][j])
    for i,p in enumerate(K['remaining7_actual_core_polynomials']):
        lam=H.scale(H.mul(H.power(du,-p['D_denominator_power']),H.power(pu,-p['P_denominator_power'])),F(p['QQ_unit_content_numerator'],p['QQ_unit_content_denominator']));T[i+5]=[H.mul(H.inv(lam),x) for x in T[i+5]]
    E3=H.ident(12);E3[5][5]=H.inv(du);E3[6][6]=H.inv(du);r=H.add(core[2],H.neg(core[0]),H.neg(core[1]))
    E3[7][5]=H.neg(H.mul(r,H.inv(H.mul(pu,du))));E3[7][7]=H.scale(H.inv(H.mul(pu,du)),F(-2));T=H.mm(E3,T);E4=H.ident(12)
    for j in range(4):
        fac=H.power(du,[2,2,2,3][j])
        for i in range(3):E4[8+j][5+i]=H.neg(H.mul(fac,H.evalpoly(C['cofactor3_by4'][i][j],core)))
    T=H.mm(E4,T);N=H.mm(T,R12)
    assert H.mm(T,hv)==[[p] for p in gs]+[[a.Z] for _ in range(4)]
    assert H.mm(N,[[p] for p in F12])==[[H.mul(kappa[1],p)] for p in gs]+[[a.Z] for _ in range(4)]
    def match1(new,old):return all(new[4*bi+t]==F(old[4*di+t]) for bi,di in LINEAR for t in range(4))
    assert all(match1(v,q) for row,old in zip(R12,D['first_source']['R12']) for v,q in zip(row,old))
    assert all(match1(v,q) for row,old in zip(T,D['first_adapter']['T']) for v,q in zip(row,old))
    assert all(match1(v,q) for row,old in zip(N,D['first_adapter']['N']) for v,q in zip(row,old))
    save('phase2-actual-N-B11.json',{'T_C_B44':T,'N_C_B44':N,'TH12_equals_G8_zero4':True,'NF12_equals_dG8_zero4':True,
      'all_actual_constants_and_firsts_match_frozen_source':True,'new_source_or_gauge_solve':False},live=[R12,hv,gs,co,Ai,BA,E3,E4,first])
    del T,R12,hv,co,Ai,BA,A5,B7,E3,E4,actual,kappa,rr,roots,core,F12,gs,first
    scalar={}
    for o,i,m in requested:
        inp,b=divmod(i,4);outrow,t=divmod(o,4);unit=[F(0)]*44;unit[b]=F(1)
        p=a.mul(N[inp][outrow],tuple(unit));bi={(7,7):2,(0,7):7,(1,7):8,(4,7):9,(5,7):10}[tuple(D['adjoint']['monomials2_36'][m])]
        scalar[o,i,m]=p[4*bi+t]
    values=[]
    for row in D['adjoint']['all619']:
        z=sum((F(c)*scalar[tuple(k)] for k,c in row['actual_L2_source_adjoint']),F(0));values.append(z)
    save('actual-L2-requested-values-and619-native-contributions.json',{'requested_actual_scalar_L2_values':[[list(k),str(v)] for k,v in sorted(scalar.items())],
      'all619_native_L2_contribution_values':values,'all_monic_nu2_basis_terms_retained':True,
      'these_are_NOT_complete619_quartic_values':True,'new_physical_source_jet_or_gauge_solve':False},live=[N])
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'actual_N2_reconstructed_from_frozen_B11_only':True,'all17_offintersection_and_R12_T_N_source_identities_verified':True,
      'new_physical_source_coefficients_or_gauge_solutions':0,'requested_L2_entries':len(scalar),
      'nonzero_native_L2_contributions':sum(bool(x) for x in values),'total619_values_or_HF4_computed':False,
      'whole_saved_plus_live_QQ_max':livepeak,'saved_QQ_slots':saved,'input_sha256':hashes,'code_sha256':H.sha(Path(__file__)),
      'note_sha256':H.sha(note),'evidence_sha256':{p.name:H.sha(p) for p in out.iterdir() if p.name!='phase.json'},'seconds':time.monotonic()-st,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
      'independent_source_value_review':'PENDING'}
    save('report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','new_physical_source_coefficients_or_gauge_solutions','requested_L2_entries','nonzero_native_L2_contributions','whole_saved_plus_live_QQ_max','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
