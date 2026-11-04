#!/usr/bin/env python3
"""Independent source, covariance and typed-DAG audit; no H1 polynomial job."""
import argparse
from fractions import Fraction as Q
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from sympy.polys.rings import ring
from s61_nonsquare_free281_circuit_review_20261104_v63 import ProofExpressions,read,sha
from s61_nonsquare_joint_Lell_free183_review_20261104_v63 import LinearProof,matrix

ROOT=Path(__file__).resolve().parents[1]
BASE='S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63'
ID=BASE+'_REVIEW'

def main(out):
    began=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    base=ROOT/'reports'/BASE;r=read(base/'report.json');c=read(base/'certificate.json')
    assert r['frozen'] and sha(base/'certificate.json')==r['certificate_sha256']
    for name,h in r['input_sha256'].items():assert sha(ROOT/name)==h,name
    fb=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63';fd=read(fb/'certificate.json')
    jb=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_FREE183_20261104_v63';jd=read(jb/'certificate.json')
    gatepaths=[ROOT/'reports'/n/'report.json' for n in [
        'S61_NONSQUARE_JOINT_LELL_FREE164_20261104_v63_REVIEW',
        'S61_NONSQUARE_LAZY_QUARTIC_ORACLE_20261103_v63_REVIEW',
        'S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63_REVIEW']]
    assert all(read(p)['all_passed'] for p in gatepaths)
    sg=c['scalar_circuit'];assert sg['nodes'][:len(fd['scalar_circuit']['nodes'])]==fd['scalar_circuit']['nodes']
    assert sg['atoms_num_den'][:len(fd['scalar_circuit']['atoms_num_den'])]==fd['scalar_circuit']['atoms_num_den']
    assert c['linear_nodes']==fd['linear_nodes']
    proof=ProofExpressions();sm=[]
    for i,n in enumerate(sg['nodes']):
        if n[0]=='z':assert i==0;v=0
        elif n[0]=='c':v=proof.atom(*sg['atoms_num_den'][n[1]])
        else:
            assert n[0] in ['a','s','m'] and all(0<=j<i for j in n[1:])
            v={'a':proof.plus,'s':proof.minus,'m':proof.times}[n[0]](sm[n[1]],sm[n[2]])
        sm.append(v)
    lin=LinearProof(proof);lm=[];supports=[]
    for i,n in enumerate(c['linear_nodes']):
        if n[0]=='z':assert i==0;v=0;sp=set()
        elif n[0]=='v':v=lin.variable(n[1]);sp={n[1]}
        else:
            assert n[0]=='l' and all(0<=co<len(sm) and 0<=j<i for co,j in n[1])
            v=lin.combine((sm[co],lm[j]) for co,j in n[1]);sp=set().union(*(supports[j] for co,j in n[1] if co))
        lm.append(v);supports.append(sp)
    R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field();unit=[[[0,0,0],'1']]
    enc=lambda p:[[list(m),str(v)] for m,v in p.items()]
    def atom(p):p=K(p);return proof.atom(enc(p.numer),enc(p.denom))
    def parse(text):return K.from_expr(s.sympify(text))
    def atomvalue(data):
        n,de=data;load=lambda t:R.from_dict({tuple(m):s.QQ(v) for m,v in t})
        return K(load(n))/K(load(de))
    C=100*a+18*e**2-9*e-2;F=500*a+90*e**2+3*e-24;B=(50*a+18*e**2-9*e-2)*d+60*e-40;J=5*a*d-6*e-4
    # Every scalar inverse remains in the advertised FULLU localization.
    seen=set();den_count=0
    for _,ts in sg['atoms_num_den']:
        key=json.dumps(ts,separators=(',',':'))
        if key in seen:continue
        seen.add(key);p=R.from_dict({tuple(m):s.QQ(v) for m,v in ts});assert p
        for known in [a,d,F,B,C,J]:
            while not p.rem(known):p=p.exquo(known)
        assert p.is_ground and p;den_count+=1
    mons=jd['new_target239_monomials'];mi={tuple(m):i for i,m in enumerate(mons)}
    kept=fd['remaining122_original_joint_target_columns'];assert kept==c['parent_target122_monomial_indices']
    classes={h:[i for i,j in enumerate(kept) if mons[j][7]==h] for h in range(3)}
    assert [len(classes[i]) for i in range(3)]==[98,23,1]
    assert {str(k):v for k,v in classes.items()}==c['target122_Hdegree_classes']
    NF=fd['full_target_projection239_to122'];assert NF==c['parent_target239_to122_operator']
    vin=NF['input_variables'];pos={v:i for i,v in enumerate(vin)}
    for i,root in enumerate(NF['output122_linear_roots']):assert all(mons[pos[v]][7]==mons[kept[i]][7] for v in supports[root])
    h2=classes[2][0];h2root=NF['output122_linear_roots'][h2]
    assert c['linear_nodes'][h2root]==['v',vin[kept[h2]]]
    Hroots=[lm[NF['output122_linear_roots'][i]] for i in classes[1]]
    defs=jd['old_source319_definitions'];C108=fd['Cnew108_original319_labels']
    assert all(defs[i][0]=='cubic_normalized' for i in fd['Bnew22_original319_labels'])
    q22=[defs[i][1] for i in C108 if defs[i][0]=='quartic'];assert q22==c['retained_quartic_original_rows22'] and len(q22)==22
    oldpair=[next(i for i,x in enumerate(defs) if x==['quartic',r]) for r in [77,131]];assert oldpair==[303,318]
    after164=[i for i in range(165) if i!=fd['separator165_coordinates'][0]]
    pair165=[57+C108.index(i) for i in oldpair]
    assert not set(pair165)&set(fd['separator165_coordinates'])
    pair164=[after164.index(i) for i in pair165];assert pair164==c['literal_free164_pair_positions_for77_131']==[150,163]
    assert c['free163_source164_changed_indices']==[i for i in range(164) if i!=pair164[0]]
    assert c['target121_original122_columns']==[i for i in range(122) if i!=h2]
    # Replay ONLY the existing924+33 archived weighted coefficient recurrences.
    hp=ROOT/'reports/S61_NONSQUARE_QUARTIC_HSYM2_PREFLIGHT_20261103_v63/preflight.json';hd=read(hp);arith=hd['arithmetic']
    opath=ROOT/'reports/S61_NONSQUARE_LAZY_QUARTIC_ORACLE_20261103_v63/oracle.json';oracle=read(opath);raw=oracle['arithmetic']
    assert arith['nodes'][:len(raw['nodes'])]==raw['nodes'] and arith['atoms'][:len(raw['atoms'])]==raw['atoms']
    nodes=arith['nodes'];weights=arith['node_weights'];atoms=[atomvalue(v) for v in arith['atoms']]
    atom_ids={str(v):i for i,v in enumerate(atoms)};node_ids={}
    for i,n in enumerate(nodes):node_ids[(n[0],tuple(tuple(t) for t in n[1])) if n[0]=='l' else tuple(n)]=i
    def constant(v):
        v=K(v)
        return 0 if not v else node_ids['c',atom_ids[str(v)]]
    def combine(ts):
        co={}
        for c,j in ts:
            if c and j:co[j]=co.get(j,K.zero)+c
        co={j:c for j,c in co.items() if c}
        if not co:return 0
        if len(co)==1:
            j,v=next(iter(co.items()))
            if v==K.one:return j
        return node_ids['l',tuple((atom_ids[str(v)],j) for j,v in sorted(co.items()))]
    cache={};vweights=arith['variable_weights']
    def coefficient(i,m):
        key=(i,m)
        if key in cache:return cache[key]
        if not i or sum(x*y for x,y in zip(vweights,m))!=weights[i]:v=0
        else:
            n=nodes[i]
            if n[0]=='v':v=constant(1) if m==tuple(int(j==n[1]) for j in range(8)) else 0
            elif n[0]=='c':v=i if not any(m) else 0
            elif n[0]=='l':v=combine((atoms[c],coefficient(j,m)) for c,j in n[1])
            else:
                assert n[0]=='m';ts=[]
                for left in itertools.product(*(range(j+1) for j in m)):
                    if sum(x*y for x,y in zip(vweights,left))!=weights[n[1]]:continue
                    right=tuple(x-y for x,y in zip(m,left));x,y=coefficient(n[1],left),coefficient(n[2],right)
                    if x and y:ts.append((K.one,node_ids['m',min(x,y),max(x,y)]))
                v=combine(ts)
        cache[key]=v;return v
    rows=oracle['genuine_quartic_original_rows'];roots=dict(zip(rows,oracle['genuine_quartic_roots']))
    assert rows==hd['source_rows'] and len(rows)==33
    light=[0,1,2,3,4,6,7];pairs=list(itertools.combinations_with_replacement(light,2))
    assert [list(p) for p in pairs]==hd['light_pair_columns']
    pairmons=[]
    for i,j in pairs:
        m=[0]*8;m[5]=1;m[i]+=1;m[j]+=1;pairmons.append(tuple(m))
    assert [list(m) for m in pairmons]==hd['H_times_light_pair_monomials']
    assert [[coefficient(roots[r],m) for m in pairmons] for r in rows]==hd['raw33_by28_coefficient_nodes']
    mH2=tuple(2 if i==5 else 0 for i in range(8))
    assert [[r,coefficient(roots[r],mH2)] for r in rows]==hd['raw_H2_coefficient_nodes']
    hs_cache={0:0}
    def hs(i):
        if i in hs_cache:return hs_cache[i]
        assert weights[i]==0;n=nodes[i]
        if n[0]=='c':v=proof.atom(*arith['atoms'][n[1]])
        elif n[0]=='m':v=proof.times(hs(n[1]),hs(n[2]))
        else:assert n[0]=='l';v=proof.total((proof.atom(*arith['atoms'][co]),hs(j)) for co,j in n[1])
        hs_cache[i]=v;return v
    rawH2={r:hs(i) for r,i in hd['raw_H2_coefficient_nodes']}
    rawH={r:[hs(i) for i in rr] for r,rr in zip(rows,hd['raw33_by28_coefficient_nodes'])}
    adpath=ROOT/'reports/S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63/adapter.json';ad=read(adpath)
    D=list(map(parse,ad['Tden']));DH=D[5];oldm2=[x[2] for x in defs[:29]];oldmi={tuple(m):i for i,m in enumerate(oldm2)}
    newm2=[m for _,m in jd['new_prime_quadric_source174_definitions'][:29]]
    TB=[{i:sm[j] for i,j in row} for row in jd['monomial_old_to_new_scalar_roots']]
    def project(values):
        supplied={v:0 for v in vin}
        for m,co in zip(newm2,values):
            if co:
                assert m[7]==0;n=m.copy();n[7]=1;supplied[vin[mi[tuple(n)]]]=co
        return lin.substitute(Hroots,supplied)
    qproject={}
    for row in rows:
        values=[0]*29
        for i,(monH,(v,w)) in enumerate(zip(pairmons,pairs)):
            m=list(monH);m[5]=0;oldcol=oldmi[tuple(m)]
            scaled=proof.times(rawH[row][i],atom(D[v]*D[w]/DH))
            for j,co in TB[oldcol].items():values[j]=proof.plus(values[j],proof.times(scaled,co))
        qproject[row]=project(values)
    assert [[r,[sm[i] for i in rr]] for r,rr in c['all33_raw_quartic_H1_post117_scalar_roots']]==[[r,qproject[r]] for r in rows]
    assert [[r,sm[i]] for r,i in c['all33_raw_H2_coefficient_scalar_roots']]==[[r,rawH2[r]] for r in rows]
    hpath=ROOT/'reports/S61_BINARY_HYPERPLANE_HEAVY_SQUARE_PIVOT_20261102.json';h=read(hpath)
    u,v=map(parse,h['row_combination_coefficients']);b0,b1=map(parse,h['coefficient_formulas']);U=K((50*a*B*C)**2)
    hcore_cache={0:K.zero}
    def hcore(i):
        if i in hcore_cache:return hcore_cache[i]
        assert weights[i]==0;n=nodes[i]
        if n[0]=='c':value=atoms[n[1]]
        elif n[0]=='m':value=hcore(n[1])*hcore(n[2])
        else:value=sum((atoms[co]*hcore(j) for co,j in n[1]),K.zero)
        hcore_cache[i]=value;return value
    h2dict=dict(hd['raw_H2_coefficient_nodes'])
    assert hcore(h2dict[77])==b0 and hcore(h2dict[131])==b1
    assert u*b0+v*b1==U and U
    M=[[u,v],[-b1,b0]];Mi=[[b0/U,-v/U],[b1/U,u/U]]
    for L,Ri in [(M,Mi),(Mi,M)]:assert all(sum(L[i][k]*Ri[k][j] for k in range(2))==int(i==j) for i in range(2) for j in range(2))
    corecache={0:K.zero}
    def core(i):
        if i in corecache:return corecache[i]
        n=sg['nodes'][i]
        if n[0]=='c':value=atomvalue(sg['atoms_num_den'][n[1]])
        elif n[0]=='a':value=core(n[1])+core(n[2])
        elif n[0]=='s':value=core(n[1])-core(n[2])
        else:value=core(n[1])*core(n[2])
        corecache[i]=value;return value
    # The stored -b1 is an explicit subtraction node, not an atomic negative.
    assert [[core(i) for i in row] for row in c['H2_source_pair_change_scalar_roots']]==M
    assert [[core(i) for i in row] for row in c['H2_source_pair_inverse_scalar_roots']]==Mi
    assert [sm[i] for i in c['H2_literal_coefficient_roots_from_reviewed_formulas']]==[atom(b0),atom(b1)]
    assert sm[c['H2_known_unit_root']]==atom(U)
    chi=[proof.plus(proof.times(atom(u),x),proof.times(atom(v),y)) for x,y in zip(qproject[77],qproject[131])]
    psi=[proof.plus(proof.times(proof.minus(0,atom(b1)),x),proof.times(atom(b0),y)) for x,y in zip(qproject[77],qproject[131])]
    assert [sm[i] for i in c['chi_H1_post117_scalar_roots']]==chi
    q21=[psi]+[[proof.minus(x,proof.times(proof.times(rawH2[row],atom(1/U)),y)) for x,y in zip(qproject[row],chi)] for row in q22 if row not in [77,131]]
    assert len(q21)==21 and [[sm[i] for i in row] for row in c['all_remaining21_quartic_H1_scalar_roots']]==q21
    # Bind every H-linear coefficient directly to actual frozen57, not just metadata.
    cp=ROOT/'reports/S61_NONSQUARE_QUARTIC_122_CIRCUIT_20261103_v63/circuit.json';cc=read(cp)
    src=read(ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json');cmons=src['normal_monomials']
    heavy=[]
    heavy_by_row={rr[0]['basis_row']:rr for rr in cc['heavy_coefficient_entries']}
    assert sorted(heavy_by_row)==list(range(41,57))
    for basis_row in range(41,57):
        rr=heavy_by_row[basis_row]
        row=[]
        for j,entry in enumerate(rr):
            assert entry['basis_row']==basis_row
            m=cmons[entry['column']];assert m[5]==1 and m[light[j]]==1 and sum(m)==2
            t=entry['sparse_encoding'].split();n=int(t[0]);assert len(t)==1+4*n
            row.append(proof.atom([[list(map(int,t[k+1:k+4])),t[k]] for k in range(1,len(t),4)],unit))
        heavy.append(row)
    assert len(heavy)==16 and [[sm[i] for i in row] for row in c['genuine_cubic16_Hlinear_old_scalar_roots']]==heavy
    raw57=ROOT/'reports/S61_NONSQUARE_ELL_SIX_BORDERS_20261104_v63/frozen57-source.txt'
    basis=read(ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json');assert sha(raw57)==basis['remaining_matrix_sha256']
    hc=[]
    for j in light:m=[0]*8;m[5]=m[j]=1;hc.append(cmons.index(m))
    actual={}
    with raw57.open() as f:
        assert f.readline().split()==['57','91','57']
        for i in range(57):
            for j in range(91):
                line=f.readline();assert line
                if i>=41 and j in hc:
                    t=line.split();n=int(t[0]);actual[i,j]=proof.atom([[list(map(int,t[k+1:k+4])),t[k]] for k in range(1,len(t),4)],unit)
        assert not f.read().strip()
    assert [[actual[41+i,j] for j in hc] for i in range(16)]==heavy
    np=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_PREFLIGHT_20261104_v63/coordinate-inverse-numerators7x7.txt';Nraw=matrix(np);delta=K(3200*a*C*d*d)/81
    N=[[proof.times(proof.atom(p,unit),atom(1/delta)) for p in row] for row in Nraw]
    oldB=[N[0],N[1],N[2],N[6],N[3],[proof.minus(proof.plus(N[5][j],N[1][j]),proof.times(atom(3),N[0][j])) for j in range(7)],N[4]]
    hnew=[[proof.total((rr[j],oldB[j][k]) for j in range(7)) for k in range(7)] for rr in heavy]
    assert [[sm[i] for i in rr] for rr in c['genuine_cubic16_Hlinear_joint_scalar_roots']]==hnew
    newmi={tuple(m):i for i,m in enumerate(newm2)};c112=[]
    for i in range(16):
        for oldj in range(7):
            values=[0]*29
            for k in range(7):
                for l in range(7):
                    m=[0]*8;m[k]+=1;m[l]+=1;j=newmi[tuple(m)]
                    values[j]=proof.plus(values[j],proof.times(hnew[i][k],oldB[oldj][l]))
            c112.append(project(values))
    assert [[sm[i] for i in rr] for rr in c['all112_genuine_cubic_H1_scalar_roots']]==c112
    assert c['all112_cubic_source319_definitions']==[x for x in defs if x[0]=='cubic_normalized']
    qb=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_QPAIR_UNITS_20261104_v63/checkpoint.json';qs=read(qb)
    pb=ROOT/'reports/S61_NONSQUARE_JOINT_LELL_OVERLAP_UNIT_20261104_v63/remaining-pair-block.json';pair=read(pb)
    assert pair['source_rows'][57]==[5,[0,0,0,0,0,0,0,1]] and 57 in qs['remaining_rows']
    qH=[proof.atom(*qs['matrix58x136'][57][qs['remaining_columns'][i]]) for i in classes[1]]
    for row in qs['remaining_rows']:
        if row!=57:assert all(not proof.atom(*qs['matrix58x136'][row][qs['remaining_columns'][i]]) for i in classes[1])
    assert [sm[i] for i in c['remaining_q5H_H1_scalar_roots']]==qH
    block=[qH]+c112+q21
    assert len(block)==134 and all(len(row)==23 for row in block)
    assert [[sm[i] for i in row] for row in c['complete_H1_image_block_scalar_roots']]==block
    # Full chi polynomial is source-bound, while no pure coefficient is requested.
    pids={('zero',):0};pdata=[('zero',)]
    def intern(key):
        if key not in pids:pids[key]=len(pdata);pdata.append(key)
        return pids[key]
    def pv(i):return intern(('var',i))
    def pc0(v):return 0 if not v else intern(('constant',v))
    def pm(x,y):return 0 if not x or not y else intern(('product',x,y))
    def pl(ts):
        ts=tuple((co,i) for co,i in ts if co and i)
        if not ts:return 0
        if len(ts)==1 and ts[0][0]==proof.one:return ts[0][1]
        return intern(('linear',ts))
    mapped=[]
    for i,n in enumerate(c['chi_polynomial_nodes']):
        if n[0]=='z':assert i==0;v0=0
        elif n[0]=='v':v0=pv(n[1])
        elif n[0]=='c':v0=pc0(sm[n[1]])
        elif n[0]=='m':assert all(j<i for j in n[1:]);v0=pm(mapped[n[1]],mapped[n[2]])
        else:assert n[0]=='l';v0=pl((sm[co],mapped[j]) for co,j in n[1])
        mapped.append(v0)
    images=[]
    for i in range(8):
        if i==5:v0=pl([(atom(DH),pv(7))])
        else:v0=pl((proof.times(atom(D[i]),oldB[light.index(i)][j]),pv(j)) for j in range(7))
        images.append(v0)
    p_cache={0:0}
    def psource(i):
        if i in p_cache:return p_cache[i]
        n=raw['nodes'][i]
        if n[0]=='v':v0=images[n[1]]
        elif n[0]=='c':v0=pc0(proof.atom(*raw['atoms'][n[1]]))
        elif n[0]=='m':v0=pm(psource(n[1]),psource(n[2]))
        else:v0=pl((proof.atom(*raw['atoms'][co]),psource(j)) for co,j in n[1])
        p_cache[i]=v0;return v0
    p77=pl([(atom(1/DH**2),psource(roots[77]))]);p131=pl([(atom(1/DH**2),psource(roots[131]))])
    pchi=pl([(atom(u),p77),(atom(v),p131)])
    assert mapped[c['chi_raw77_Qhat_polynomial_root']]==p77 and mapped[c['chi_raw131_Qhat_polynomial_root']]==p131
    assert mapped[c['chi_full_joint_weight4_polynomial_root']]==pchi
    assert c['pure98_coefficient_queries_performed']==0 and r['pure98_coefficient_queries']==0
    out.mkdir(parents=True,exist_ok=False)
    ep=out/'independent-proof.json'
    evidence={'literal_old319_pair':[303,318],'literal_free164_pair':pair164,'target_classes':[98,23,1],
        'source_pair_determinant':'(50*a*B*C)^2','all228_unique_scalar_denominators_known_FULLU_units':den_count==228,
        'all924_archived_Hlightpair_and33_H2_weighted_coefficient_recurrences_verified':True,
        'all112_actual_frozen57_Hcoefficients_bound':True,'all33_quartic_T34_to_Tn_joint_NF117_transports_verified':True,
        'all112_cubic_and21_quartic_plus_q5H_134by23_circuit_roots_match':True,
        'full_chi_77_131_weight4_polynomial_circuit_bound_without_pure98_queries':True,
        'actual_same_cokernel_free163_to121':True,'H1_redundant_generators_not_free_source163_basis':True,
        'no_H1_rank_or_full239_or_pure98_expansion_or_newpoint_claim':True}
    ep.write_text(json.dumps(evidence,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),base/'report.json',base/'certificate.json',fb/'certificate.json',jb/'certificate.json',hp,opath,adpath,hpath,cp,np,qb,pb,raw57,*gatepaths,
        ROOT/'code/s61_nonsquare_free281_circuit_review_20261104_v63.py',ROOT/'code/s61_nonsquare_joint_Lell_free183_review_20261104_v63.py']
    result={'id':ID,'all_passed':True,'literal_H2_knownunit_pair_and_global_samecokernel163_to121_accepted':True,
        'all924_old_H_queries_and_actual112_heavy_coefficients_source_bound':True,
        'complete_H1_projection_134by23_circuit_and_covariance_accepted':True,
        'full_lazy_chi_keeps_actual_pure98_coefficients_unqueried':True,
        'redundant134_rows_do_not_automatically_give_free163_minus_number_of_pivots':True,
        'H1_rank_full121_full239_or_new_kernel_count_claim':False,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{ep.name:sha(ep)},
        'seconds':time.monotonic()-began,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:result[k] for k in ['all_passed','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
