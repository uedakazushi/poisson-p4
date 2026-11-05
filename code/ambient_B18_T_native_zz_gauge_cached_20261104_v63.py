#!/usr/bin/env python3
"""Nine-zz canonical gauge extension from composed105 source2 cache only."""
import argparse,collections,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_comparison_X_origin_20261104_v63 as H
from ambient_B18_native128_first_20261104_v63 import PATHS as OLD_PATHS,Algebra as FirstAlgebra
from ambient_B18_T_source2_cache_composition_20261104_v63 import MON,IDX,mul as sparsemul
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63'
def counts(x):
    if isinstance(x,F):return 1
    if isinstance(x,dict):return sum(counts(v) for v in x.values())
    if isinstance(x,(list,tuple)):return sum(counts(v) for v in x)
    return 0
def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [encode(v) for v in x]
    return x
class Finite:
    def __init__(self,m=0,nus=()):
        self.m=m;self.b=1<<m;self.n=20*self.b;self.Z=(F(0),)*self.n;self.ONE=self.c(1);self.nus=nus;self.table={}
        for i in range(self.n):
            ai,ti=divmod(i,self.b)
            for j in range(self.n):
                aj,tj=divmod(j,self.b);v=sparsemul({ai:F(1)},{aj:F(1)});carry=ti&tj;mask=ti^tj
                for k in range(m):
                    if carry&(1<<k):v=sparsemul(v,{a:c for a,c in enumerate(nus[k]) if c})
                if v:self.table[i,j]=[(self.b*a+mask,c) for a,c in v.items()]
    def c(self,x):return (F(x),)+(F(0),)*(self.n-1)
    def add(self,*ps):
        v=[F(0)]*self.n
        for p in ps:
            for i,c in enumerate(p):
                if c:v[i]+=c
        return tuple(v)
    def neg(self,p):return tuple(-c for c in p)
    def scale(self,p,c):return tuple(c*a for a in p)
    def mul2(self,p,q):
        if p==self.Z or q==self.Z:return self.Z
        if p==self.ONE:return q
        if q==self.ONE:return p
        v=[F(0)]*self.n;qq=[(j,c) for j,c in enumerate(q) if c]
        for i,c in enumerate(p):
            if c:
                for j,d in qq:
                    for k,a in self.table.get((i,j),()):v[k]+=c*d*a
        return tuple(v)
    def mul(self,*ps):
        v=self.ONE
        for p in ps:v=self.mul2(v,p)
        return v
    def power(self,p,n):
        if n<0:return self.power(self.inv(p),-n)
        if n>=9 and not p[0]:return self.Z
        v=self.ONE
        while n:
            if n&1:v=self.mul2(v,p)
            n//=2
            if n:p=self.mul2(p,p)
        return v
    def inv(self,p):
        assert p[0];q=list(p);q[0]=F(0);q=self.scale(tuple(q),-1/p[0]);v=self.ONE;t=self.ONE
        for _ in range(8):t=self.mul(t,q);v=self.add(v,t)
        v=self.scale(v,1/p[0]);assert self.mul(p,v)==self.ONE;return v
    def dual(self,p):
        v=[F(0)]*self.n
        for i,c in enumerate(p):v[self.b*i]=F(c)
        return tuple(v)
    def tau(self,k):
        v=[F(0)]*self.n;v[1<<k]=F(1);return tuple(v)
    def install(self):
        for k in ['c','add','neg','scale','mul','power','inv']:setattr(H,k,getattr(self,k))
        H.Z=self.Z;H.ONE=self.ONE;H.TA=self.tau(0) if self.m else self.Z;H.TB=self.tau(1) if self.m>1 else self.Z
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False);saved=0;peak=0
    paths={k:OLD_PATHS[k] for k in ['source','interface','proper','second5','preparation']}
    paths.update({
      'source20':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63/phase1-actual105-source2-and-copy.json',
      'physical11':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase1-source-B11-values.json',
      'gauge11':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase2-gauge30-B11-correction.json',
      'proper11':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase3-actual-proper-G8-B11.json',
      'native11':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/actual-native128-effective-S2-B11.json',
      'B11':'reports/AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT_20261104_v63/actual-B11-source-and-gauge-preflight.json',
      'J30':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase2-fixed30-Jacobian.json',
      'firstpreflight':'reports/AMBIENT_B18_NATIVE128_FIRST_SOURCE_PREFLIGHT_20261104_v63/actual-native128-first-input-and-mixed-index-profile.json',
      'adapter1':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase5-native-adapter-firsts.json',
      'Soperator':OLD_PATHS['Q_operator'],
      'lower':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63_OLDB_AD/phase4-all619-lower-knownparts-and-symbols.json'})
    def save(name,data,live=()):
        nonlocal saved,peak
        q=counts(data);peak=max(peak,saved+q+counts(live));assert peak<150000,('QQ_GUARD',name,peak)
        p=out/name;p.write_text(json.dumps(encode(data),separators=(',',':'))+'\n');saved+=q;assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
        (out/'phase.json').write_text(json.dumps({'phase':name,'saved_QQ':saved,'peak_saved_live_QQ':peak,'seconds':time.monotonic()-start})+'\n')
    hashes={p:H.sha(ROOT/p) for p in paths.values()};save('source-manifest.json',{'input_sha256':hashes,'no_new_physical_implicit_solve':True})
    D={k:H.read(ROOT/p) for k,p in paths.items()};S,I,P,K,prep=[D[k] for k in ['source','interface','proper','second5','preparation']]
    C=D['source20'];assert C['base20_monomials_original8']==[list(m) for m in MON] and C['all105_equations_verified']
    def dense(p):
        v=[F(0)]*20
        for i,c in p:v[i]=F(c)
        return tuple(v)
    first={i:dense(v) for i,v in C['all105_actual_source2_values']};bm=[tuple(i for i,n in enumerate(m) for _ in range(n)) for m in D['B11']['B11_basis_monomials8']]
    oldidx=[IDX[m] for m in bm];oldset=set(oldidx);newquad=[i for i,m in enumerate(MON) if len(m)==2 and i not in oldset];assert len(newquad)==9
    def embed11(p):
        v=[F(0)]*20
        for i,c in enumerate(p):v[oldidx[i]]=F(c)
        return tuple(v)
    nodes=S['DAG_scalar_nodes'];roots=[i for f in S['families'].values() for p in f['complete_normalized_physical_w_coefficients'][0] for i in p]+[f['actual_tracefree_nu'] for f in S['families'].values()]
    needed=set();todo=list(roots)
    while todo:
        i=todo.pop()
        if i in needed:continue
        needed.add(i)
        if i in first:continue
        n=nodes[i]
        if n[0] in ['a','m','n','i']:todo.extend(n[1:])
        else:assert n[0]=='c',('UNCACHED_PHYSICAL_SOURCE',i,n)
    uses=collections.Counter()
    for i in needed:
        if i not in first and nodes[i][0] in ['a','m','n','i']:uses.update(nodes[i][1:])
    uses.update(roots);base=Finite();vals={};liveness=0;sourcepeak=0
    for i in sorted(needed):
        n=nodes[i];o=n[0]
        if i in first:v=first[i]
        elif o=='c':v=base.c(F(n[1],n[2]))
        elif o=='a':v=base.add(vals[n[1]],vals[n[2]])
        elif o=='m':v=base.mul(vals[n[1]],vals[n[2]])
        elif o=='n':v=base.neg(vals[n[1]])
        elif o=='i':v=base.inv(vals[n[1]])
        else:raise AssertionError(n)
        vals[i]=v;liveness+=20;sourcepeak=max(sourcepeak,liveness);assert sourcepeak+15000<150000
        if i not in first:
            for j in n[1:] if o in ['a','m','n','i'] else []:
                uses[j]-=1
                if not uses[j]:del vals[j];liveness-=20
    pairs={f:[[vals[x],vals[y]] for x,y in data['complete_normalized_physical_w_coefficients'][0]] for f,data in S['families'].items()}
    nu=[vals[f['actual_tracefree_nu']] for f in S['families'].values()];assert list(S['families'])==['CA','CB']
    assert all([v[i] for i in oldidx]==list(map(F,q)) for f in pairs for p,qpair in zip(pairs[f],D['physical11']['physical150_tracepairs_B11'][f]) for v,q in zip(p,qpair))
    assert all([v[i] for i in oldidx]==list(map(F,q)) for v,q in zip(nu,D['physical11']['nu_A_B_B11']))
    save('phase1-actual-physical-source-base20.json',{'base20_monomials':MON,'physical150_tracepairs_base20':pairs,'nu_A_B_base20':nu,'all_B11_projections_match':True},live=first)
    del vals
    ev=[embed11(v) for v in D['gauge11']['even15_B11']];od=[embed11(v) for v in D['gauge11']['odd15_B11']]
    GL=D['firstpreflight']['GL18_source_parameter_positions'];positions=I['actual_group15_positions'];rows=P['zero_coordinate_rows36'];ar=P['CA_constant_left_inverse_rows19'];fixed=P['CA_fiber10_fixed_coefficient_indices'];free=P['CA_free_fiber9_coefficient_indices'];centers=list(map(F,P['actual_CA_source19_center']))
    def frame_family(a,family,tauindex):
        a.install();f=S['families'][family];p=H.add(a.tau(tauindex),H.scale(a.dual(first[f['p_quadratic_relation']['s']]),F(1,2)));G=H.ident(5)
        for name,(aa,bb) in f['coefficient_unknown_pairs'].items():
            if name.startswith('gamma') and name!='gamma0':
                i,j=divmod(GL[int(name[5:])-1],5);G[i][j]=H.add(G[i][j],a.dual(first[aa]),H.mul(a.dual(first[bb]),p))
        core=H.ident(5)
        if family=='CA':core[0][3]=H.scale(p,F(-2));core[1][0]=p;core[1][3]=H.neg(H.mul(p,p))
        else:core[0][3]=H.scale(p,F(1,2));core[1][0]=p;core[1][3]=H.scale(H.mul(p,p),F(1,4))
        return H.mm(G,core)
    def physical(a,family,tauindex):
        a.install();return [H.add(a.dual(x),H.mul(a.dual(y),a.tau(tauindex))) for x,y in pairs[family]]
    def gauge(a):
        a.install();G=H.ident(5)
        for (i,j),x,y in zip(positions,ev,od):G[i][j]=H.add(G[i][j],a.dual(x),H.mul(a.dual(y),a.tau(0)))
        return G
    def constraints(a,G,HA,PA):
        a.install();FA=H.frame(H.mm(G,HA));tr=H.mm(FA,G);moved=H.physical_push(PA,tr,H.mi(tr));own=H.mm(H.constmat(P['CA_left_inverse19']),[[moved[rows[j]]] for j in ar]);own=[q[0] for q in own]
        return [FA[1][0],FA[1][3],FA[1][4],FA[0][3],FA[0][4]]+[H.add(own[j],H.c(-centers[j])) for j in fixed]
    a40=Finite(1,[nu[0]]);HA=frame_family(a40,'CA',0);PA=physical(a40,'CA',0);res=constraints(a40,gauge(a40),HA,PA)
    assert all(not q[2*b+t] for q in res for b in range(20) if b not in newquad for t in range(2))
    rhs=[[res[i][2*b+t] for b in newquad] for t in range(2) for i in range(15)];Ji=[[F(c) for c in row] for row in D['J30']['J30_inverse']]
    corr=[[-sum((Ji[i][j]*rhs[j][k] for j in range(30)),F(0)) for k in range(9)] for i in range(30)]
    for i in range(15):
        x=list(ev[i]);y=list(od[i])
        for k,b in enumerate(newquad):x[b]+=corr[i][k];y[b]+=corr[15+i][k]
        ev[i]=tuple(x);od[i]=tuple(y)
    assert all(q==a40.Z for q in constraints(a40,gauge(a40),HA,PA))
    save('phase2-canonical270-zz-gauge.json',{'new_quadratic_indices9':newquad,'new_quadratic_pairs9':[MON[i] for i in newquad],
      'RHS30_by9':rhs,'new_correction30_by9':corr,'actual_even15_base20':ev,'actual_odd15_base20':od,
      'all30_scalar_gauge_constraints_exact':True,'all_B11_values_retained':True,'fixed_J30_inverse_reused':True},live=[pairs,first])
    del HA,PA,res,rhs,Ji,corr
    a=Finite(2,nu);G=gauge(a);HA=frame_family(a,'CA',0);HB=frame_family(a,'CB',1);assert H.frame(H.mm(G,HA))==H.ident(5);FB=H.frame(H.mm(G,HB))
    PA=physical(a,'CA',0);PB=physical(a,'CB',1);moved=H.physical_push(PA,G,H.mi(G));own=H.mm(H.constmat(P['CA_left_inverse19']),[[moved[rows[j]]] for j in ar]);own=[q[0] for q in own];assert all(own[j]==H.c(centers[j]) for j in fixed)
    actual=[H.add(own[j],H.c(-centers[j])) for j in free]+[FB[1][0],FB[1][3],FB[1][4],FB[0][3],FB[0][4]];kappa=[actual[j] for j in I['retained_source9_indices']]
    F12=[H.add(PA[j],H.neg(PB[j])) for j in [18,19,28,29,60,67,72,63,64,68,73,74]]
    du=H.add(a.ONE,H.neg(H.mul(kappa[1],kappa[1])));pu=H.add(H.c(5),H.mul(kappa[1],kappa[1]));core=[kappa[j] for j in [0,1,2,6]];inds=[K['source9_variable_order'].index(n) for n in K['selected5_pivot_order']];rr=[]
    for n in K['selected5_pivot_order']:
        p=K['five_actual_rational_roots'][n];rr.append(H.mul(H.evalpoly(p['numerator'],core),H.power(du,-p['D_power']),H.power(pu,-p['P_power'])))
    gs=[H.add(kappa[j],H.neg(v)) for j,v in zip(inds,rr)]+[H.evalpoly(p,core) for p in prep['generators3']]
    for name,new,old in [('actual14_C_B44',actual,D['proper11']['actual14_C_B44']),('kappa9_C_B44',kappa,D['proper11']['kappa9_C_B44']),('physicalF12_C_B44',F12,D['proper11']['physicalF12_C_B44']),('actualG8_C_B44',gs,D['proper11']['actualG8_C_B44'])]:
        assert all([p[4*b+t] for b in oldidx for t in range(4)]==list(map(F,q)) for p,q in zip(new,old)),name
    save('phase3-actual-proper-G8-base20.json',{'actual14_C80':actual,'kappa9_C80':kappa,'physicalF12_C80':F12,'actualG8_C80':gs,
      'all_B11_and_first_projections_match':True,'all_canonical_CA_source_and_frame_conditions_exact':True},live=[pairs,first])
    del PA,PB,HA,HB,FB,G,moved,own,actual,kappa,F12,core,rr,du,pu,ev,od
    # Scalar L0+L1 stays exactly the frozen actual first matrix.
    pf=D['firstpreflight'];ad=FirstAlgebra(8,2,[[F(c) for c in p] for p in pf['nu_A_B_first9']]);oldN=D['adapter1']['N'];L={}
    for j in range(48):
        jn,b=divmod(j,4);unit=[F(0)]*36;unit[b]=F(1)
        for i in range(12):
            v=ad.mul(tuple(map(F,oldN[jn][i])),tuple(unit))
            for t in range(4):
                p=[F(0)]*20;p[0]=v[t]
                for k in [0,1,4,5,7]:p[IDX[(k,)]]=v[4*(k+1)+t]
                if any(p):L[4*i+t,j]=tuple(p)
    a.install();Sn=[[a.Z]*32 for _ in range(12)]
    for i,j,g,sign in D['Soperator']['native_S12_by32']['koszul28_in_first8_rows']:Sn[i][j]=H.scale(gs[g],F(sign))
    for i,j,c in D['Soperator']['native_S12_by32']['free4_zero_source_columns']:Sn[i][j]=H.c(c)
    lower=D['lower'];pairs36=list(__import__('itertools').combinations_with_replacement(range(8),2));req={tuple(k) for row in lower['all619_uncomputed_native_second_source_functionals'] for k,c in row}
    oldnv=[{row:{int(i):F(c) for i,c in p} for row,p in col} for col in D['native11']['source_columns128_rows48_sparse_B11_effective_second']]
    oldpair={MON[oldidx[i]]:i for i in range(11)};lookup={};newreq=collections.defaultdict(set)
    for kind,j,row,m in req:
        pair=pairs36[m]
        if pair in oldpair:lookup[kind,j,row,m]=oldnv[j].get(row,{}).get(oldpair[pair],F(0))
        else:newreq[j].add(row)
    for j,rowsneeded in newreq.items():
        jn,b=divmod(j,4);unit=[F(0)]*80;unit[b]=F(1);sv={}
        for i in range(12):
            v=a.mul(Sn[i][jn],tuple(unit))
            for t in range(4):
                p=tuple(v[4*k+t] for k in range(20))
                if any(p):sv[4*i+t]=p
        for row in rowsneeded:
            p=base.Z
            for k,q in sv.items():
                if (row,k) in L:p=base.add(p,base.mul(L[row,k],q))
            for key in req:
                if key[1]==j and key[2]==row and pairs36[key[3]] not in oldpair:lookup[key]=p[IDX[pairs36[key[3]]]]
    assert set(lookup)==req
    values=[]
    assert lower['all_ordinary_honest_second_sensitivities_zero']
    for known,row in zip(lower['all619_oldB_known_parts'],lower['all619_uncomputed_native_second_source_functionals']):values.append(F(known)+sum((F(c)*lookup[tuple(k)] for k,c in row),F(0)))
    save('all619-lower-contribution-values.json',{'new_native_second_keys':[[list(k),str(v)] for k,v in sorted(lookup.items()) if pairs36[k[3]] not in oldpair],
      'seeded_B11_key_count':sum(pairs36[k[3]] in oldpair for k in req),'all619_actual_lower_contribution_values':values,
      'effective_native_values_used_only_after_zero_origin_source_combinations':True,'individual_actual_nativeX2_not_claimed':True,
      'not_complete619_quartic_values_or_HF4':True},live=[pairs,first,gs,L,Sn])
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'new_gauge_second_QQ_coefficients':270,'new_physical_implicit_coefficients':0,
      'all30_canonical_gauge_constraints_and_B11_projections_exact':True,'all619_lower_contributions_now_fixed':True,
      'new_native_second_queries':sum(pairs36[k[3]] not in oldpair for k in req),'nonzero_lower_contributions':sum(bool(v) for v in values),
      'actual619_complete_quartic_values_or_HF4_computed':False,'source3_or6_or_N2zz_computed':False,
      'saved_QQ_slots':saved,'peak_saved_plus_live_QQ':peak,'input_sha256':hashes,'code_sha256':H.sha(Path(__file__)),
      'note_sha256':H.sha(note),'evidence_sha256':{p.name:H.sha(p) for p in out.iterdir() if p.name!='phase.json'},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'independent_review':'PENDING'}
    save('report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','new_gauge_second_QQ_coefficients','all619_lower_contributions_now_fixed','new_native_second_queries','nonzero_lower_contributions','peak_saved_plus_live_QQ','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
