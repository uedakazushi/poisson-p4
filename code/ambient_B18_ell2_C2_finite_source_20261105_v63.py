#!/usr/bin/env python3
"""One actual finite33 source105 evaluation and right-exact ell² C2 quotient."""
import argparse
import collections
import hashlib
import json
import resource
import signal
import time
import traceback
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,product_columns,eye,srec,blockdiag

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(v):return {i:F(c) for i,c in v}
def guard(v):
    assert all(max(abs(c.numerator).bit_length(),c.denominator.bit_length())<=384 for c in v.values()),'384bit coefficient guard'
    return v
class Algebra33:
    def __init__(self,d):
        self.basis=d['basis33'];self.degrees=d['basis_degrees']
        self.table=[[parse(v) for v in M] for M in d['table_columns_for_each_left_basis']]
        self.one={0:F(1)}
    def mul(self,a,b):
        if not a or not b:return {}
        if len(a)==1 and 0 in a:return guard(scale(b,a[0]))
        if len(b)==1 and 0 in b:return guard(scale(a,b[0]))
        out={}
        for i,c in a.items():
            for j,e in b.items():
                if self.degrees[i]+self.degrees[j]<5:out=add(out,self.table[i][j],c*e)
        return guard(out)
    def inv(self,a):
        c=a.get(0,F(0));assert c,'inverse only at actual nonzero constant'
        n=dict(a);n.pop(0,None);n=scale(n,-1/c)
        out=dict(self.one);v=dict(self.one)
        for _ in range(4):v=self.mul(v,n);out=add(out,v)
        out=guard(scale(out,1/c));assert self.mul(a,out)==self.one
        return out
    def order(self,a):return min((self.degrees[i] for i in a),default=5)

def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);hashes={}
    def read(path):
        p=ROOT/path;hashes[path]=sha(p);return json.loads(p.read_text())
    def write(name,obj):
        p=out/name;p.write_text(json.dumps(obj,separators=(',',':'))+'\n')
        assert sum(q.stat().st_size for q in out.iterdir() if q.is_file())<32*1024*1024,'actual aggregate32MiB guard'
        return sha(p)
    progress={'phase':'input source gate','completed_source_stages':0};write('phase.json',progress)
    pref=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_PREFLIGHT_20261105_v63/report.json');assert pref['all_passed']
    assert pref['whole_slot_bound']==104147
    P=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_PREFLIGHT_20261105_v63/whole-C2-source-and-budget-gate.json')
    table=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_PREFLIGHT_20261105_v63/exact-existing-raw33-multiplication-table.json')
    raw=read('reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json')
    rp=read('reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63_REVIEW/report.json');assert rp['all_passed']
    source_path='reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    joint_path='reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json'
    oldpref_path='reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_SOURCE_PREFLIGHT_20261104_v63/actual-Artin20-source-and-cost-preflight.json'
    S=read(source_path);J=read(joint_path);oldP=read(oldpref_path)
    for p,h in pref['input_sha256'].items():assert sha(ROOT/p)==h
    art=Algebra33(table);basis=table['basis33'];assert basis==raw['actual33_basis'] and table['m5_zero']
    actions=[[parse(c) for c in A] for A in raw['coordinate_actions9_sparse_columns']]
    one=art.one;zvars=[A[0] for A in actions[:7]];r=actions[7][0];ell=actions[8][0];s=add(r,ell,F(1,5))
    targets=zvars+[r,s]
    nodes=S['DAG_scalar_nodes'];ids=P['source_unknown105'];eqids=P['source_equation105'];order=P['source_closure_order'];roots=P['source_requested164_roots']
    assert ids==J['all105_unknown_original_node_order'] and eqids==J['all105_equation_original_node_order']
    uses0=collections.Counter(roots)
    for i in order:
        if nodes[i][0] in ['a','m','n','i']:uses0.update(nodes[i][1:])
    Jinv=[[(i,F(c)) for i,c in enumerate(row) if F(c)] for row in J['fixed105_Jacobian_inverse']]
    assert sum(map(len,Jinv))==866
    max_live=0;max_live_QQ=0
    def evaluate(unknown):
        nonlocal max_live,max_live_QQ
        values={};uses=uses0.copy();inverses=[];qq=0
        for i in order:
            n=nodes[i];op=n[0]
            if op=='v':v=unknown[i]
            elif op=='c':
                c=F(n[1],n[2]) if len(n)==3 else F(n[1]);v={0:c} if c else {}
            elif op=='a':v=add(values[n[1]],values[n[2]])
            elif op=='m':v=art.mul(values[n[1]],values[n[2]])
            elif op=='n':v=scale(values[n[1]],-1)
            elif op=='i':
                v=art.inv(values[n[1]]);inverses.append([i,n[1],srec(values[n[1]]),srec(v)])
            else:raise ValueError(n)
            guard(v);values[i]=v;qq+=len(v);max_live=max(max_live,len(values));max_live_QQ=max(max_live_QQ,qq)
            assert len(values)<=414 and qq<=13662
            if op in ['a','m','n','i']:
                for j in n[1:]:
                    uses[j]-=1
                    if not uses[j]:qq-=len(values[j]);del values[j]
        return {i:values[i] for i in roots},inverses
    def residual(values):
        res=[dict(values[i]) for i in eqids]
        for k,t in enumerate(targets):res[48+k]=guard(add(res[48+k],t,-1))
        return res
    unknown={i:{} for i in ids};stage_hashes={}
    for stage in range(1,5):
        values,inverses=evaluate(unknown);res=residual(values);assert min(map(art.order,res))>=stage
        correction=[]
        for row in Jinv:
            v={}
            for i,c in row:v=add(v,res[i],-c)
            correction.append(guard(v))
        unknown={i:guard(add(unknown[i],v)) for i,v in zip(ids,correction)}
        name=f'stage{stage}-source-residual-and-correction.json'
        stage_hashes[name]=write(name,{'stage':stage,'residual105':[srec(v) for v in res],
          'correction105':[srec(v) for v in correction],'minimum_residual_order':min(map(art.order,res))})
        progress={'phase':'source corrections','completed_source_stages':stage};write('phase.json',progress)
    values,inverses=evaluate(unknown);res=residual(values);assert all(not v for v in res)
    assert all(not u.get(0,0) for u in unknown.values())
    sourcevalues={'actual_CA_F9_targets':[srec(v) for v in targets],'joint105_unknown_order':ids,'joint105_equation_order':eqids,
      'all105_unknown_final_values':[srec(unknown[i]) for i in ids],
      'all_requested164_root_values':[[i,srec(values[i])] for i in roots],
      'all105_final_residuals':[srec(v) for v in res],'all105_equations_exactly_zero':True,
      'all6_unit_inverse_products_exact':inverses,'physical_R0_rows9':S['physical_zero_base_rows9'],
      'physical_R0_readback':[srec(unknown[i]) for i in range(2,11)],'stage_sha256':stage_hashes,
      'finite_actual_function_uniqueness':'same original unit J105, nilpotent algebra m5=0; four fixed-J corrections plus exact equations',
      'not_a_cached_value_or_full_multivariable_jet':True}
    write('actual-source105-ell2-raw33-values.json',sourcevalues)
    # Mandatory exact projection to the independently accepted old ell=0 source.
    H=read('reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json')
    oldsource=read('reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63/actual-source105-Artin20-values.json')
    oldactions=read('reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63/actual-physical12-actions-and-C-special.json')
    oldA=[[parse(v) for v in A] for A in H['multiplication8_column_sparse']]
    Q20=[]
    for t,a,l in basis:
        v={0:F(1)}
        if l:v={}
        for _ in range(a):v=mv(oldA[7],v)
        if t:v=mv(oldA[t-1],v)
        Q20.append(v)
    assert all(product_columns(Q20,A)==product_columns(B,Q20) for A,B in zip(actions,oldA+[[{} for _ in range(20)]]))
    assert all(mv(Q20,unknown[i])==parse(v) for i,v in zip(ids,oldsource['all105_unknown_final_values']))
    oldroot=dict(oldsource['all_requested_original_source_root_values'])
    assert all(mv(Q20,values[i])==parse(oldroot[i]) for i in roots)
    # The canonical quotient maps and their common monomial lifts.
    branches=raw['canonical_branch_quotients6'];Q=[list(map(parse,b['quotient_map6_by33_columns'])) for b in branches]
    branches_actions=[[[parse(c) for c in A] for A in b['coordinate_actions9_sparse_columns']] for b in branches]
    b6=branches[0]['basis_plane_monomials'];assert b6==branches[1]['basis_plane_monomials']
    basisidx={tuple(b):i for i,b in enumerate(basis)};lift=[{basisidx[(0,a,l)]:F(1)} for a,l in b6]
    assert all(mv(q,v)=={i:F(1)} for q in Q for i,v in enumerate(lift))
    def project(k,v):return mv(Q[k],v)
    def block(a,b,c):return dict(list(a.items())+[(i+6,x) for i,x in b.items()]+[(i+12,x) for i,x in c.items()])
    conductor=H['source_conductor_generators_rows_and_ids'];IA=[values[i] for row,i in conductor['CA']];IB=[values[i] for row,i in conductor['CB']]
    assert all(not project(0,v) for v in IA) and all(not project(1,v) for v in IB)
    assert all(not art.mul(a,b) for a in IA for b in IB)
    nuA,nuB=[values[i] for i in oldP['nu_roots_CA_CB']]
    physical=[];physical_data=[]
    for row in oldP['physical12_coefficient_roots']:
        aA,aB,bA,bB=[values[row[k]] for k in ['aA_node','aB_node','bA_node','bB_node']];delta=guard(add(aA,aB,-1))
        assert not project(0,bA) and not project(1,bB)
        assert all(not art.mul(v,bB) for v in IA) and all(not art.mul(v,bA) for v in IB)
        M=[]
        for v in lift:M.append(block(project(0,art.mul(delta,v)),{},scale(art.mul(bB,v),-1)))
        for v in lift:M.append(block({},project(1,art.mul(delta,v)),art.mul(bA,v)))
        pA=scale(art.mul(nuB,bB),-1);pB=art.mul(nuA,bA)
        for j in range(33):
            v={j:F(1)};M.append(block(project(0,art.mul(pA,v)),project(1,art.mul(pB,v)),art.mul(delta,v)))
        assert len(M)==45;physical.append(M)
        physical_data.append({'physical_row':row['physical_row'],'roots':row,'delta':srec(delta),'bA':srec(bA),'bB':srec(bB)})
    rawactions=[[parse(c) for c in A] for A in raw['raw45_coordinate_actions9_sparse_columns']]
    assert all(product_columns(A,B)==product_columns(B,A) for A in physical for B in rawactions)
    assert all(product_columns(A,B)==product_columns(B,A) for i,A in enumerate(physical) for B in physical[i+1:])
    # Raw45 -> raw26 is the direct sum of the canonical ell0 quotient maps.
    B3=[{i:F(1)} if l==0 else {} for i,(a,l) in enumerate(b6)]
    assert b6==[[0,0],[1,0],[2,0],[0,1],[1,1],[2,1]]
    Q26=[]
    Q26.extend([{a:F(1)} if not l else {} for a,l in b6])
    Q26.extend([{3+a:F(1)} if not l else {} for a,l in b6])
    Q26.extend([{6+i:c for i,c in v.items()} for v in Q20])
    oldF=[[parse(v) for v in A] for A in oldactions['physical_actions12_column_sparse']]
    assert all(product_columns(Q26,A)==product_columns(B,Q26) for A,B in zip(physical,oldF))
    columns=[c for A in physical for c in A];piv={};receipts=[]
    for j,orig in enumerate(columns):
        v=dict(orig);steps=[]
        while v:
            p=min(v)
            if p not in piv:
                d=v[p];piv[p]=guard(scale(v,1/d));receipts.append({'source_column':j,'pivot':p,'divide_by':str(d),'subtract':[[k,str(c)] for k,c in steps],'normalized':srec(piv[p])});break
            d=v[p];steps.append((p,d));v=guard(add(v,piv[p],-d))
    free=[i for i in range(45) if i not in piv];pos={i:j for j,i in enumerate(free)};red={}
    for i in reversed(range(45)):
        if i in pos:red[i]={pos[i]:F(1)}
        else:
            v={}
            for j,c in piv[i].items():
                if j!=i:v=add(v,red[j],-c)
            red[i]=guard(v)
    projection=[red[i] for i in range(45)];assert all(not mv(projection,c) for c in columns)
    Cactions=[[mv(projection,A[j]) for j in free] for A in rawactions]
    assert all(product_columns(A,B)==product_columns(B,A) for A in Cactions for B in Cactions)
    assert 12<=len(free)<=16,'mandatory rank4/torsionmu4 C2 dimension bound'
    # D/ell² Jordan data are exact from the supplied ell action.
    E=Cactions[8];assert all(not v for v in product_columns(E,E))
    epiv={}
    for orig in E:
        v=dict(orig)
        while v:
            p=min(v)
            if p not in epiv:epiv[p]=scale(v,1/v[p]);break
            v=add(v,epiv[p],-v[p])
    assert len(free)-len(epiv)==8 and len(epiv)==len(free)-8
    evidence={'physical12_coefficients':physical_data,'physical_actions12_sparse_columns':[[srec(v) for v in A] for A in physical],
      'raw33_to_H20_ell0_columns':[srec(v) for v in Q20],'raw45_to_raw26_ell0_columns':[srec(v) for v in Q26],
      'all105_and164_source_ell0_crosschecks':True,'all12_physical_ell0_intertwiners':True,
      'all64_source_product_checks':True,'all12_branch_lift_and_annihilator_checks':True,
      'all108_coordinate_intertwiners':True,'all66_physical_commutators':True,
      'all540_relation_columns':[srec(v) for v in columns],'rank':len(piv),'fixed_rank_receipts':receipts,
      'actual_C2_dimension':len(free),'C2_basis_raw45_indices':free,
      'all45_projection_columns':[srec(v) for v in projection],
      'C2_coordinate_actions9_sparse_columns':[[srec(v) for v in A] for A in Cactions],
      'ell_action_rank':len(epiv),'ell2_cyclic_summands':len(epiv),'ell_cyclic_summands':len(free)-2*len(epiv),
      'actual_C_DVR_torsion_exponents_at_least2':len(free)-12,
      'no_actual_T2_or_H0_length_claim':True}
    write('actual-physical12-and-C2-right-exact-quotient.json',evidence)
    progress={'phase':'complete','completed_source_stages':4,'actual_C2_dimension':len(free)};write('phase.json',progress)
    # Actual numerical owned arrays are far below the conservative gate.
    qslots=sum(len(v) for v in unknown.values())+sum(len(v) for v in values.values())+sum(len(v) for A in physical for v in A)+sum(len(v) for v in piv.values())
    assert qslots+18584+13662+27720+3465+6930+602+866+185<104147
    report={'id':STEM,'all_passed':True,'actual_C2_dimension':len(free),'actual_C2_physical_relation_rank':len(piv),
      'actual_C_DVR_torsion_cyclic_summands':4,'actual_C_DVR_torsion_exponents_at_least2':len(free)-12,
      'new_job_kind':'finite actual33 unit105 function evaluation, four fixed-J corrections; not a formal jet or gauge solve',
      'all105_equations_exact':True,'all_source_ell0_crosschecks':True,'all_intrinsic_and_physical_actions_exact':True,
      'honest_T2_H0_length_or_K_exponents_computed':False,'source_live_nodes':max_live,'source_live_QQ':max_live_QQ,
      'conservative_whole_slot_bound':104147,'actual_named_owned_QQ':qslots,
      'input_sha256':hashes,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='phase.json'},
      'seconds':time.monotonic()-start,'peak_RSS_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    write('report.json',report)
    print(json.dumps({k:report[k] for k in ['all_passed','actual_C2_dimension','actual_C_DVR_torsion_exponents_at_least2','source_live_QQ','actual_named_owned_QQ','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args()
    try:main(a.output)
    except Exception as e:
        if a.output.exists():(a.output/'failure.json').write_text(json.dumps({'error':str(e),'traceback':traceback.format_exc(),'no_automatic_retry':True})+'\n')
        raise
