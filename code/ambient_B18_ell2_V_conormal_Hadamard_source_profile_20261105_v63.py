#!/usr/bin/env python3
"""Exact selected conormal source adjoints; all V Taylor callbacks remain formal."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL2_V_CONORMAL_HADAMARD_SOURCE_PROFILE_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(v):return {i:F(c) for i,c in v}
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);hashes={}
    def read(p):hashes[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
    r=read('reports/AMBIENT_B18_ELL2_C2_FIRST_BETTI_SOURCE_PROFILE_20261105_v63/report.json');assert r['all_passed'] and r['actual_beta1_C2_over_B']==26
    d2=read('reports/AMBIENT_B18_ELL2_C2_FIRST_BETTI_SOURCE_PROFILE_20261105_v63/actual-d2-origin-rank-and-source-factorization.json')
    res=read('reports/AMBIENT_B18_ELL2_C2_RESOLUTION_SOURCE_PREFLIGHT_20261105_v63/actual-C2-resolution-and11-unit-source-cancellation.json')
    c=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json')
    chain=read('reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json')
    # A fixed independent100-column basis, not a new rank search.
    anchors=res['source_anchors11'];selected=d2['fixed_original_minor_columns'];allcols=list(map(parse,d2['all_d2_constant_columns']))
    fixed=[{i:F(-1)} for i in anchors]+[allcols[i] for i in selected]
    piv={};receipts=[]
    for j,orig in enumerate(fixed):
        v=dict(orig);steps=[]
        while v and min(v) in piv:
            p=min(v);a=v[p];steps.append([p,str(a)]);v=add(v,piv[p],-a)
        assert v,'given fixed100 basis must stay independent'
        p=min(v);a=v[p];piv[p]=scale(v,1/a)
        receipts.append({'fixed_basis_column':j,'pivot':p,'divide_by':str(a),'subtract':steps,'normalized':srec(piv[p])})
    assert len(piv)==100
    free=[i for i in range(126) if i not in piv];assert len(free)==26;fi={p:i for i,p in enumerate(free)}
    reduced={}
    for i in reversed(range(126)):
        if i in fi:reduced[i]={fi[i]:F(1)}
        else:
            v={}
            for j,a in piv[i].items():
                if j!=i:v=add(v,reduced[j],-a)
            reduced[i]=v
    Pi=[reduced[i] for i in range(126)]
    assert all(not mv(Pi,v) for v in fixed+allcols)
    A=[[parse(v) for v in M] for M in c['C2_coordinate_actions9_sparse_columns']]
    eta=res['physical_target3_basis_positions'];assert eta==[0,5,10]
    def mon_action(m,v):
        for i,n in enumerate(m):
            for _ in range(n):v=mv(A[i],v)
        return v
    monomials=[]
    for total in range(1,5):
        for ell in range(2):
            for inds in itertools.combinations_with_replacement(range(8),total-ell):
                m=[0]*9
                for i in inds:m[i]+=1
                m[8]=ell;monomials.append(tuple(m))
    assert len(monomials)==659
    library=[];ordercounts={};active_monomials=set()
    for s,pos in enumerate(eta):
        for m in monomials:
            if m[8]:
                v=mon_action(m[:8],{pos:F(1)});gamma=v
            else:
                i=max(i for i,n in enumerate(m[:8]) if n);alpha=list(m[:8]);alpha[i]-=1
                v=mon_action(alpha,{pos:F(1)});gamma={14+14*i+k:a for k,a in v.items()}
            value=mv(Pi,gamma)
            if value:
                library.append({'physical_target':s,'monomial9':list(m),'conormal26_value':srec(value)})
                active_monomials.add(m);ordercounts[sum(m)]=ordercounts.get(sum(m),0)+1
    # Literal V: first48 actual normalizer columns followed by -U9[1:].
    V=[[(x,F(1)) for x in column[1:]] for column in chain['exact_existing_Q48_C4']]
    V += [[(x,F(-1)) for x in column[1:]] for column in chain['exact_existing_U9_C4'][1:]]
    assert len(V)==56 and all(len(v)==3 for v in V)
    assert chain['exact_existing_U9_C4'][0]==[1,0,0,0]
    ops=chain['operator_scalar_nodes'];assert ops[0][0]=='c' or str(ops[0]).lower().find('zero')>=0
    requests={}
    for j,col in enumerate(V):
        for v in library:
            root,sign=col[v['physical_target']]
            if root==0:continue
            key=(j,root,tuple(v['monomial9']))
            value=scale(parse(v['conormal26_value']),sign)
            requests[key]=add(requests.get(key,{}),value)
    requests={k:v for k,v in requests.items() if v}
    roots=sorted({root for j,root,m in requests})
    total=sum(map(len,Pi))+sum(len(v['conormal26_value']) for v in library)+sum(map(len,requests.values()))
    assert total+8660+659*9+126*26<150000
    proof={'basis':'B=(k[ell]/ell²)[[zA0,...,zA6,rA]]; x orderz0..6,rA, ell last',
      'notation_clarification':'the earlier FIRST_BETTI note shorthand zA7 means the seven coordinates zA0,...,zA6, not an extra variable',
      'projection':'F1(0)/(span11a0+image d2(0)) -> k26; actual conormal J/mJ',
      'ordered_source_Hadamard':'f(X)-f(A)=sum_i (XiI-Ai)*Hi(f), evaluating prefixes0,...,i-1 atAi and leaving suffixformal',
      'no_extra_order_gap':'constant Hi uses A^alpha-ei on scalar Taylor degree|alpha|; C2 m4=0 kills |alpha|>=5. Periodic block usesell coefficient with|alpha|<=3, totaldegree<=4. Exactly one divided-difference order is included.',
      'horizontal_gamma':'for ell0 monomialalpha, i=maxsupport(alpha), sourceblocki vector A^(alpha-ei)*eta_s',
      'periodic_gamma':'for ell1 monomialalpha, periodic sourceblock vector A^alpha*eta_s',
      'regular_V_assumption':'literal fullR0 Q48 and-U9[1:] regular germs restrictedell²; their origin values vanish sincephysical3targetisminimal',
      'all_callbacks_remain_symbolic':True,
      'physical_implicit_and_honest_division_orders':'bound4 applies to REGULAR V scalar germs. A coefficient of an honest r/g function may require higher numerator data; no such coefficient is inferred from finiteC2values.',
      'still_required':'source-specific ancestor/coherentcopy expansion and existing-cache comparison of requested scalar Taylorcallbacks before evaluatingactual26x56 comparisonmatrix',
      'no_actual_m0_n_h_T2_or_H0_value':True}
    evidence={'fixed100_basis_sources':{'anchors11':anchors,'d2_selected89':selected},'fixed100_elimination_receipts':receipts,
      'Pi26_source126_columns':[srec(v) for v in Pi],'conormal_standard_source_indices26':free,
      'all518_d2_and11_anchor_vectors_annihilated':True,'scalar_monomial_adjoint_library':library,
      'actual_V56_source_columns':[[[root,str(sign)] for root,sign in col] for col in V],
      'actual_V_scalar_function_reference_artifact':'same allorder chain scalar nodes, without ell0 wrapper',
      'requested_actual_V_source_coefficients':[
        {'V_column':j,'chain_scalar_root':root,'monomial9':list(m),'conormal26_weight':srec(v)} for (j,root,m),v in sorted(requests.items())],
      'active_chain_scalar_roots':roots,'monomial_order_counts':ordercounts,'source_proof':proof}
    ep=out/'actual26-conormal-and-symbolic-V-Hadamard-source-adjoints.json';ep.write_text(json.dumps(evidence,separators=(',',':'))+'\n')
    report={'id':STEM,'all_passed':True,'actual_C2_conormal_dimension':26,'regular_scalar_Taylor_order_bound':4,
      'nonzero_scalar_monomial_adjoint_patterns':len(library),'requested_V_scalar_callback_count':len(requests),
      'active_chain_scalar_roots':len(roots),'active_orders':ordercounts,'named_QQ_slots':total,
      'actual_V_source_coefficients_evaluated':False,'input_sha256':hashes,'code_sha256':sha(Path(__file__)),
      'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start}
    (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps({k:report[k] for k in ['all_passed','nonzero_scalar_monomial_adjoint_patterns','requested_V_scalar_callback_count','active_chain_scalar_roots','active_orders','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
