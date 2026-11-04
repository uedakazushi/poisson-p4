#!/usr/bin/env python3
"""ONE necessary two-family source4 correction, only four quotient firsts.

Reuses exact saved source3 and fixed57 inverses. No CB48 refit, no source>4,
full honest2 matrix, point, prime, GB, or mutation of old regular functions.
"""
import argparse
from fractions import Fraction as F
import hashlib,json
from pathlib import Path
import resource,signal,time
from ambient_B18_top_graph_telescoping_review_20261104_v63 import closure
from ambient_B18_T3_zero_structural_20261104_v63 import insert_basis
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63'
MON=[(d-i,i) for d in range(5) for i in range(d+1)]
POS={m:i for i,m in enumerate(MON)}
CONV=[[(i,j) for i,a in enumerate(MON) for j,b in enumerate(MON) if (a[0]+b[0],a[1]+b[1])==m] for m in MON]
ZERO=(F(0),)*15
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def jm(a,b):
    if not any(a[1:]):return tuple(a[0]*x for x in b)
    if not any(b[1:]):return tuple(b[0]*x for x in a)
    return tuple(sum((a[i]*b[j] for i,j in pairs if a[i] and b[j]),F(0)) for pairs in CONV)
def ji(a):
    assert a[0];v=[1/a[0]]
    for k in range(1,15):v.append(-sum((a[i]*v[j] for i,j in CONV[k] if i and a[i] and v[j]),F(0))/a[0])
    return tuple(v)
def evaluate(nodes,needed,bind):
    v={}
    for i in needed:
        n=nodes[i];op=n[0]
        if op=='c':q=(F(n[1],n[2]),)+(F(0),)*14
        elif op=='v':q=bind[i]
        elif op=='a':q=tuple(x+y for x,y in zip(v[n[1]],v[n[2]]))
        elif op=='n':q=tuple(-x for x in v[n[1]])
        elif op=='m':q=jm(v[n[1]],v[n[2]])
        elif op=='i':q=ji(v[n[1]])
        else:raise ValueError(op)
        v[i]=q
    return v
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    jp=ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json'
    ap=ROOT/'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json'
    bp=ROOT/'reports/AMBIENT_B18_CB_G_DIVISION_ORIGINS_20261104_v63/degree3-certificate.json'
    cp=ROOT/'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json'
    pp=ROOT/'reports/AMBIENT_B18_T0_FIRSTJET_SENSITIVITY_PREFLIGHT_20261104_v63/formal-sensitivities-and-source4-size.json'
    fp=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json'
    source,jet,CA3,CB3,C0,pre,full=map(read,[sp,jp,ap,bp,cp,pp,fp]);nodes=source['DAG_scalar_nodes']
    assert CA3['monomials']==CB3['ordered_monomials']==list(map(list,MON[:10]))
    assert CB3['single57_unknown_source_node_ids']==jet['both57_source2']['CB']['unknown57']
    assert CB3['original48_equation_nodes']+CB3['selected_F9_original_nodes']==jet['both57_source2']['CB']['equations57']
    actual={};overrides={};quartic_columns=range(10,15)
    for name,cache in [('CA',dict(CA3['all105_unknown_jets'])),('CB',dict(CB3['all57_unknown_jets']))]:
        frame=jet['both57_source2'][name];u=frame['unknown57'];E=frame['equations57']
        assert set(u)<=set(cache)
        fixedinv=[list(map(F,r)) for r in frame['fixed_J57_inverse']]
        g=frame['f63'];broots=[source['families'][name]['complete_normalized_physical_w_coefficients'][0][row][1] for row in [64,68,73,74]]
        roots=E+[g]+broots;needed=sorted(closure(nodes,roots))
        assert {i for i in needed if nodes[i][0]=='v'}<=set(u)
        sol={i:tuple(map(F,cache[i]))+(F(0),)*5 for i in u}
        before=evaluate(nodes,needed,sol)
        target=[]
        for k in range(57):
            a=[F(0)]*15
            if k==55:a[1]=1
            if k==56:a[2]=1
            target.append(tuple(a))
        assert all(before[node][:10]==target[k][:10] for k,node in enumerate(E))
        residual=[[before[node][j]-target[k][j] for j in quartic_columns] for k,node in enumerate(E)]
        correction=[[-sum((fixedinv[i][k]*residual[k][j] for k in range(57)),F(0)) for j in range(5)] for i in range(57)]
        sol4={node:sol[node][:10]+tuple(correction[k]) for k,node in enumerate(u)}
        after=evaluate(nodes,needed,sol4)
        assert all(after[node]==target[k] for k,node in enumerate(E))
        gg=after[g];assert all(x==0 for x in gg[:6]) and any(gg[6:10])
        if name=='CA':assert list(map(str,gg[:10]))==C0['cached_actual_CA_gA_all_degree3_coefficients']
        else:assert list(map(str,gg[:10]))==CB3['g_CB_all10_coefficients']
        mat=[[gg[6],F(0)],[gg[7],gg[6]],[gg[8],gg[7]],[gg[9],gg[8]],[F(0),gg[9]]]
        assert gg[6] and gg[9]
        first=[];nums=[]
        for b in broots:
            bb=after[b];assert all(x==0 for x in bb[:10])
            h=[bb[10]/gg[6],bb[14]/gg[9]]
            assert all(sum((a*v for a,v in zip(row,h)),F(0))==value for row,value in zip(mat,bb[10:]))
            first.append(h);nums.append(bb)
        actual[name]={'actual_unknown57_order':u,'actual_equations57':E,'source_needed_node_ids':needed,
          'fixed_J57':frame['fixed_J57'],'fixed_J57_inverse':frame['fixed_J57_inverse'],'fixed_J57_determinant':frame['determinant'],
          'cached_source3_all57':[[i,list(map(str,sol[i][:10]))] for i in u],
          'quartic_residual57_by5':[[str(x) for x in a] for a in residual],
          'quartic_correction57_by5':[[str(x) for x in a] for a in correction],
          'all57_source4_jets':[[i,list(map(str,sol4[i]))] for i in u],
          'all57_equations_all15_coefficients_verified':True,
          'g_original_node':g,'g_all15_coefficients':list(map(str,gg)),
          'four_numerator_original_nodes':broots,'four_numerator_all15_coefficients':[list(map(str,a)) for a in nums],
          'quotient_firsts4_by2':[list(map(str,a)) for a in first],
          'endpoint2_by2_minor':str(gg[6]*gg[9]),'all5_coefficients_each_quotient_identity_verified':True}
        overrides[name]=first
        checkpoint=out/f'{name}-quartic-stage.json';checkpoint.write_text(json.dumps(actual[name],separators=(',',':'))+'\n')
    adapter=[list(map(F,r)) for r in jet['YB_to_YA_linear']]
    grads=[];diff=[]
    for a,b in overrides['CB']:
        grads.append([a*adapter[7][k]+b*adapter[8][k] for k in range(9)])
    for a,b in overrides['CA']:grads.append([a*F(k==7)+b*F(k==8) for k in range(9)])
    for i in range(4):diff.append([grads[i][k]-grads[4+i][k] for k in range(9)])
    matrix=[list(map(F,a)) for a in pre['all640_known_firstjet81_columns']]
    for rec in pre['all_nonzero_opaque_function_sensitivities']:
        col=matrix[rec['relation']];sens=[list(map(F,a)) for a in rec['sensitivities8_by9']]
        for row in range(9):
            for k in range(9):col[9*row+k]+=sum((sens[i][row]*grads[i][k] for i in range(8)),F(0))
    basis=[];pivots=[]
    for i,a in enumerate(matrix):
        if insert_basis(basis,a) is not None:pivots.append(i)
    def reduce(a):
        a={i:x for i,x in enumerate(a) if x}
        for p,b in basis:
            if p in a:
                c=a[p]
                for i,x in b.items():
                    v=a.get(i,F(0))-c*x
                    if v:a[i]=v
                    else:a.pop(i,None)
        return a
    survivor_basis=[];class_tests=[]
    for k in range(9):
        a=[F(0)]*81;a[72+k]=1;r=reduce(a)
        insert_basis(survivor_basis,[r.get(i,F(0)) for i in range(81)])
        class_tests.append({'CA_F9_direction':k,'gA_basis_class_multiple_in_linear_relation_span':not bool(r),
                            'remainder':[[i,str(v)] for i,v in sorted(r.items())]})
    evidence={'monomials_bivariate_through4':list(map(list,MON)),'families':actual,
      'all8_regular_function_first_gradients_in_FULL9_CA_F9':[list(map(str,a)) for a in grads],
      'four_oldCB_minus_newCA_first_gradients_in_FULL9':[list(map(str,a)) for a in diff],
      'original_regular_function_nodes_unmodified':{'old_CB':[q['node'] for q in full['all4_known_regular_g_divisions']],
                                                   'new_CA':[q['node'] for q in C0['all4_known_regular_gA_divisions']]},
      'separate_exact_first_derivative_override_only':True,
      'all640_exact_D0_firstjet81_columns':[[str(x) for x in a] for a in matrix],
      'exact_firstspan_rank':len(basis),'independent_column_indices':pivots,
      'exact_firstspan_basis':[[p,[[i,str(x)] for i,x in sorted(a.items())]] for p,a in basis],
      'visible_gA_class_first_multiplication_tests':class_tests,
      'visible_class_multiplication_first_surviving_rank':len(survivor_basis),
      'no_socle_or_depth_from_firstjet_inferred':True}
    ep=out/'quotient-firsts-and-actual-D0-firstjets.json';ep.write_text(json.dumps(evidence,separators=(',',':'))+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',sp,jp,ap,bp,cp,pp,fp]
    report={'id':STEM,'all_passed':True,'status':'exact_two_family_necessary_quartic_contraction_and_D0_firstjets_only',
      'each_source57_equations_verified_all15_coefficients':True,
      'all8_quotient_firsts_verified_against_all5_quartic_coefficients':True,
      'D0_exact_firstspan_rank':len(basis),'visible_gA_class_first_multiplication_surviving_rank':len(survivor_basis),
      'actual_T0_nonzero_already_proved_elsewhere':True,'T0_socle_cyclicity_or_depth_claimed':False,
      'no_new_source_degree_above4_CB48_refit_full_Honest2_GB_points_or_primes':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start,
      'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['all_passed','D0_exact_firstspan_rank','visible_gA_class_first_multiplication_surviving_rank','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
