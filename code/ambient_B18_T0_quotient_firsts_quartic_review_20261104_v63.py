#!/usr/bin/env python3
"""Independent fixed-certificate source4/opaque-sensitivity audit.

No producer worker, source solve, new jet, point or GB is run.  Exact sparse
bivariate polynomials interpret the stored source3/source4 substitutions.
"""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_T0_cokernel_Z_review_20261104_v63 import (
    add as dadd,mul as dmul,inv as dinv,scale as dscale,
    quotient_first,DescendingSpan)

ROOT=Path(__file__).resolve().parents[1]
AUTHOR='AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63'
PRE='AMBIENT_B18_T0_FIRSTJET_SENSITIVITY_PREFLIGHT_20261104_v63'
ID=AUTHOR+'_REVIEW'
ZERO=(0,0)
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pconst(c):return {ZERO:Q(c)} if c else {}
def padd(a,b):
    out=dict(a)
    for m,x in b.items():
        out[m]=out.get(m,Q(0))+x
        if not out[m]:del out[m]
    return out
def pscale(a,c):return {m:x*c for m,x in a.items() if x*c}
def pmul(a,b):
    if not a or not b:return {}
    if set(a)=={ZERO}:return pscale(b,a[ZERO])
    if set(b)=={ZERO}:return pscale(a,b[ZERO])
    out={}
    for m,x in a.items():
        for n,y in b.items():
            k=(m[0]+n[0],m[1]+n[1])
            if sum(k)<=4:out[k]=out.get(k,Q(0))+x*y
    return {m:x for m,x in out.items() if x}
def pinv(a):
    c=a.get(ZERO,Q(0));assert c
    h=pscale({m:x for m,x in a.items() if m!=ZERO},-1/c)
    ans=pconst(1);power=pconst(1)
    for _ in range(4):power=pmul(power,h);ans=padd(ans,power)
    return pscale(ans,1/c)
def poly(values,mon):return {tuple(m):Q(x) for m,x in zip(mon,values) if Q(x)}
def closure(nodes,roots,stops=frozenset()):
    seen=set();todo=list(roots)
    while todo:
        i=todo.pop()
        if i in seen:continue
        seen.add(i)
        if i in stops:continue
        op,*a=nodes[i]
        if op in ('a','m'):todo.extend(a)
        elif op in ('n','i'):todo.append(a[0])
        else:assert op in ('c','v'),(i,nodes[i])
    return sorted(seen)
def peval(nodes,needed,binding):
    out={}
    for i in needed:
        op,*a=nodes[i]
        if op=='c':v=pconst(Q(a[0],a[1]))
        elif op=='v':v=binding[i]
        elif op=='a':v=padd(out[a[0]],out[a[1]])
        elif op=='n':v=pscale(out[a[0]],Q(-1))
        elif op=='m':v=pmul(out[a[0]],out[a[1]])
        else:assert op=='i';v=pinv(out[a[0]])
        out[i]=v
    return out
def exact_det(A):
    A=[list(r) for r in A];v=Q(1)
    for k in range(len(A)):
        p=next(i for i in range(k,len(A)) if A[i][k])
        if p!=k:A[k],A[p]=A[p],A[k];v=-v
        c=A[k][k];v*=c
        for i in range(k+1,len(A)):
            q=A[i][k]/c
            for j in range(k+1,len(A)):A[i][j]-=q*A[k][j]
            A[i][k]=0
    return v

def main(out):
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/AUTHOR;report=read(folder/'report.json')
    prebase=ROOT/'reports'/PRE;pr=read(prebase/'report.json')
    hashes={}
    for folder_,r in ((folder,report),(prebase,pr)):
        assert r['all_passed']
        for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h;hashes[p]=h
        for p,h in r['evidence_sha256'].items():assert sha(folder_/p)==h;hashes[str((folder_/p).relative_to(ROOT))]=h
        hashes[str((folder_/'report.json').relative_to(ROOT))]=sha(folder_/'report.json')
    gates=[
      'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json',
      'reports/AMBIENT_B18_T0_COKERNEL_Z_20261104_v63_REVIEW/report.json',
      'reports/AMBIENT_B18_T1_T2_QUADRATIC_PHI_20261104_v63_REVIEW/report.json',
      'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW/report.json',
      'reports/AMBIENT_B18_CB_G_DIVISION_ORIGINS_20261104_v63_REVIEW/report.json',
    ]
    for p in gates:assert read(ROOT/p)['all_passed'];hashes[p]=sha(ROOT/p)
    source=read(ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json')
    cache=read(ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json')
    ca3=read(ROOT/'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json')
    cb3=read(ROOT/'reports/AMBIENT_B18_CB_G_DIVISION_ORIGINS_20261104_v63/degree3-certificate.json')
    C0=read(ROOT/'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json')
    full=read(ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json')
    pre=read(prebase/'formal-sensitivities-and-source4-size.json')
    evidence=read(folder/'quotient-firsts-and-actual-D0-firstjets.json')
    mon=[(d-j,j) for d in range(5) for j in range(d+1)]
    assert evidence['monomials_bivariate_through4']==list(map(list,mon))
    assert ca3['monomials']==cb3['ordered_monomials']==list(map(list,mon[:10]))
    nodes=source['DAG_scalar_nodes'];firsts={};source_counts={}
    for name,raw in (('CA',dict(ca3['all105_unknown_jets'])),('CB',dict(cb3['all57_unknown_jets']))):
        stage=read(folder/f'{name}-quartic-stage.json')
        assert stage==evidence['families'][name]
        frame=cache['both57_source2'][name];unknown=frame['unknown57'];eq=frame['equations57']
        assert len(unknown)==len(eq)==57 and len(set(unknown))==57
        assert stage['actual_unknown57_order']==unknown and stage['actual_equations57']==eq
        assert eq==source['families'][name]['all48_zero_IFT_quadratic_coefficient_equations']+frame['F9_nodes']
        if name=='CB':
            assert cb3['single57_unknown_source_node_ids']==unknown
            assert cb3['original48_equation_nodes']+cb3['selected_F9_original_nodes']==eq
        g=stage['g_original_node'];broots=stage['four_numerator_original_nodes']
        assert g==frame['f63']
        assert broots==[source['families'][name]['complete_normalized_physical_w_coefficients'][0][r][1] for r in (64,68,73,74)]
        needed=closure(nodes,eq+[g]+broots)
        assert needed==stage['source_needed_node_ids']
        assert {i for i in needed if nodes[i][0]=='v'}<=set(unknown)
        before={i:poly(raw[i],mon[:10]) for i in unknown}
        assert [[i,[str(before[i].get(m,Q(0))) for m in mon[:10]]] for i in unknown]==stage['cached_source3_all57']
        assert all(len(v)==15 for i,v in stage['all57_source4_jets'])
        after={i:poly(v,mon) for i,v in stage['all57_source4_jets']}
        assert set(after)==set(unknown)
        for k,i in enumerate(unknown):
            assert {m:x for m,x in after[i].items() if sum(m)<4}==before[i]
            assert [after[i].get(m,Q(0)) for m in mon[10:]]==list(map(Q,stage['quartic_correction57_by5'][k]))
        b=peval(nodes,needed,before);a=peval(nodes,needed,after)
        residual=[]
        for k,node in enumerate(eq):
            target={} if k<55 else {(1,0) if k==55 else (0,1):Q(1)}
            assert {m:x for m,x in b[node].items() if sum(m)<4}==target
            assert a[node]==target,(name,k,a[node])
            residual.append([b[node].get(m,Q(0)) for m in mon[10:]])
        assert residual==[list(map(Q,v)) for v in stage['quartic_residual57_by5']]
        J=[list(map(Q,v)) for v in stage['fixed_J57']]
        assert stage['fixed_J57']==frame['fixed_J57'] and stage['fixed_J57_inverse']==frame['fixed_J57_inverse']
        correction=[list(map(Q,v)) for v in stage['quartic_correction57_by5']]
        for k in range(57):
            for j in range(5):
                assert sum((J[k][i]*correction[i][j] for i in range(57)),Q(0))==-residual[k][j]
        assert a[g]==poly(stage['g_all15_coefficients'],mon)
        g3={m:x for m,x in a[g].items() if sum(m)==3}
        assert len(g3) and all(sum(m)>=3 for m in a[g])
        assert g3[(3,0)]*g3[(0,3)]==Q(stage['endpoint2_by2_minor'])
        hs=[]
        for root,stored,expected in zip(broots,stage['four_numerator_all15_coefficients'],stage['quotient_firsts4_by2']):
            assert a[root]==poly(stored,mon)
            assert all(sum(m)==4 for m in a[root])
            h=(a[root].get((4,0),Q(0))/g3[(3,0)],a[root].get((0,4),Q(0))/g3[(0,3)])
            assert h==tuple(map(Q,expected))
            assert pmul(g3,poly(h,((1,0),(0,1))))==a[root]
            hs.append(h)
        firsts[name]=hs;source_counts[name]=len(needed)
    assert firsts['CA']==[(Q(3,8),Q(3,4))]+[(Q(0),Q(0))]*3
    assert firsts['CB']==[(Q(1),Q(-1,2))]+[(Q(0),Q(0))]*3

    # Full9 CB->CA conversion independently checked through common physical
    # coordinates, not by identifying the two source parameter pairs.
    vars={n[1]:i for i,n in enumerate(nodes) if n[0]=='v'}
    sols={name:{i:list(map(Q,v)) for i,v in cache['both57_source2'][name]['all57_source2jets']} for name in ('CA','CB')}
    physical=[vars[f'x{r}'] for r in source['physical_zero_base_rows9']]
    XA=[sols['CA'][i][1:10] for i in physical];XB=[sols['CB'][i][1:10] for i in physical]
    adapter=[list(map(Q,v)) for v in cache['YB_to_YA_linear']]
    assert exact_det(XB) and exact_det(adapter)==Q(1,4096)
    assert [[sum((XB[i][k]*adapter[k][j] for k in range(9)),Q(0)) for j in range(9)] for i in range(9)]==XA
    gradients=[[x*adapter[7][k]+y*adapter[8][k] for k in range(9)] for x,y in firsts['CB']]
    gradients += [[x*Q(k==7)+y*Q(k==8) for k in range(9)] for x,y in firsts['CA']]
    assert gradients==[list(map(Q,v)) for v in evidence['all8_regular_function_first_gradients_in_FULL9_CA_F9']]
    differences=[[gradients[i][k]-gradients[4+i][k] for k in range(9)] for i in range(4)]
    assert not any(x for a in differences for x in a)
    assert differences==[list(map(Q,v)) for v in evidence['four_oldCB_minus_newCA_first_gradients_in_FULL9']]

    # Independently rebuild C0's COMPLETE640 firstjets with eight distinct
    # opaque positive-order function symbols. No g derivative is set to zero.
    allnodes=C0['DAG_scalar_nodes'];N=len(nodes);values={}
    for name in ('CA','CB'):
        for i,v in sols[name].items():
            d=v[1:10]
            if name=='CB':d=[sum((d[k]*adapter[k][j] for k in range(9)),Q(0)) for j in range(9)]
            q=(Q(0),{k:x for k,x in enumerate(d) if x})
            if i in values:assert values[i]==q
            values[i]=q
    for i,n in enumerate(nodes):
        op,*ar=n
        if op=='v':assert i in values
        elif op=='c':values[i]=(Q(ar[0],ar[1]),{})
        elif op=='a':values[i]=dadd(values[ar[0]],values[ar[1]])
        elif op=='m':values[i]=dmul(values[ar[0]],values[ar[1]])
        elif op=='n':values[i]=dscale(values[ar[0]],Q(-1))
        else:assert op=='i';values[i]=dinv(values[ar[0]])
    special={int(i):(Q(v[0]),{k:Q(x) for k,x in enumerate(v[1:]) if Q(x)}) for key in ('all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA') for i,v in cache[key].items()}
    opaque=[q['node'] for q in full['all4_known_regular_g_divisions']]+[q['node'] for q in C0['all4_known_regular_gA_divisions']]
    assert evidence['original_regular_function_nodes_unmodified']=={'old_CB':opaque[:4],'new_CA':opaque[4:]}
    for k,i in enumerate(opaque):special[i]=(Q(0),{9+k:Q(1)})
    special[full['CB_g_restriction_node']]=(Q(0),{})
    special[C0['actual_C0_basis8_scalar_nodes'][-1]]=(Q(0),{})
    mon2=list(map(tuple,cache['source2_monomials']))
    for rec in C0['all35_regular_coordinate_divisions']:
        p=cache['both57_source2']['CA']['pure_function2jets'][str(rec['original_pure_CA_source_function'])]
        for k,q in enumerate(rec['all7_regular_coordinate_divisions']):
            h=quotient_first(p,mon2,k,True);assert h[0]==Q(q['origin']);special[q['node']]=h
    relationcols=C0['honest_D0_relation_columns'];assert len(relationcols)==640
    needed=closure(allnodes,[i for col in relationcols for i in col],set(range(N))|set(special))
    for i in needed:
        if i<N:continue
        if i in special:values[i]=special[i];continue
        op,*ar=allnodes[i]
        if op=='c':q=(Q(ar[0],ar[1]),{})
        elif op=='a':q=dadd(values[ar[0]],values[ar[1]])
        elif op=='m':q=dmul(values[ar[0]],values[ar[1]])
        elif op=='n':q=dscale(values[ar[0]],Q(-1))
        else:assert op=='i';q=dinv(values[ar[0]])
        values[i]=q
    records=[];matrix=[];fixed=DescendingSpan();enlarged=DescendingSpan()
    for j,col in enumerate(relationcols):
        assert all(values[i][0]==0 for i in col)
        known=[values[i][1].get(k,Q(0)) for i in col for k in range(9)]
        assert known==list(map(Q,pre['all640_known_firstjet81_columns'][j]))
        sens=[[values[i][1].get(9+q,Q(0)) for i in col] for q in range(8)]
        if not any(x for a in sens for x in a):fixed.insert(known)
        else:records.append((j,sens))
        actual=[values[i][1].get(k,Q(0))+sum((values[i][1].get(9+q,Q(0))*gradients[q][k] for q in range(8)),Q(0)) for i in col for k in range(9)]
        matrix.append(actual);enlarged.insert(known)
        for q in range(8):
            for t in (7,8):
                direction=adapter[t] if q<4 else [Q(k==t) for k in range(9)]
                enlarged.insert([c*direction[k] for c in sens[q] for k in range(9)])
    assert len(records)==4
    for q,(index,sens) in enumerate(records):
        assert index==604+9*q
        expect=[[Q(0)]*9 for _ in range(8)];expect[q][8]=1;expect[q+4][8]=-1
        assert sens==expect
        stored=pre['all_nonzero_opaque_function_sensitivities'][q]
        assert stored['relation']==index and sens==[list(map(Q,a)) for a in stored['sensitivities8_by9']]
    assert len(fixed.rows)==34==pre['fixed_relation_firstspan_rank']
    assert len(enlarged.rows)==40==pre['max_enlarged_span_rank']
    assert matrix==[list(map(Q,v)) for v in evidence['all640_exact_D0_firstjet81_columns']]
    actualspace=DescendingSpan()
    for col in matrix:actualspace.insert(col)
    assert len(actualspace.rows)==38==evidence['exact_firstspan_rank']
    recorded=DescendingSpan()
    for p,ar in evidence['exact_firstspan_basis']:
        v={k:Q(x) for k,x in ar};assert not actualspace.reduce(v);assert recorded.insert(v)
    assert len(recorded.rows)==38
    survivors=DescendingSpan()
    for k,test in enumerate(evidence['visible_gA_class_first_multiplication_tests']):
        v=[Q(0)]*81;v[72+k]=1;r=actualspace.reduce(v)
        assert bool(r)==(k>=7) and test['gA_basis_class_multiple_in_linear_relation_span']==(not bool(r))
        survivors.insert(r)
    assert len(survivors.rows)==2
    proof=dict(both_cached_source3_frames_bound_no_refit=True,
      independent_sparse_QQ15_jet_interpreter=True,actual_ancestor_counts=source_counts,
      all114_equations_all15_coefficients_verified=True,
      all570_saved_quartic_corrections_and_linear_residual_identities_verified=True,
      all8_regular_quotient_firsts_verified_all5_quartic_coefficients=True,
      CA_H64_first=['3/8','3/4'],CB_G64_first=['1','-1/2'],other6_firsts_zero=True,
      full9_physical_adapter_checked_and_all4_gradient_differences_zero=True,
      all640_dual_sensitivity_columns_independently_reconstructed=True,
      exactly4_unknown_GminusH_sensitivities_in_gA_row=True,
      all8_opaque_function_derivatives_retained_not_assumed_zero=True,
      exact_D0_firstspan_rank=38,visible_gA_zero7_multiples_killed=True,
      visible_gA_r_s_multiples_independently_nonzero=True,
      no_all_order_lift_with_nonzero_visible_gA_leading_class_can_be_socle=True,
      hidden_finite_socle_or_depth_not_determined=True,
      no_author_worker_source_solve_new_jet_point_prime_or_GB=True)
    (out/'independent-quartic-and-sensitivity-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    helper=ROOT/'code/ambient_B18_T0_cokernel_Z_review_20261104_v63.py'
    for p in (Path(__file__),helper,ROOT/f'research/notes/{ID}.md'):hashes[str(p.relative_to(ROOT))]=sha(p)
    review=dict(id=ID,all_passed=True,input_sha256=hashes,
      evidence_sha256={'independent-quartic-and-sensitivity-proof.json':sha(out/'independent-quartic-and-sensitivity-proof.json')},
      source4_quotient_firsts_and_D0_firstspan38_independently_accepted=True,
      visible_leading_class_not_socle_any_lift_hidden_socle_open=True,
      seconds=time.monotonic()-started,max_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(review,indent=2)+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    print(json.dumps(review))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
