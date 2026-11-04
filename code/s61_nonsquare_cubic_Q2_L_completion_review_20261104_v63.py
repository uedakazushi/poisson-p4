#!/usr/bin/env python3
"""Independent actual three-dot, GL3 and full lazy transport audit, no solver."""
import argparse,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from s61_nonsquare_cubic_zw25_pair_completion_review_20261104_v63 import read,sha,plus,times,neg,sub,var,eye,total,mm,ZERO,ONE
ROOT=Path(__file__).resolve().parents[1]
BASE='S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63';ID=BASE+'_REVIEW'
def constant(v):return {():Q(v)} if v else ZERO

def universal():
    b0,b1,b2,v=[var('universal',i) for i in range(4)];c=constant
    U0=[[v,ONE,ZERO],[c(-18),ZERO,ONE],[c(-9),ZERO,ZERO]]
    Ui=[[ZERO,ZERO,c(Q(-1,9))],[ONE,ZERO,times(v,c(Q(1,9)))],[ZERO,ONE,c(-2)]]
    E=[[ONE,neg(b0),neg(b1)],[ZERO,ONE,ZERO],[ZERO,ZERO,ONE]];Ei=[[ONE,b0,b1],[ZERO,ONE,ZERO],[ZERO,ZERO,ONE]]
    assert mm(U0,Ui)==mm(Ui,U0)==eye(3) and mm(E,Ei)==mm(Ei,E)==eye(3)
    U=mm(U0,E);V=mm(Ei,Ui);assert mm(U,V)==mm(V,U)==eye(3)
    t=plus(times(b0,v),plus(times(c(-18),b1),times(c(-9),b2)))
    assert mm([[b0,b1,b2]],U)[0]==[t,times(sub(ONE,t),b0),times(sub(ONE,t),b1)]
    determinant=ZERO
    for p in itertools.permutations(range(3)):
        x=ONE
        for i in range(3):x=times(x,U0[i][p[i]])
        determinant=plus(determinant,neg(x) if sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))%2 else x)
    assert determinant==c(-9)
    return dict(both_inverse_products_U0_E_U=True,all3_row_output_identities=True,det_U0='-9',independent_free_polynomial_engine=True)

def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/BASE;pr=read(folder/'report.json');cert=read(folder/'certificate.json');assert pr['frozen'] and sha(folder/'certificate.json')==pr['certificate_sha256'];hashes={}
    for p,h in pr['input_sha256'].items():assert sha(ROOT/p)==h;hashes[p]=h
    parentfolder=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_UNIT_COMPLETION_20261104_v63';parent=read(parentfolder/'certificate.json');oldgate=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_UNIT_COMPLETION_20261104_v63_REVIEW/report.json';og=read(oldgate);assert og['all_passed']
    assert og['input_sha256'][str((parentfolder/'certificate.json').relative_to(ROOT))]==sha(parentfolder/'certificate.json');hashes[str(oldgate.relative_to(ROOT))]=sha(oldgate)
    cp=read(ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json');minorfolder=ROOT/'reports/S61_NONSQUARE_CUBIC_Q2_L_MINORS_20261104_v63';mc=read(minorfolder/'certificate.json');mr=read(minorfolder/'report.json');assert sha(minorfolder/'certificate.json')==mr['certificate_sha256']
    for p,h in mr['input_sha256'].items():assert sha(ROOT/p)==h;hashes[p]=h
    R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field()
    def poly(p):return R.from_dict({tuple(m):s.QQ(c) for m,c in p})
    def decode(q):return K(poly(q[0]))/K(poly(q[1]))
    rows=parent['new16_source_original25_rows'];cols=parent['new23_target_original32_columns'];Aroots=parent['remaining16x23_scalar_roots'];group=[21,22,25];pnodes=parent['scalar_DAG']['nodes']
    block=[[decode(cp['matrix25x32'][i][j]) for j in group] for i in [5,9]]
    for r,i in enumerate([5,9]):
        for j,col in enumerate(group):assert pnodes[Aroots[rows.index(i)][cols.index(col)]]==['checkpoint',i,col]
    assert block==[[decode(q) for q in row] for row in mc['actual2x3']]
    # Three already fixed minors; no search and no assertion about their ideal.
    named=[a,d,500*a+90*e**2+3*e-24,(50*a+18*e**2-9*e-2)*d+60*e-40,100*a+18*e**2-9*e-2,5*a*d-6*e-4]
    for outcome,(j,k) in zip(mc['exact_three_minor_outcomes'],itertools.combinations(range(3),2)):
        delta=block[0][j]*block[1][k]-block[0][k]*block[1][j];assert delta and delta==decode(outcome['minor']);assert outcome['original_target_column_pair']==[group[j],group[k]]
        residual=poly(outcome['untouched_unknown_residual_numerator']);product=residual
        for f,n in zip(named,outcome['named_unit_numerator_factor_powers_a_d_F_B_C_J']):assert n>=0;product*=f**n
        assert delta==K(product)/K.from_expr(s.sympify(outcome['nonzero_QQ_denominator']))
        assert residual and not residual.is_ground and all(residual.rem(f) for f in named)
        assert not outcome['named_unit_verdict']
    b=block[0];T=900*a-126*e-121;dual=[K(6*(18*e+7))/K(5*a),K(-18),K(-9)]
    assert b[0]==K(a*T)/150 and 2*b[1]+b[2]==K((18*e+7)*T)/1125-K.one/9
    assert sum((b[i]*dual[i] for i in range(3)),K.zero)==1
    actual=cert['actual_three_dot'];assert actual['old25_source_row']==5 and actual['old32_target_columns']==group
    assert b==list(map(decode,actual['coefficients'])) and dual==list(map(decode,actual['dual'])) and decode(actual['dot'])==1
    assert dual[0].denom==5*a and dual[1].denom==dual[2].denom==R.one
    univ=universal()
    nodes=cert['scalar_DAG']['nodes'];assert nodes[:len(pnodes)]==pnodes and cert['scalar_DAG']['parent_nodes_unchanged_prefix_length']==len(pnodes)==2138;assert len(nodes)==3914
    atoms=list(map(decode,cert['scalar_DAG']['new_unit_atoms_num_den']));assert atoms==dual+[K(-1)/9,K.one/9,K(-2)]
    values=[]
    for i,n in enumerate(nodes):
        o,*t=n
        if o=='z':q=ZERO
        elif o=='one':q=ONE
        elif o=='checkpoint':assert cp['matrix25x32'][t[0]][t[1]][0];q=var('checkpoint',*t)
        elif o=='accepted_pair_second':q=var('dual',t[0])
        elif o=='new_unit_atom':q=var('vGL3') if t[0]==0 else constant(Q(str(atoms[t[0]])))
        else:
            assert o in ['a','s','m'] and all(j<i for j in t)
            q=plus(values[t[0]],values[t[1]]) if o=='a' else sub(values[t[0]],values[t[1]]) if o=='s' else times(values[t[0]],values[t[1]])
        values.append(q)
    def get(A):return [[values[q] for q in row] for row in A]
    A=get(Aroots);b0,b1,b2=[A[0][cols.index(j)] for j in group];v=var('vGL3');c=constant
    U0=[[v,ONE,ZERO],[c(-18),ZERO,ONE],[c(-9),ZERO,ZERO]];Ui=[[ZERO,ZERO,c(Q(-1,9))],[ONE,ZERO,times(v,c(Q(1,9)))],[ZERO,ONE,c(-2)]]
    E=eye(3);Ei=eye(3);E[0][1]=neg(b0);E[0][2]=neg(b1);Ei[0][1]=b0;Ei[0][2]=b1
    U=mm(U0,E);V=mm(Ei,Ui)
    for key,M in [('U0_roots',U0),('U0_inverse_roots',Ui),('E_roots',E),('E_inverse_roots',Ei),('U3_roots',U),('U3_inverse_roots',V)]:assert get(cert[key])==M
    target=eye(23);targeti=eye(23);pos=[cols.index(j) for j in group]
    for i in range(3):
        for j in range(3):target[pos[i]][pos[j]]=U[i][j];targeti[pos[i]][pos[j]]=V[i][j]
    outside=[j for j in cols if j not in group and Aroots[0][cols.index(j)]];assert len(outside)==15 and cert['tail_clear_original32_columns']==outside
    tail=eye(23);taili=eye(23);pivot=pos[0]
    for j in outside:tail[pivot][cols.index(j)]=neg(A[0][cols.index(j)]);taili[pivot][cols.index(j)]=A[0][cols.index(j)]
    assert mm(tail,taili)==mm(taili,tail)==eye(23)
    target=mm(target,tail);targeti=mm(taili,targeti);assert get(cert['target23_roots'])==target and get(cert['target23_inverse_roots'])==targeti
    AT=mm(A,target);assert [values[i] for i in cert['pivot_full23_row_roots_before_exact_simplification']]==AT[0]
    source=eye(16)
    for i in range(1,16):source[i][0]=neg(AT[i][pivot])
    assert get(cert['source16_roots'])==source
    newcols=[j for j in cols if j!=21];expected=[[AT[i][cols.index(j)] for j in newcols] for i in range(1,16)]
    assert get(cert['remaining15x22_scalar_roots'])==expected
    assert cert['old16_source_original25_rows']==rows and cert['old23_target_original32_columns']==cols
    assert cert['new15_source_original25_labels']==rows[1:] and cert['new22_target_original32_columns']==newcols
    assert cert['pivot_row_exact_simplified']==[int(i==pivot) for i in range(23)]
    for j in [27,28,29,30,31]:
        k=cols.index(j);assert not A[0][k]
        for i in range(23):assert target[i][k]==targeti[i][k]==(ONE if i==k else ZERO)
    assert cert['source_change_determinant']=='1' and cert['target_change_determinant']=='-9' and cert['total_unit_pivots_in_transformed_presentation']==10
    assert not any(n[0] in ['inv','i'] for n in nodes)
    proof=dict(old2138_prefix_unchanged=True,all3914_free_polynomial_nodes_replayed=True,all330_final_entries_bound=True,actual_three_dot_source_checkpoint_bound=True,visible_T_identities_QQ=True,new_denominators_only_QQ_and_a=True,universal=univ,all15_tail_clearings_and_source16_bound=True,all5_heavy_targets_unchanged=True,actual_same_cokernel_shapes=[[16,23],[15,22]],total_transformed_unit_pivots=10,three_fixed_minor_frontier_reconstructed=True,no_named_unit_minor_does_not_imply_rankdrop_or_proper_ideal=True,scope='One global GL3 completion only, no wholeFULLU57/full239 or kernel count change')
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in [Path(__file__),ROOT/'code/s61_nonsquare_cubic_zw25_pair_completion_review_20261104_v63.py',ROOT/f'research/notes/{ID}.md',folder/'report.json',folder/'certificate.json']:hashes[str(p.relative_to(ROOT))]=sha(p)
    r=dict(id=ID,all_passed=True,actual_GL3_unit_completion_and_same_cokernel_verified=True,all3914_DAG_nodes_and330_remaining_entries_verified=True,no_author_solver_search_newpoint_prime_GB_or_unknown_inverse=True,input_sha256=hashes,evidence_sha256={'independent-proof.json':sha(out/'independent-proof.json')},seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=r['seconds'],rss=r['max_rss_kib'],report_sha256=sha(out/'report.json'))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
