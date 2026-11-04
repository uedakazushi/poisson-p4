#!/usr/bin/env python3
"""Independent bottom-up cofactor DAG evaluation and saved4-step audit; no search."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def accepted(p):
    r=read(p);assert r['all_passed']
    for n,h in r.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n
    for n,h in r.get('evidence_sha256',{}).items():
        q=ROOT/n if (ROOT/n).exists() else p.parent/n
        assert sha(q)==h,n
    return r

def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    assert not out.exists();out.mkdir(parents=True)
    base=ROOT/'reports/S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63'
    pp=ROOT/'reports/S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63'
    parentreview=pp.with_name(pp.name+'_REVIEW')/'report.json';pr=accepted(parentreview)
    assert pr['all15_exact354_vectors_and35unit_deletions_match_all319_coefficients']
    assert pr['all11_final_vectors_have_literal_identity_quartic_block']
    qp=ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63_REVIEW/report.json';qr=accepted(qp)
    assert qr['pointwise_22_directions_independent_on_FULLU']
    ap=base/'report.json';author=accepted(ap);cp=base/'certificate.json';c=read(cp)
    sp=base/'scope-receipt.json';scope=read(sp)
    assert scope['all_passed'] and sha(cp)==author['certificate_sha256']==scope['raw_certificate_sha256']
    assert sha(ap)==scope['raw_probe_report_sha256']
    source=ROOT/'code/s61_nonsquare_bianchi_residual_four_20261104_v63.py'
    note=ROOT/'research/notes/S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63.md'
    assert sha(source)==scope['source_sha256'] and sha(note)==scope['note_sha256']
    pcp=pp/'certificate.json';pivpath=pp/'block-pivots.json';p=read(pcp);piv=read(pivpath)
    assert sha(pcp)==c['parent_certificate_sha256']
    R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field()
    parse=lambda v:K.from_expr(s.sympify(v))
    U=[[parse(v) for v in row] for row in piv['elementary_transform']]
    ids=[0,4,13,14]
    assert U[11:]==[[K(int(i==j)) for j in range(15)] for i in ids]
    assert c['residual_transform_rows15to4']==[[str(v) for v in row] for row in U[11:]]
    assert c['residual_original_block_rows']==ids
    assert c['residual_parent_piece_indices']==[p['raw15_source_lifts'][i]['source_piece_index'] for i in ids]
    definitions=p['actual319_source_definitions']
    cubic=[j for j,v in enumerate(definitions) if v[0]=='cubic_normalized']
    quartic=[j for j,v in enumerate(definitions) if v[0]=='quartic']
    assert len(cubic)==112 and len(quartic)==33
    assert cubic==c['cubic_source319_indices'] and quartic==c['quartic_source319_indices']
    originals=[p['raw15_source_lifts'][i]['after35_deletions_coefficients319'] for i in ids]
    cn=[[v[j] for j in cubic] for v in originals];qn=[[v[j] for j in quartic] for v in originals]
    assert cn==c['cubic_coefficient_nodes'] and qn==c['quartic_coefficient_nodes']
    A=p['arithmetic'];nodes=A['nodes'];weights=A['node_weights'];needed={0};stack=[n for row in cn+qn for n in row]
    # Different evaluator: collect closure first, then evaluate in topological order.
    while stack:
        n=stack.pop()
        if n in needed:continue
        needed.add(n);assert weights[n]==0
        nd=nodes[n]
        if nd[0]=='l':stack.extend(i for coeff,i in nd[1])
        elif nd[0]=='m':stack.extend(nd[1:])
        else:assert nd[0]=='c'
    acache={}
    def atom(i):
        if i not in acache:
            num,den=A['atoms'][i]
            poly=lambda terms:R.from_dict({tuple(ex):s.QQ(v) for ex,v in terms})
            acache[i]=K(poly(num))/K(poly(den))
        return acache[i]
    values={0:K.zero}
    for i in sorted(needed-{0}):
        nd=nodes[i]
        if nd[0]=='c':values[i]=atom(nd[1])
        elif nd[0]=='m':values[i]=values[nd[1]]*values[nd[2]]
        else:values[i]=sum((atom(k)*values[n] for k,n in nd[1]),K.zero)
    assert all(not values[n] for row in qn for n in row)
    M=[[values[n] for n in row] for row in cn]
    assert [[str(v) for v in row] for row in M]==c['all448_exact_cubic_values']
    assert sum(bool(v) for row in M for v in row)==119
    B=[row[:] for row in M];T=[[K(int(i==j)) for j in range(4)] for i in range(4)]
    known=[a,500*a+90*e**2+3*e-24,(50*a+18*e**2-9*e-2)*d+60*e-40,100*a+18*e**2-9*e-2]
    def unit(v):
        powers=[];ground=[]
        for poly in [v.numer,v.denom]:
            ns=[]
            for f in known:
                n=0
                while not poly.rem(f):poly=poly.exquo(f);n+=1
                ns.append(n)
            assert poly.is_ground and poly;ground.append(poly);powers.append(ns)
        return dict(constant=str(K(ground[0])/K(ground[1])),a_F_B_C_powers=[x-y for x,y in zip(*powers)])
    selected=[]
    for i,step in enumerate(c['unit_pivot_steps']):
        assert step['row']==i
        k,j=step['swap_row'],step['cubic_column'];v=parse(step['pivot'])
        assert B[k][j]==v and step['source319_column']==cubic[j]
        assert unit(v)==step['unit']
        B[i],B[k]=B[k],B[i];T[i],T[k]=T[k],T[i]
        B[i]=[x/v for x in B[i]];T[i]=[x/v for x in T[i]]
        el=[[r,str(B[r][j])] for r in range(4) if r!=i and B[r][j]];assert el==step['eliminations']
        for r,vstr in el:
            v=parse(vstr);B[r]=[x-v*y for x,y in zip(B[r],B[i])];T[r]=[x-v*y for x,y in zip(T[r],T[i])]
        selected.append(j)
    assert len(selected)==4
    assert [[str(v) for v in row] for row in T]==c['pivot_row_transform']
    assert [[str(v) for v in row] for row in B]==c['post_pivot_cubic_matrix']
    assert [[sum((T[i][k]*M[k][j] for k in range(4)),K.zero) for j in range(112)] for i in range(4)]==B
    assert [[B[i][j] for j in selected] for i in range(4)]==[[K(int(i==j)) for j in range(4)] for i in range(4)]
    assert [cubic[j] for j in selected]==[176,177,237,181]
    assert [(definitions[cubic[j]][2],definitions[cubic[j]][3]) for j in selected]==[(18,2),(18,3),(86,0),(20,0)]
    composite=[[sum((T[i][k]*U[11+k][j] for k in range(4)),K.zero) for j in range(15)] for i in range(4)]
    assert [[str(v) for v in row] for row in composite]==scope['new4_relation_coefficients_as_combinations_of_parent15']
    assert len(nodes)==10710 and not c['additional_DAG_nodes'] and c['parent_DAG_node_count']==10710
    assert scope['combined_uniform_kernel_directions_on_FULLU']==37 and scope['missing_for_80_gate']==43
    evidence=dict(required_scalar_DAG_nodes_including_zero=len(needed),core_atoms_evaluated=len(acache),
        independent_final4_by15_composite=[[str(v) for v in row] for row in composite],
        selected319_columns=[cubic[j] for j in selected],all132_quartic_zero=True,all448_exact_cubic_match=True,
        nonzero_cubic_entries=119)
    ep=out/'independent-certificate.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',source,note,ap,cp,sp,pcp,pivpath,parentreview,qp]
    result=dict(id=STEM,all_passed=True,
        accepted_parent15_source_identities_and_transport_reused_no_source175_replay=True,
        own_bottom_up_scalar_DAG_evaluator_no_Oracle_or_producer_invocation=True,
        all132_quartic_coefficients_exact_zero_and448_cubic_match=True,nonzero_cubic_coefficients=119,
        saved4_unit_pivot_operations_and_full_matrix_transform_verified=True,
        all4_by15_explicit_composites_and_selected_I4_verified=True,
        no_parent_DAG_mutation=True,new_direct15_independent_everywhere_on_aFBC=True,
        combined_with_old22_only_on_FULLU=37,remaining_for80_gate=43,
        raw_null_field_is_unexecuted_cubic_zero_branch_not_positive37_result=True,
        full80_actualC122_rank64_full239_and_affine_exhaustion_claimed=False,
        no_expanded_quadratic_cofactors_239_matrix_new_point_or_GB=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:result[k] for k in ['all_passed','combined_with_old22_only_on_FULLU','remaining_for80_gate','seconds']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);main(p.parse_args().output)
